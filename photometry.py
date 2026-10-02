import csv
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import ICRS
from astropy import units as u
from astropy.stats import SigmaClip
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry
from astroquery.vizier import Vizier

SAT_THRESHOLD = 62000

target_cats = {
    'B': {'Bmag':  ['I/322A/out', 'II/336/apass9']},
    'R': {'rmag':  ['I/322A/out', 'II/349/ps1'],
          'r_mag': ['II/336/apass9']},
    'I': {'imag':  ['I/322A/out', 'II/349/ps1'],
          'i_mag': ['II/336/apass9']},
}

priority_order = ['I/322A/out', 'II/336/apass9', 'II/349/ps1']

def sort_by_priority(cat_list):
    return sorted(cat_list, key=lambda c: priority_order.index(c) if c in priority_order else len(priority_order))

for band in target_cats:
    for col_name in target_cats[band]:
        target_cats[band][col_name] = sort_by_priority(target_cats[band][col_name])

def load_stack(path):
    with fits.open(path) as hdul:
        image = hdul[0].data
        header = hdul[0].header
    wcs = WCS(header)
    return image, header, wcs

def display_and_select(image, roles=None):
    if roles is None:
        roles = ["target", "ref1", "ref2", "ref3", "nonref1", "nonref2", "nonref3"]
    fig, ax = plt.subplots(figsize=(10, 10))
    vmin, vmax = np.percentile(image, [5, 99])
    ax.imshow(image, cmap='gray', vmin=vmin, vmax=vmax, origin='lower')

    sat_overlay = np.where(image >= SAT_THRESHOLD, 1.0, np.nan)
    ax.imshow(sat_overlay, cmap='Reds', vmin=0, vmax=1, alpha=0.5, origin='lower')

    selections = {}
    for role in roles:
        ax.set_title(f"Click {role}")
        fig.canvas.draw()
        x, y = plt.ginput(1, timeout=0)[0]
        selections[role] = (x, y)
        ax.plot(x, y, 'o', mfc='none', mec='cyan', ms=15)
        fig.canvas.draw()

    plt.close(fig)
    return selections

def query_catalog(ra, dec, filter_band, radius_arcsec=5.0):
    coord = ICRS(ra * u.degree, dec * u.degree)
    for col_name, cat_list in target_cats[filter_band].items():
        for cat in cat_list:
            print(cat)
            try:
                result = Vizier.query_region(coord, radius=radius_arcsec * u.arcsec, catalog=cat)
                if result is not None:
                    time.sleep(0.1)
            except Exception as err:
                print(f"Query failed for {cat}: {err}")
                continue
            if result and len(result) > 0:
                value = float(result[0][col_name][0])
                if not np.isnan(value):
                    print(f"{filter_band}: {cat}")
                    return value, cat
    return None, None

def measure_fwhm(image, x, y, cutout_size=21):
    half = cutout_size // 2
    x_int, y_int = int(round(x)), int(round(y))
    cutout = image[y_int - half:y_int + half + 1, x_int - half:x_int + half + 1]
    cutout_sub = np.maximum(cutout - np.median(cutout), 0)

    yc, xc = np.indices(cutout_sub.shape)
    total = cutout_sub.sum()
    mom1_x = (xc * cutout_sub).sum() / total
    mom1_y = (yc * cutout_sub).sum() / total
    mom2_x = ((xc - mom1_x)**2 * cutout_sub).sum() / total
    mom2_y = ((yc - mom1_y)**2 * cutout_sub).sum() / total

    sigma_x, sigma_y = mom2_x**0.5, mom2_y**0.5
    sigma_to_fwhm = np.sqrt(8 * np.log(2))
    fwhm_x, fwhm_y = sigma_x * sigma_to_fwhm, sigma_y * sigma_to_fwhm
    fwhm = (fwhm_x * fwhm_y)**0.5

    x_centroid = (x_int - half) + mom1_x
    y_centroid = (y_int - half) + mom1_y

    return fwhm, x_centroid, y_centroid

def measure_aperture(image, noise, x, y, fwhm):
    pos = (x, y)
    src_ap = CircularAperture(pos, r=2 * fwhm)
    sky_an = CircularAnnulus(pos, r_in=3 * fwhm, r_out=5 * fwhm)

    sky_stats = ApertureStats(image, sky_an, sigma_clip=SigmaClip(sigma=3.0))
    sky_per_pix = sky_stats.median
    sky_std = sky_stats.std

    phot = aperture_photometry(image, src_ap, error=noise)
    raw_sum = phot['aperture_sum'][0]
    raw_err = phot['aperture_sum_err'][0]

    ap_area = src_ap.area
    net_flux = raw_sum - sky_per_pix * ap_area
    net_err = np.sqrt(raw_err**2 + ap_area * sky_std**2)

    return net_flux, net_err, sky_per_pix, sky_std

def compute_zp(refs, exptime):
    zps = [cat_mag + 2.5 * np.log10(flux / exptime) for flux, cat_mag in refs]
    zp_mean = np.mean(zps)
    zp_std = np.std(zps, ddof=1)
    return zp_mean, zp_std

def calibrate(flux, flux_err, exptime, zp_mean, zp_std):
    instr_mag = -2.5 * np.log10(flux / exptime)
    instr_err = (2.5 / np.log(10)) * (flux_err / flux)
    calib_mag = instr_mag + zp_mean
    calib_err = np.sqrt(instr_err**2 + zp_std**2)
    return calib_mag, calib_err, instr_mag

def photometry(stack_path, noise_path, exptime, filter_band):
    image, _, wcs = load_stack(stack_path)
    noise = fits.getdata(noise_path)

    print(f"\n=== {stack_path} ===")
    print(f"Filter: {filter_band}, EXPTIME: {exptime}s")

    selections = display_and_select(image)
    measurements = {}

    for role, (x, y) in selections.items():
        ra, dec = wcs.pixel_to_world_values(x, y)
        if role == 'target':
            cat_mag, cat_id = None, None
        else:
            cat_mag, cat_id = query_catalog(float(ra), float(dec), filter_band)
            while cat_mag is None:
                print(f"No catalog match for {role}. Re-click.")
                new_sel = display_and_select(image, roles=[role])
                x, y = new_sel[role]
                ra, dec = wcs.pixel_to_world_values(x, y)
                cat_mag, cat_id = query_catalog(float(ra), float(dec), filter_band)
        measurements[role] = {'x': x, 'y': y, 'ra': float(ra), 'dec': float(dec),
                              'cat_mag': cat_mag, 'cat_id': cat_id}

    for role, star in measurements.items():
        fwhm, x_centroid, y_centroid = measure_fwhm(image, star['x'], star['y'])
        star.update({'fwhm': fwhm, 'x_centroid': x_centroid, 'y_centroid': y_centroid})

    fwhms = [star['fwhm'] for star in measurements.values()]
    fwhm_mean = float(np.mean(fwhms))
    fwhm_std = float(np.std(fwhms, ddof=1))

    for role, star in measurements.items():
        flux, flux_err, sky, sky_std = measure_aperture(image, noise, star['x_centroid'], star['y_centroid'], fwhm_mean)
        x_int, y_int = int(round(star['x_centroid'])), int(round(star['y_centroid']))
        peak = float(image[max(y_int - 2, 0):y_int + 3, max(x_int - 2, 0):x_int + 3].max())
        star.update({'flux': flux, 'flux_err': flux_err, 'sky': sky,
                     'sky_std': sky_std, 'peak': peak, 'ap_radius': 2 * fwhm_mean})

    refs = [(measurements[role]['flux'], measurements[role]['cat_mag']) for role in ('ref1', 'ref2', 'ref3')]
    zp_mean, zp_std = compute_zp(refs, exptime)

    for role, star in measurements.items():
        calibrated_mag, calibrated_err, instrumental_mag = calibrate(star['flux'], star['flux_err'], exptime, zp_mean, zp_std)
        star.update({'instr_mag': instrumental_mag, 'calib_mag': calibrated_mag, 'calib_err': calibrated_err})

    print(f"Mean FWHM: {fwhm_mean:.2f} ± {fwhm_std:.2f} pix")
    print(f"ZP: {zp_mean:.3f} ± {zp_std:.3f}")
    print(f"Target calib mag: {measurements['target']['calib_mag']:.3f} ± {measurements['target']['calib_err']:.3f}")
    print("Sanity check (nonrefs):")
    for role in ('nonref1', 'nonref2', 'nonref3'):
        star = measurements[role]
        print(f"  {role}: derived={star['calib_mag']:.3f} cat={star['cat_mag']:.3f} Δ={star['calib_mag']-star['cat_mag']:+.3f}")

    cols = ['role', 'x', 'y', 'ra', 'dec', 'fwhm', 'peak', 'ap_radius',
            'flux', 'flux_err', 'sky', 'sky_std',
            'cat_mag', 'cat_id', 'instr_mag', 'calib_mag', 'calib_err']
    csv_path = stack_path.replace('.fits', '_photometry.csv')
    with open(csv_path, 'w', newline='') as csv_file:
        csv_file.write(f"# stack={stack_path}\n")
        csv_file.write(f"# filter={filter_band}\n")
        csv_file.write(f"# exptime={exptime}\n")
        csv_file.write(f"# fwhm_mean={fwhm_mean:.4f}\n")
        csv_file.write(f"# fwhm_std={fwhm_std:.4f}\n")
        csv_file.write(f"# zp_mean={zp_mean:.4f}\n")
        csv_file.write(f"# zp_std={zp_std:.4f}\n")
        writer = csv.writer(csv_file)
        writer.writerow(cols)
        for role, star in measurements.items():
            writer.writerow([role] + [star.get(col) for col in cols[1:]])
    print(f"Wrote {csv_path}")

    return measurements, zp_mean, zp_std


stacks = [
    ('J0813_2026-04-14/J08_14_B_solved.fits',
     'J0813_2026-04-14/J0813_2026-04-14_B_300.0_noise.fits',
     300.0, 'B'),
    ('J0813_2026-04-14/J08_14_R_solved.fits',
     'J0813_2026-04-14/J0813_2026-04-14_R_300.0_noise.fits',
     300.0, 'R'),
    ('J0813_2026-04-14/J08_14_I_solved.fits',
     'J0813_2026-04-14/J0813_2026-04-14_I_300.0_noise.fits',
     300.0, 'I'),
    ('J0813_2026-04-23/J08_23_R_solved.fits',
     'J0813_2026-04-23/J0813_2026-04-23_R_180.0_noise.fits',
     180.0, 'R'),
    ('J0813_2026-04-23/J08_23_I_solved.fits',
     'J0813_2026-04-23/J0813_2026-04-23_I_300.0_noise.fits',
     300.0, 'I'),
]

for stack_path, noise_path, exptime, filter_band in stacks:
    photometry(stack_path, noise_path, exptime, filter_band)
