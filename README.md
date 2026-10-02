# CCD Photometry Pipeline — Quasar UV Variability

A calibration, photometry, and SED-fitting pipeline designed for small-telescope CCD data,
together with the dataset and results it produced for a two-epoch study of the
quasar SDSS J081331.28+254503.0 (z = 1.512).

See the writeup:
**[Constraining Quasar UV Variability for the L_X–L_UV Distance Ladder](https://jacksonlyle.dev/papers/quasar-uv-variability-2026.pdf)**

## The project

At z = 1.512 the observed SDSS r band maps to a rest-frame wavelength near 2500 Å —
the continuum point used in the Risaliti & Lusso L_X–L_UV relation, which is being
explored as a standard candle in the cosmological distance ladder. The dispersion in that relation is
sensitive to UV variability, so multi-epoch r-band photometry of a single quasar
probes the wavelength that matters.

The observing proposal was written for AST3722C (Techniques of Observational
Astronomy) and was ~4 of 40 selected for telescope time through the class proposal rounds. The
target was observed on two nights nine days apart with the 14" Meade LX200GPS at
Rosemary Hill Observatory, in SDSS r, SDSS i, and Johnson B. This repository holds
the code that reduced those frames, the frames themselves, and the measured results.

The measured variability was consistent with zero at 0.88σ — a null result, but with
photometric precision sufficient to detect quasar variability at this magnitude from
a small-aperture telescope. Ideally, more epochs of observations would allow for much better
contraints on the variability.

## Pipeline

| Script | What it does |
|---|---|
| `data_reduction.py` | Builds master bias, darks (per exposure time), and flats (per filter); calibrates and aligns light frames; stacks them; propagates per-pixel noise into a matching noise map. |
| `photometry.py` | Interactive aperture photometry: click the target and reference stars, measure FWHM, query reference magnitudes from Vizier (UCAC4 → APASS9 → PS1 by priority), compute a zero point, and output calibrated magnitudes with uncertainties. |
| `sed_fitting.py` | Fits a power law across the B/r/i photometry, interpolates the rest-frame 2500 Å flux, and compares the two epochs. |

Noise is propagated through every step rather than estimated at the end, so the final
magnitude uncertainties trace back to read noise, dark current, and flat-field error.

## Running it

You can use any environment manager you like, but I will assume uv (my favorite) for this demo.

The data reduction pipeline utilizes an interactive selector for your target and calibration stars,
It will attempt to automatically select the same stars for every frame, so you should only have to 
use the selector once per filter you use. For the sake of this demo, all of the necessary photometric
data is already hardcoded in sed_fitting.py, so running the reduction and photometry pipelines are
not necessary to reproduce the results.

```bash
uv sync
python data_reduction.py     # produces masters/ and the stacked frames
# plate-solve the stacks (see below)
python photometry.py         # interactive: click target, then reference stars
# enter photometry values in sed_fitting.py
python sed_fitting.py        # power-law fit and 2500 Å interpolation
```

### Plate solving

`photometry.py` needs WCS coordinates in the headers so it can query catalogs, and
`data_reduction.py` does not add them (image solving was out-of-scope). The stacks were solved with
[astrometry.net](https://nova.astrometry.net):

1. Upload a `*_stack.fits` file from the data reduction pipeline (solved frames are already included in this repo)
2. Set a scale hint — these frames are 2.82″/pixel (0.94″/pixel native, 3×3 binned).
3. Download the "new-image.fits" result, which is your frame with WCS keywords added.
4. Save it in the same dir as the stack.fits frames

The solved frames for all five stacks are already included, so `photometry.py` runs
without this step.

### Entering photometry data

`photometry.py` returns csv files containing extensive statistics on the sampled frames. To correctly 
run `sed_fitting.py` for different data, you will need to change the hardcoded `magnitudes` dict 
to reflect your own values from the `photometry.py` outputs. Note that the keys '4_14' and '4_23'
are used throughout this pipeline, so running the sed fitting pipeline will require modifications 
if you choose to use different or more epochs of observations.

As is, the pipeline should run with the actual data from the paper without any modifications.

## Running with different data

Each script has its configuration in a block at the top: the parent dir (here, using YYYY-MM-DD), the object name,
and the glob patterns for each frame type:

```python
night = "2026-04-23"
obj_name = f"J0813_{night}"
dark_paths = glob.glob(f"{night}/Darks/*.fits")
```
Running the scripts with your own data requires changing these configs to point to your intended frames and directories.

## What's in the repository

```
2026-04-23/              Raw sample frames (SDSS r path only — see note)
  Biases/  Darks/  Flats/  Lights/
J0813_2026-04-14/        Plate-solved stacks + noise maps (B, r, i)
J0813_2026-04-23/        Plate-solved stacks + noise maps (r, i)
figures/                 Result figures from the paper
```

This repository only gives a subset of the actual data with the SDSS r frames (mainly for file size)
along with the necessary calibration frames. The actual results can be found in the paper above.

## Analysis

The goal of this lab was to see if SED fitting could be used to detect and constrain the variability
of quasars. This is important for the characterization of their brightness, allowing for further 
study of their potential as a standard candle. Ideally, this means reproducing a statistically  significant
difference in 2500 angstrom flux values on a cadence of ~4 weeks (depending on the target quasar).

This report proved the concept to be viable for amateur telescopes to provide meaningful photometric
contraints on variability for objects at or even above redshifts of z = 1.5, however, due to specific
contraints with time and lunar phase this project on its own is not enough to more generally contrain 
the variability. Further research with more complete data along with this pipeline could be critical 
for furthering the characterization of these objects.

## License

MIT — see [LICENSE](LICENSE).
