import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astroalign as aa
import shutil
import os
import glob

night = "2026-04-23"
obj_name = f"J0813_{night}"
dark_paths = glob.glob(f"{night}/Darks/*.fits")
bias_paths = glob.glob(f"{night}/Biases/*.fits")
flat_paths = glob.glob(f"{night}/Flats/*.fits")
light_paths = glob.glob(f"{night}/Lights/*.fits")

def create_master_darks(dark_paths):
    darks = {}
    master_darks = {}
    dark_noise = {}

    for frame in dark_paths:
        header = fits.getheader(frame)
        exp = header['EXPTIME']
        if exp not in darks:
            darks[exp] = []
        darks[exp].append(frame)

    for exp, paths in darks.items():
        stack = np.array([fits.getdata(path) for path in paths])
        master_darks[exp] = np.median(stack, axis=0)
        dark_noise[exp] = np.std(stack, axis=0) / np.sqrt(len(stack))

    return master_darks, dark_noise

def create_master_bias(bias_paths):
    stack = np.array([fits.getdata(path) for path in bias_paths])
    master_bias = np.median(stack, axis=0)
    bias_noise = np.std(stack, axis=0) / np.sqrt(len(stack))
    return master_bias, bias_noise

def create_master_flats(flat_paths, master_darks, dark_noise, bias_paths):
    flats = {}
    master_flats = {}
    flat_noise = {}
    total_flat_noise = {}

    master_bias, bias_noise = create_master_bias(bias_paths)

    for frame in flat_paths:
        header = fits.getheader(frame)
        filter = header['FILTER']
        if filter not in flats:
            flats[filter] = []
            flat_noise[filter] = []
        flats[filter].append(frame)

    for filter, paths in flats.items():
        stack = []
        for path in paths:
            frame = fits.getdata(path)
            header = fits.getheader(path)
            exp = header["EXPTIME"]

            if exp in master_darks:
                calib_frame = master_darks[exp]
                calib_noise = dark_noise[exp]
            else:
                #If no exact darks are found for an exposure time, use the closest possible dark
                dark_exp = min(master_darks.keys(), key=lambda x: abs(x - exp))

                dark_current = (master_darks[dark_exp] - master_bias) / dark_exp
                dark_current_noise = np.sqrt(dark_noise[dark_exp]**2 + bias_noise**2) / dark_exp 

                calib_frame = master_bias + dark_current * exp
                calib_noise = np.sqrt(bias_noise**2 + (exp * dark_current_noise)**2)

            flat_noise[filter].append(np.sqrt(frame + calib_noise**2))
            stack.append(frame - calib_frame)

        noise_stack = np.array(flat_noise[filter])
        total_flat_noise[filter] = np.sqrt(np.sum(noise_stack**2, axis=0)) / len(noise_stack)

        master_flat = np.mean(np.array(stack), axis=0)
        master_flats[filter] = master_flat / np.mean(master_flat)

        total_flat_noise[filter] = total_flat_noise[filter] / np.mean(master_flat)

    return master_flats, total_flat_noise

def create_master_lights(light_paths, master_darks, dark_noise, master_flats, flat_noise, master_bias, bias_noise, obj_name):

    if os.path.exists(obj_name):
        shutil.rmtree(obj_name)
    os.mkdir(obj_name)
    lights = {}
    light_noise = {}
    calib_lights = {}
    
    #sort frames by filter and exposure time
    for path in light_paths:
        header = fits.getheader(path)
        filter = header['FILTER']
        exp = header['EXPTIME']

        if filter not in lights:
            lights[filter] = {}
            calib_lights[filter] = {}
            light_noise[filter] = {}

        if exp not in lights[filter]:
            lights[filter][exp] = []
            calib_lights[filter][exp] = []
            light_noise[filter][exp] = []

        lights[filter][exp].append(path)

    for filter, data in lights.items():
        for exp, paths in data.items():
            for path in paths:
                frame = fits.getdata(path)

                if exp in master_darks:
                    calib_frame = master_darks[exp]
                    calib_noise = dark_noise[exp]
                else:
                    #If no exact darks are found for an exposure time, use the closest possible dark
                    dark_exp = min(master_darks.keys(), key=lambda x: abs(x - exp))

                    dark_current = (master_darks[dark_exp] - master_bias) / dark_exp
                    dark_current_noise = np.sqrt(dark_noise[dark_exp]**2 + bias_noise**2) / dark_exp

                    calib_frame = master_bias + dark_current * exp
                    calib_noise = np.sqrt(bias_noise**2 + (exp * dark_current_noise)**2)

                dark_sub = frame - calib_frame
                dark_sub_noise = np.sqrt(frame + calib_noise**2)
                calib_light = dark_sub / master_flats[filter]
                calib_lights[filter][exp].append(calib_light)

                shot_est = np.median((np.sqrt(frame) / master_flats[filter])[dark_sub > 0])
                dark_est = np.median((calib_noise / master_flats[filter])[dark_sub > 0])
                flat_est = np.median((np.abs(calib_light) * flat_noise[filter] / master_flats[filter])[dark_sub > 0])

                print(f"{filter} exp={exp}: shot={shot_est:.3f}, dark={dark_est:.3f}, flat={flat_est:.3f}")

                frac_dark = np.zeros_like(dark_sub, dtype=float)
                frac_flat = np.zeros_like(master_flats[filter], dtype=float)

                mask_dark = dark_sub > 0
                mask_flat = master_flats[filter] > 0

                frac_dark[mask_dark] = (dark_sub_noise[mask_dark] / dark_sub[mask_dark])**2
                frac_flat[mask_flat] = (flat_noise[filter][mask_flat] / master_flats[filter][mask_flat])**2

                calib_light_noise = np.abs(calib_light) * np.sqrt(frac_dark + frac_flat)
                light_noise[filter][exp].append(calib_light_noise)

    #Due to the large number of files, I will be using the astroalign package for alignment
    aligned_lights = {}
    aligned_noise = {}
    master_lights = {}
    total_light_noise = {}

    for filter, data in calib_lights.items():
        aligned_lights[filter] = {}
        aligned_noise[filter] = {}
        
        #here I am choosing the first frame in each filter as the reference
        #only really matters if there are multiple exposure times for the same filter, which there shouldn't be.
        ref_exp = list(data.keys())[0]
        reference = data[ref_exp][0]

        for exp, frames in data.items():
            aligned_lights[filter][exp] = []
            aligned_noise[filter][exp] = []

            for frame_num, frame in enumerate(frames):
                noise = light_noise[filter][exp][frame_num]

                if exp == ref_exp and frame_num == 0:
                    aligned_lights[filter][exp].append(frame)
                    aligned_noise[filter][exp].append(noise)
                else:
                    try:
                        transform, _ = aa.find_transform(frame, reference)
                        aligned_frame, footprint = aa.apply_transform(transform, frame, reference)
                        aligned_noise_frame, _ = aa.apply_transform(transform, noise, reference)

                        aligned_lights[filter][exp].append(aligned_frame)
                        aligned_noise[filter][exp].append(aligned_noise_frame)
                    except:
                        print(f"Failed to find the transform for: {filter}, {exp}, {frame_num}.")

        all_frames = []
        all_noise = []

        #Here, I am combining all of the different exposure times per filter into the final stacked frame
        #Again, it shouldn't matter if the exposure times are all correct, but since I've allowed for
        #differing exposure times in my pipeline, I will stay self-consistent

        for exp in aligned_lights[filter]:
            all_frames.extend(aligned_lights[filter][exp])
            all_noise.extend(aligned_noise[filter][exp])

        frame_stack = np.array(all_frames)
        noise_stack = np.array(all_noise)

        master_lights[filter] = np.median(frame_stack, axis=0)
        total_light_noise[filter] = np.sqrt(np.sum(noise_stack**2, axis=0)) / len(noise_stack)
        output_filename = f"{obj_name}/{obj_name}_{filter}_{ref_exp}_stack.fits"
        noise_filename = f"{obj_name}/{obj_name}_{filter}_{ref_exp}_noise.fits"
        fits.writeto(output_filename, master_lights[filter], overwrite=True)
        fits.writeto(noise_filename, total_light_noise[filter], overwrite=True)
    return master_lights, total_light_noise


masters_dir = f"masters/{obj_name}"
os.makedirs(masters_dir, exist_ok=True)

master_darks, dark_noise = create_master_darks(dark_paths)
master_bias, bias_noise   = create_master_bias(bias_paths)
master_flats, flat_noise  = create_master_flats(flat_paths, master_darks, dark_noise, bias_paths)
master_lights, total_light_noise = create_master_lights(light_paths, master_darks, dark_noise, master_flats, flat_noise, master_bias, bias_noise, obj_name)

fits.writeto(f"{masters_dir}/master_bias.fits", master_bias, overwrite=True)
for exp, dark in master_darks.items():
    fits.writeto(f"{masters_dir}/master_dark_{exp}s.fits", dark, overwrite=True)
for filt, flat in master_flats.items():
    fits.writeto(f"{masters_dir}/master_flat_{filt}.fits", flat, overwrite=True)

read_noise = np.std(master_bias)
print(f"Calibration Stats: {obj_name}")
print(f"Read noise (std of master bias): {read_noise:.4f} counts/pixel")
for exp in sorted(master_darks.keys()):
    dark_current = np.mean((master_darks[exp] - master_bias) / exp)
    print(f"Dark current ({exp:6.3f}s dark):  {dark_current:.6f} counts/pixel/sec")
