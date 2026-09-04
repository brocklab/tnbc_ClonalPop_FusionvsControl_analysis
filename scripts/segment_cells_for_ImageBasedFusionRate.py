## Use environment_fusionrate_segment.yml

import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cellpose import models
from skimage import measure
from scipy.spatial import KDTree
import time
import re

# Limit threads to 20 (shared server)
os.environ["OMP_NUM_THREADS"] = "12"
os.environ["OPENBLAS_NUM_THREADS"] = "12"
os.environ["MKL_NUM_THREADS"] = "12"
os.environ["VECLIB_MAXIMUM_THREADS"] = "12"
os.environ["NUMEXPR_NUM_THREADS"] = "12"

# === CONFIGURATION ===
input_dir = '/stor/work/Brock/kennedy/SC_repo/data/ImageBasedFusionRateAnalysis/representative_MDA-MB-231_images'
output_dir = '/stor/work/Brock/kennedy/SC_repo/scripts/outputs'
######## Note: change experiment header for filename under def extract_metadata #########

# === FUNCTIONS ===
def remove_scale_bar(image, percentage=0.1):
    h = image.shape[0]
    image[int(h * (1 - percentage)):, :] = 0
    return image

def segment_nuclei(image, model):
    gray = image if len(image.shape) == 2 else cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    enhanced = clahe.apply(gray)
    masks, *_ = model.eval(enhanced, diameter=None, channels=[0, 0],
                           flow_threshold=1.0, cellprob_threshold=-6, augment=False)
    return masks

def calculate_iou(mask1, mask2):
    intersection = np.logical_and(mask1, mask2).sum()
    union = np.logical_or(mask1, mask2).sum()
    return intersection / union if union > 0 else 0

def region_intensity(region, img):
    return np.mean([img[coord[0], coord[1]] for coord in region.coords])

def extract_metadata(filename):
    pattern = r"KH2514_(DOX|VIN|PAX|CTRL)_(GFP|MC|overlap)_([A-Z]+\d+)_\d+_(\d{2})d\d{2}h\d{2}m"
    match = re.search(pattern, filename)
    if match:
        return match.groups()
    return None

# === LOAD AND GROUP FILES ===
grouped = {}
for fname in os.listdir(input_dir):
    meta = extract_metadata(fname)
    if not meta:
        continue
    drug, channel, well, day = meta
    image_key = f"{drug}_{well}_{day}"
    grouped.setdefault(image_key, {})[channel.lower()] = os.path.join(input_dir, fname)

# === START MODEL ===
model = models.Cellpose(gpu=False, model_type='nuclei')

skipped_keys = []

# === MAIN LOOP ===
for image_key, paths in grouped.items():
    if not all(k in paths for k in ('gfp', 'mc', 'overlap')):
        skipped_keys.append(image_key)
        continue

    drug, well, day = image_key.split('_')
    well_dirname = f"{drug}_{well}"
    well_dir = os.path.join(output_dir, well_dirname)
    os.makedirs(well_dir, exist_ok=True)

    out_csv = os.path.join(well_dir, f'{image_key}_classified.csv')
    if os.path.exists(out_csv):
        continue

    print(f"\nProcessing {image_key} at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    start_time = time.time()

    try:
        # === Load and preprocess ===
        gfp = remove_scale_bar(cv2.imread(paths['gfp'], cv2.IMREAD_GRAYSCALE))
        mc = remove_scale_bar(cv2.imread(paths['mc'], cv2.IMREAD_GRAYSCALE))
        overlap_rgb = remove_scale_bar(cv2.cvtColor(cv2.imread(paths['overlap'], cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB))

        masks_gfp = segment_nuclei(gfp, model)
        masks_mc = segment_nuclei(mc, model)

        regions_gfp = measure.regionprops(masks_gfp)
        regions_mc = measure.regionprops(masks_mc)

        # === Classify ===
        mc_centroids = np.array([r.centroid for r in regions_mc])
        mc_region_map = {r.label: r for r in regions_mc}
        matched_mc_labels = set()
        classified = []

        if len(mc_centroids) > 0:
            mc_tree = KDTree(mc_centroids)

        for gfp_region in regions_gfp:
            gfp_mask = np.zeros_like(gfp, dtype=bool)
            gfp_mask[tuple(gfp_region.coords.T)] = True
            gfp_centroid = gfp_region.centroid
            matched = False
            gfp_int = region_intensity(gfp_region, gfp)
            mc_int = region_intensity(gfp_region, mc)

            if len(mc_centroids) > 0:
                nearby_idx = mc_tree.query_ball_point(gfp_centroid, r=10)
                for idx in nearby_idx:
                    mc_region = mc_region_map[regions_mc[idx].label]
                    if mc_region.label in matched_mc_labels:
                        continue
                    mc_mask = np.zeros_like(mc, dtype=bool)
                    mc_mask[tuple(mc_region.coords.T)] = True
                    if calculate_iou(gfp_mask, mc_mask) > 0.8:
                        matched_mc_labels.add(mc_region.label)
                        classified.append({'x': gfp_centroid[1], 'y': gfp_centroid[0],
                                           'label': 'fusion',
                                           'gfp_intensity': gfp_int,
                                           'mcherry_intensity': mc_int})
                        matched = True
                        break

            if not matched:
                classified.append({'x': gfp_centroid[1], 'y': gfp_centroid[0],
                                   'label': 'gfp',
                                   'gfp_intensity': gfp_int,
                                   'mcherry_intensity': mc_int})

        for mc_region in regions_mc:
            if mc_region.label not in matched_mc_labels:
                mc_centroid = mc_region.centroid
                gfp_int = region_intensity(mc_region, gfp)
                mc_int = region_intensity(mc_region, mc)
                classified.append({'x': mc_centroid[1], 'y': mc_centroid[0],
                                   'label': 'mcherry',
                                   'gfp_intensity': gfp_int,
                                   'mcherry_intensity': mc_int})

        df = pd.DataFrame(classified)
        df.to_csv(out_csv, index=False)

        # === Save images ===
        def save_img(df_sub, fname, title, color='yellow'):
            plt.figure(figsize=(12, 12))
            plt.imshow(overlap_rgb)
            for _, row in df_sub.iterrows():
                rect = plt.Rectangle((row['x'] - 5, row['y'] - 5), 10, 10,
                                     linewidth=0.8, edgecolor=color, facecolor='none')
                plt.gca().add_patch(rect)
            plt.title(title)
            plt.axis('off')
            plt.savefig(os.path.join(well_dir, fname), dpi=300)
            plt.close()

        save_img(df, f'{image_key}_all.png', "All Classified Nuclei", color='white')
        save_img(df[df['label'] == 'fusion'], f'{image_key}_fusion.png', "Fusion Only", color='yellow')

        # === Update summary CSV ===
        summary_path = os.path.join(well_dir, f"{drug}_{well}_summary.csv")
        counts = {
            'day': int(day),
            'GFP': (df['label'] == 'gfp').sum(),
            'MC': (df['label'] == 'mcherry').sum(),
            'FUSION': (df['label'] == 'fusion').sum()
        }

        if os.path.exists(summary_path):
            df_summary = pd.read_csv(summary_path)
            df_summary = df_summary[df_summary['day'] != counts['day']]
            df_summary = pd.concat([df_summary, pd.DataFrame([counts])], ignore_index=True)
        else:
            df_summary = pd.DataFrame([counts])

        df_summary = df_summary.sort_values(by='day')
        df_summary.to_csv(summary_path, index=False)

        elapsed = time.time() - start_time
        print(f"Finished {image_key} in {elapsed:.2f} sec")

    except Exception as e:
        print(f"Error processing {image_key}: {e}")
        skipped_keys.append(image_key)
        continue

# === Report Skipped ===
print("\n=== Script complete ===")
if skipped_keys:
    print("Skipped image keys:")
    for k in skipped_keys:
        print(" -", k)
else:
    print("All image keys processed without error.")