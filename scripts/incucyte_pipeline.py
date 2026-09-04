# incucyte_pipeline.py
# Engine for Incucyte morphology analysis (Cells & Nuclei)
# Filename format: <root>_<chan>_<sample>_<well>_<site>_<time>.tif
# Example: KH2506_GFP_F7_B6_1_01d00h00m.tif

import os, re, math, warnings
from dataclasses import dataclass
from typing import Optional, Dict, List
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from skimage import filters, morphology, measure, segmentation, util, feature, color, exposure, restoration
from skimage.morphology import remove_small_objects, remove_small_holes, white_tophat, disk
from scipy import ndimage as ndi
import json

# ---- Optional: Cellpose for best cell segmentation on phase ----
try:
    from cellpose import models
    _HAS_CELLPOSE = True
except Exception:
    _HAS_CELLPOSE = False
    warnings.warn("Cellpose not found. Classic segmentation will be used for cells.", RuntimeWarning)

# ---- Safe thread defaults (you can override via env) ----
os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("MKL_NUM_THREADS", "6")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "6")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "6")


# =========================
# Config dataclass
# =========================
@dataclass
class PipelineConfig:
    parent_dir: str = "./"
    outputs_dir: str = "./outputs"

    # <root>_<chan>_<sample>_<well>_<site>_<time>.tif
    filename_regex: str = r"^(?P<root>.+?)_(?P<chan>GFP|MC|PHASE|OVERLAP)_(?P<sample>[A-Za-z0-9-]+)_(?P<well>[A-Za-z0-9-]+)_(?P<site>[A-Za-z0-9]+)_(?P<time>\d{2}d\d{2}h\d{2}m)(?:[^.]*)\.(?:tif|tiff|png)$"

    # --- Inner-circle crop (exclude well edge) ---
    use_inner_circle_crop: bool = True     # turn on to enable
    crop_radius_frac: float = 0.70          # 0<r<=1, fraction of min(H,W)/2
    crop_center_xy: Optional[tuple] = None  # (cx, cy) in pixels; None => image center
    hard_mask_images: bool = True           # zero pixels outside circle before seg

    # Scale
    microns_per_pixel: Optional[float] = None
    scale_bar_microns: float = 800.0
    scale_bar_search_height_pct: float = 0.18
    overlap_image_hint: Optional[str] = None

    # Selection & reproducibility
    random_seed: int = 1337
    target_cells: int = 50
    grid_rows: int = 5
    grid_cols: int = 5
    exclude_phase_edge_touching: bool = True
    exclude_nucleus_edge_touching: bool = True

    # Size filters (µm²)
    min_cell_area_um2: float = 50.0
    max_cell_area_um2: float = 5_000.0
    min_nuc_area_um2: float = 15.0
    max_nuc_area_um2: float = 2_000.0

    # FOV scan limit
    max_fovs_per_sample: int = 20

    # Cellpose (CP4)
    use_cellpose: bool = True
    cellpose_model: str = "cyto2"
    cellpose_diameter: Optional[float] = None
    cellpose_flow_threshold: float = 0.4
    cellpose_cellprob_threshold: float = 0.0
    cellpose_batch_size: int = 8
    normalize_cellpose: bool = True
    invert_cellpose: bool = False

    # Classic fallback (phase)
    phase_gaussian_sigma: float = 1.0
    phase_thresh_method: str = "otsu"  # "otsu"|"yen"|"tri"
    phase_min_size_px: int = 200

    # ---- Nuclei (basic) ----
    nuclei_gaussian_sigma: float = 1.0
    nuclei_min_size_px: int = 50
    nuclei_watershed: bool = True
    nuclei_watershed_compactness: float = 0.0

    # ---- Nuclei (advanced) ----
    nuc_top_hat_radius: int = 9        # 5–11 is typical
    nuc_clahe_clip: float = 0.005       # 0.005–0.02
    nuc_log_sigma_min: float = 1.2     # px
    nuc_log_sigma_max: float = 4.2     # px
    nuc_log_threshold: float = 0.018    # LoG response threshold
    refine_nuclei_within_cell: bool = True  # refine inside ROI for the selected cell
    nuc_threshold_offset: float = 0.0      # + makes masks smaller/tighter
    nuc_sauvola_window: int = 25
    nuc_sauvola_k: float = 0.2
    nuc_min_pixels_after_thresh: int = 50
    nuc_shrink_px: int = 0   # 0 = off; 1 = ~1-px shrink before watershed
    # size/cleanup
    nuc_fill_holes_area: int = 128        # fill small interior holes
    nuc_open_radius: int = 0              # 0 disables; 1 is a gentle trim
    # seeding / splitting
    nuc_min_distance_px: int = 8          # fallback peak spacing (prevents merges)
    # Mask smoothing
    mask_smooth_radius_px: int = 0         # 0 = off; 1–2 are gentle, >2 gets blobby
    subpixel_smooth_sigma: float = 1.0   # set 0 to disable
    smooth_for_display_only: bool = False  # <<< set False to affect measurements

    # Viz
    qc_panel_crop_pad_px: int = 20
    qc_panels_per_row: int = 5

    # --- QC saving controls ---
    save_individual_qc_panels: bool = False   # per-cell PNGs
    save_combined_qc_panel: bool = True      # big gallery PNG
    save_selection_maps: bool = True         # NEW: per-FOV "selection map" PNGs

    # --- Selection manifest (replay your reviewed picks) ---
    save_selection_manifest: bool = True
    use_selection_manifest_path: Optional[str] = None  # path to JSON to replay selection (optional)

    # Selection map styling
    selection_map_box_color: tuple = (1.0, 0.0, 0.0)  # RGB red
    selection_map_box_thickness: int = 2
    selection_map_pad_px: int = 6                    # pad around each bbox
    selection_map_label_boxes: bool = True           # show 1..50 on each box

    # Runtime
    save_intermediate_masks: bool = True
    verbose: bool = True


# =========================
# Parsing & IO
# =========================
@dataclass(frozen=True)
class ImageKey:
    root: str
    sample: str
    well: str
    site: str
    timecode: str

def _fname_re(cfg: PipelineConfig):
    return re.compile(cfg.filename_regex)

def parse_filename(p: Path, cfg: PipelineConfig):
    m = _fname_re(cfg).match(p.name)
    if not m:
        return None
    d = m.groupdict()
    return {
        "root": d["root"], "chan": d["chan"], "sample": d["sample"], "well": d["well"],
        "site": d["site"], "timecode": d["time"], "path": str(p)
    }

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def group_fovs_by_key(parent: Path, cfg: PipelineConfig, sample_filter: Optional[str]=None, max_fovs: Optional[int]=None):
    """Return dict ImageKey -> {chan->path}, requiring PHASE, GFP, MC. OVERLAP optional."""
    files = sorted([x for x in parent.iterdir() if x.is_file() and x.suffix.lower() in (".tif",".tiff",".png")])
    records = []
    for f in files:
        info = parse_filename(f, cfg)
        if not info:
            continue
        if sample_filter and info["sample"] != sample_filter:
            continue
        records.append(info)
    fovs = {}
    for r in records:
        key = ImageKey(root=r["root"], sample=r["sample"], well=r["well"], site=r["site"], timecode=r["timecode"])
        fovs.setdefault(key, {})
        fovs[key][r["chan"]] = r["path"]
    need = {"PHASE","GFP","MC"}  # OVERLAP optional
    kept = {k:v for k,v in fovs.items() if need.issubset(set(v.keys()))}
    if max_fovs is not None and len(kept) > max_fovs:
        keys = list(kept.keys())
        keys.sort(key=lambda x: (x.well, x.site, x.timecode))
        rng = np.random.RandomState(cfg.random_seed)
        sel_idx = rng.choice(len(keys), size=max_fovs, replace=False)
        kept = {keys[i]: kept[keys[i]] for i in sel_idx}
    return kept

def load_image(path: str):
    return tiff.imread(path)

def rescale_for_display(img):
    p0, p99 = np.percentile(img, [0.5, 99.5])
    if p99 <= p0:
        return np.zeros_like(img, dtype=np.float32)
    out = (img - p0) / (p99 - p0)
    return np.clip(out, 0, 1).astype(np.float32)

# ========================
# CROP OUTER CIRCLE
# ========================

def _circle_mask(shape, radius_frac=0.9, center_xy=None):
    """Return boolean mask with a filled circle (True = keep)."""
    H, W = shape[:2]
    cy = center_xy[1] if (center_xy is not None) else H/2.0
    cx = center_xy[0] if (center_xy is not None) else W/2.0
    R = radius_frac * 0.5 * min(H, W)
    yy, xx = np.ogrid[:H, :W]
    dist2 = (yy - cy)**2 + (xx - cx)**2
    return dist2 <= (R*R)

def _apply_roi_mask_to_images(mask, phase, g, r):
    """Zero pixels outside ROI (if requested)."""
    phase2 = phase.copy()
    g2 = g.copy()
    r2 = r.copy()
    phase2[~mask] = 0
    g2[~mask] = 0
    r2[~mask] = 0
    return phase2, g2, r2

def _keep_labels_inside_mask(lbl, mask):
    """Remove labels whose centroid is outside the mask."""
    out = np.zeros_like(lbl, dtype=lbl.dtype)
    for rp in measure.regionprops(lbl):
        y, x = rp.centroid
        if 0 <= int(y) < mask.shape[0] and 0 <= int(x) < mask.shape[1] and mask[int(y), int(x)]:
            out[lbl == rp.label] = rp.label
    return measure.label(out, connectivity=2)

# =========================
# Scale bar calibration
# =========================
from skimage import color as _color, morphology as _morph, filters as _filters, measure as _measure

def detect_scale_bar_pixels(overlap_img, search_height_pct: float=0.18):
    H, W = overlap_img.shape[:2]
    h0 = int(H * (1.0 - search_height_pct))
    strip = overlap_img[h0:, :]
    g = _color.rgb2gray(strip) if strip.ndim == 3 else util.img_as_float(strip)
    thr = _filters.threshold_otsu(g)
    bw1 = _morph.binary_opening(g > thr, _morph.rectangle(3,7))
    bw2 = _morph.binary_opening(g < thr, _morph.rectangle(3,7))
    best = 0
    for bw in (bw1, bw2):
        lab = _measure.label(bw, connectivity=2)
        for rp in _measure.regionprops(lab):
            minr, minc, maxr, maxc = rp.bbox
            height = maxr - minr
            width = maxc - minc
            if height < 0.15*strip.shape[0] and width > best:
                best = width
    return best if best>0 else None

def _find_overlap_hint_for_fovs(parent: Path, fovs: dict):
    first_key = next(iter(fovs.keys()))
    patt = re.compile(rf"^{first_key.root}_OVERLAP_.*\.(?:tif|tiff|png)$")
    for p in parent.iterdir():
        if p.is_file() and patt.match(p.name):
            return str(p)
    return None

def calibrate_microns_per_pixel(cfg: PipelineConfig, fovs: dict) -> float:
    if cfg.microns_per_pixel is not None:
        return float(cfg.microns_per_pixel)
    overlap_path = cfg.overlap_image_hint or _find_overlap_hint_for_fovs(Path(cfg.parent_dir), fovs)
    if overlap_path is None:
        raise ValueError("Provide microns_per_pixel or overlap_image_hint (path to a Designer OVERLAP image).")
    img = load_image(overlap_path)
    bar_px = detect_scale_bar_pixels(img, cfg.scale_bar_search_height_pct)
    if not bar_px:
        raise RuntimeError("Scale bar not detected. Provide microns_per_pixel manually or adjust search params.")
    mpp = cfg.scale_bar_microns / float(bar_px)
    if cfg.verbose:
        print(f"Scale bar: {bar_px:.1f} px for {cfg.scale_bar_microns} µm => {mpp:.4f} µm/px")
    return mpp


# =========================
# Segmentation
# =========================
def segment_cells_phase(phase_img, cfg: PipelineConfig):
    if cfg.use_cellpose and _HAS_CELLPOSE:
        # CP4: use CellposeModel; let it pick CUDA 0 if available
        model = models.CellposeModel(
            gpu=True,
            device=None,
            model_type=cfg.cellpose_model
        )
        masks, flows, styles = model.eval(
            phase_img,
            channels=[0, 0],
            diameter=cfg.cellpose_diameter,
            flow_threshold=cfg.cellpose_flow_threshold,
            cellprob_threshold=cfg.cellpose_cellprob_threshold,
            batch_size=cfg.cellpose_batch_size,
            normalize=cfg.normalize_cellpose,
            invert=cfg.invert_cellpose
        )
        return measure.label(masks > 0, connectivity=2).astype(np.int32)
    # Classic fallback
    img = phase_img.astype(np.float32)
    if cfg.phase_gaussian_sigma > 0:
        img = filters.gaussian(img, cfg.phase_gaussian_sigma)
    def _bin(im):
        method = cfg.phase_thresh_method
        if method == "yen": thr = filters.threshold_yen(im)
        elif method == "tri": thr = filters.threshold_triangle(im)
        else: thr = filters.threshold_otsu(im)
        bw = im > thr
        bw = remove_small_objects(bw, cfg.phase_min_size_px)
        bw = morphology.binary_closing(bw, morphology.disk(2))
        return bw
    bw1 = _bin(img)
    bw2 = _bin(1.0 - rescale_for_display(img))
    bw = bw1 if bw1.sum() >= bw2.sum() else bw2
    dist = ndi.distance_transform_edt(bw)
    peaks = feature.peak_local_max(dist, labels=bw, footprint=np.ones((3,3)), exclude_border=False)
    markers = np.zeros_like(dist, dtype=np.int32)
    for i,(r,c) in enumerate(peaks, start=1):
        markers[r,c] = i
    lbl = segmentation.watershed(-dist, markers, mask=bw)
    return lbl.astype(np.int32)

def _robust_01(x, lo=0.5, hi=99.7):
    x = x.astype(np.float32)
    p0, p1 = np.percentile(x, [lo, hi])
    if p1 <= p0:
        return np.zeros_like(x, dtype=np.float32)
    y = (x - p0) / (p1 - p0)
    return np.clip(y, 0, 1).astype(np.float32)

def _adaptive_nuc_threshold(f, cfg):
    """
    Adaptive nuclei threshold.
    - Primary: Sauvola (good for uneven background)
    - Fallback: global Li (then Yen) if Sauvola yields too few pixels
    Applies cfg.nuc_threshold_offset to tighten/loosen masks.
    """
    # Primary: Sauvola
    try:
        from skimage.filters import threshold_sauvola
        win = getattr(cfg, "nuc_sauvola_window", 25)
        k   = getattr(cfg, "nuc_sauvola_k", 0.2)
        thr_map = threshold_sauvola(f, window_size=win, k=k)
        thr_map = thr_map + getattr(cfg, "nuc_threshold_offset", 0.0)   # << tighten
        bw = f > thr_map
    except Exception:
        bw = None

    # Fallback if Sauvola missing/too small
    min_px = max(50, int(getattr(cfg, "nuc_min_pixels_after_thresh", 50)))
    if bw is None or bw.sum() < min_px:
        try:
            t = filters.threshold_li(f)
        except Exception:
            t = filters.threshold_yen(f)
        t = t + getattr(cfg, "nuc_threshold_offset", 0.0)                # << tighten
        bw = f > t

    return bw


def segment_nuclei_single_channel(img, cfg: PipelineConfig):
    # 0) Robust [0,1] normalize (handles 32-bit floats)
    f = _robust_01(img)

    # 1) CLAHE boosts dim nuclei without blowing background
    f = exposure.equalize_adapthist(f, clip_limit=cfg.nuc_clahe_clip)

    # 2) Rolling-ball background subtraction (fallback: large Gaussian)
    try:
        bg = restoration.rolling_ball(f, radius=max(20, cfg.nuc_top_hat_radius*3))
        f = np.clip(f - bg, 0, 1)
    except Exception:
        f = np.clip(f - filters.gaussian(f, sigma=cfg.nuc_top_hat_radius), 0, 1)

    # 3) White tophat to pop small bright blobs
    try:
        f = white_tophat(f, footprint=disk(cfg.nuc_top_hat_radius))
    except TypeError:
        f = white_tophat(f, selem=disk(cfg.nuc_top_hat_radius))

    # 4) Adaptive threshold (with smart fallback)  -------------------------
    #    NOTE: we pass cfg so nuc_threshold_offset / window / k are applied.
    bw = _adaptive_nuc_threshold(f, cfg)

    #    >>> ADD: fix donut holes + gentle size trim BEFORE seeding
    bw = remove_small_holes(bw, getattr(cfg, "nuc_fill_holes_area", 128))   # e.g., 64–256
    if getattr(cfg, "nuc_open_radius", 0) > 0:
        bw = morphology.opening(bw, morphology.disk(cfg.nuc_open_radius))   # try 1
    # ----------------------------------------------------------------------

    # 5) Clean + multiscale LoG seeds → watershed split
    bw = remove_small_objects(bw, max(1, cfg.nuclei_min_size_px))
    bw = remove_small_holes(bw, 16)

    sigmas = np.linspace(cfg.nuc_log_sigma_min, cfg.nuc_log_sigma_max, 4)
    log_resp = np.stack([filters.laplace(filters.gaussian(f, s)) * (s**2) for s in sigmas], axis=0)
    seeds = (log_resp.min(axis=0) < -cfg.nuc_log_threshold)
    seeds = morphology.label(morphology.binary_opening(seeds, morphology.disk(1)))

    # Fallback seeding if LoG yields nothing
    if seeds.max() == 0:
        dist = ndi.distance_transform_edt(bw)
        peaks = feature.peak_local_max(
            dist, labels=bw,
            min_distance=getattr(cfg, "nuc_min_distance_px", 6),  # <<< NEW: prevents merges
            exclude_border=False
        )
        seeds = np.zeros_like(bw, dtype=np.int32)
        for i, (r, c) in enumerate(peaks, start=1):
            seeds[r, c] = i

    # 6) REPLACE distance-based watershed with GRADIENT-based watershed -----
    grad = filters.sobel(f)  # edges around nuclei
    lbl = segmentation.watershed(
        grad,                     # <<< was -distance
        seeds,
        mask=bw,
        compactness=cfg.nuclei_watershed_compactness,
        watershed_line=False
    )
    # ----------------------------------------------------------------------

    return lbl.astype(np.int32)


def merge_nuclei_classes(lbl_g, lbl_r):
    # Merge by union, then assign class by overlap area fraction
    union = (lbl_g > 0) | (lbl_r > 0)
    lbl = measure.label(union, connectivity=2)
    classes = {}
    for rp in measure.regionprops(lbl):
        mask = (lbl == rp.label)
        g_area = (lbl_g > 0)[mask].sum()
        r_area = (lbl_r > 0)[mask].sum()
        if g_area > 0 and r_area > 0:
            cls = "DOUBLE"
        elif g_area > 0:
            cls = "GFP"
        elif r_area > 0:
            cls = "MC"
        else:
            continue
        classes[rp.label] = cls
    return lbl.astype(np.int32), classes


# =========================
# Candidate selection
# =========================
def label_overlaps(lbl_a, lbl_b):
    overlaps = {}
    for a in np.unique(lbl_a)[1:]:
        hits = np.unique(lbl_b[lbl_a==a])
        hits = [int(h) for h in hits if h!=0]
        overlaps[int(a)] = hits
    return overlaps

def find_single_cells(lbl_cells, lbl_nuclei, cfg: PipelineConfig):
    ov = label_overlaps(lbl_cells, lbl_nuclei)
    singles = [cid for cid,hits in ov.items() if len(hits)==1]
    if cfg.exclude_phase_edge_touching:
        H,W = lbl_cells.shape
        _out = []
        for cid in singles:
            rr,cc = np.where(lbl_cells==cid)
            if rr.min()==0 or cc.min()==0 or rr.max()==H-1 or cc.max()==W-1:
                continue
            _out.append(cid)
        singles = _out
    return singles

def grid_coords(x, y, H, W, rows, cols):
    r_sz = H/rows; c_sz = W/cols
    r = min(int(y//r_sz), rows-1); c = min(int(x//c_sz), cols-1)
    return (r,c)

def sample_spatially_diverse_cells(lbl_cells, cell_ids, target_n, rows, cols, seed=1337):
    rng = np.random.RandomState(seed)
    props = {r.label:r for r in measure.regionprops(lbl_cells) if r.label in cell_ids}
    H,W = lbl_cells.shape
    bins = {(r,c):[] for r in range(rows) for c in range(cols)}
    for cid, rp in props.items():
        y,x = rp.centroid
        bins[grid_coords(x,y,H,W,rows,cols)].append(cid)
    for k in bins: rng.shuffle(bins[k])
    selection = []
    backups = {k: bins[k][:] for k in bins}
    while len(selection)<target_n and any(len(v)>0 for v in bins.values()):
        for k in sorted(bins.keys()):
            if len(selection)>=target_n: break
            if bins[k]: selection.append(bins[k].pop())
    return selection, backups


# =========================
# QC & Measurements
# =========================
def crop_around_label(img, lbl, label_id, pad=20):
    rr,cc = np.where(lbl==label_id)
    if rr.size==0:
        return img, (0,img.shape[0],0,img.shape[1])
    r0,r1 = rr.min(), rr.max()
    c0,c1 = cc.min(), cc.max()
    r0 = max(0, r0-pad); c0 = max(0, c0-pad)
    r1 = min(img.shape[0]-1, r1+pad); c1 = min(img.shape[1]-1, c1+pad)
    return img[r0:r1+1, c0:c1+1], (r0, r1+1, c0, c1+1)

def make_qc_panel(cell_id, phase_img, lbl_cells, lbl_nuclei, nuc_classes, green_img, red_img, cfg: PipelineConfig, overlap_img=None):
    # crop around the cell on PHASE
    crop_phase, bbox = crop_around_label(phase_img, lbl_cells, cell_id, cfg.qc_panel_crop_pad_px)
    r0,r1,c0,c1 = bbox

    # cropped masks
    crop_cellmask = (lbl_cells[r0:r1, c0:c1]==cell_id).astype(np.uint8)
    crop_nuclbl   = lbl_nuclei[r0:r1, c0:c1]
    crop_nucmask  = (crop_nuclbl > 0).astype(np.uint8)

    # raw fluorescence composite (max of R/G)
    crop_green = green_img[r0:r1, c0:c1]
    crop_red   = red_img[r0:r1, c0:c1]
    g_disp = rescale_for_display(crop_green)
    r_disp = rescale_for_display(crop_red)
    nuc_raw_rgb = np.dstack([r_disp, g_disp, np.zeros_like(g_disp)])

    # real OVERLAP crop if provided (else synthesize from r/g + phase)
    if overlap_img is not None:
        crop_overlap = overlap_img[r0:r1, c0:c1]
        if crop_overlap.ndim == 2:
            crop_overlap = np.stack([crop_overlap]*3, axis=-1)
        ov_disp = rescale_for_display(crop_overlap)
    else:
        ph = rescale_for_display(crop_phase)
        ov_disp = np.dstack([np.maximum(r_disp, ph*0.2), np.maximum(g_disp, ph*0.2), ph*0.2])

    outline = segmentation.find_boundaries(crop_cellmask, mode="inner")
    ov2 = ov_disp.copy()
    ov2[outline] = [1,1,1]

    fig, axes = plt.subplots(1,5, figsize=(12,3))
    axes[0].imshow(rescale_for_display(crop_phase), cmap="gray"); axes[0].set_title(f"Phase (raw) #{cell_id}")
    axes[1].imshow(crop_cellmask, cmap="viridis"); axes[1].set_title("Phase Seg")
    axes[2].imshow(nuc_raw_rgb); axes[2].set_title("Raw Nuclei (R/G)")
    axes[3].imshow(crop_nucmask, cmap="magma"); axes[3].set_title("Nuc Seg (merged)")
    axes[4].imshow(ov2); axes[4].set_title("OVERLAP + Cell outline")
    for ax in axes: ax.axis("off")
    fig.tight_layout()
    return fig

def _write_selection_manifest(selected_items, out_path: Path):
    """Save reviewed selection with centroid so we can robustly reload."""
    out = []
    for idx, it in enumerate(selected_items, start=1):
        k   = it["key"]
        cid = int(it["cid"])
        lbl = it["lbl_cells"]
        rp  = next((r for r in measure.regionprops(lbl) if r.label == cid), None)
        cy, cx = (float(rp.centroid[0]), float(rp.centroid[1])) if rp is not None else (None, None)
        area = int(rp.area) if rp is not None else None
        out.append({
            "selected_rank": idx,
            "root": k.root,
            "sample": k.sample,
            "well": k.well,
            "site": k.site,
            "timecode": k.timecode,
            "cell_id": cid,
            "centroid_rc": [cy, cx],     # << NEW (row, col) in full-image coords
            "area_px": area              # optional, helpful for debugging
        })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

def _indices_from_manifest(candidates, manifest_path: str, contain_check=True, nearest_px=10):
    """
    Return candidate indices in the manifest's order.
    Strategy:
      1) Exact key + cid match (fast path)
      2) If missing: try "label that contains the saved centroid"
      3) If still missing: nearest centroid within `nearest_px`
    """
    with open(manifest_path, "r") as f:
        items = json.load(f)

    # Build exact-match index: (root, sample, well, site, timecode, cid) -> cand_idx
    exact = {}
    # Also precompute centroids for fallback
    cand_key = []
    cand_centroids = []
    for i, c in enumerate(candidates):
        k = c["key"]; cid = int(c["cid"])
        exact[(k.root, k.sample, k.well, k.site, k.timecode, cid)] = i
        # centroid for fallback
        rp = next((r for r in measure.regionprops(c["lbl_cells"]) if r.label == cid), None)
        if rp is not None:
            cand_centroids.append((float(rp.centroid[0]), float(rp.centroid[1])))
        else:
            cand_centroids.append((np.nan, np.nan))
        cand_key.append((k.root, k.sample, k.well, k.site, k.timecode))

    out_indices = []
    misses = 0

    for it in items:
        root, sample, well, site, tcode = it["root"], it["sample"], it["well"], it["site"], it["timecode"]
        cid = int(it.get("cell_id", -1))
        cy, cx = (it.get("centroid_rc") or [None, None])

        # (1) exact lookup first
        tup = (root, sample, well, site, tcode, cid)
        if tup in exact:
            out_indices.append(exact[tup])
            continue

        # Restrict to same FOV
        fov_idxs = [i for i, key in enumerate(cand_key) if key == (root, sample, well, site, tcode)]
        if not fov_idxs:
            misses += 1
            continue

        # (2) label containing saved centroid
        if contain_check and cy is not None and cx is not None:
            found = None
            for i in fov_idxs:
                lbl = candidates[i]["lbl_cells"]
                rr, cc = int(round(cy)), int(round(cx))
                if 0 <= rr < lbl.shape[0] and 0 <= cc < lbl.shape[1]:
                    # if the same object was re-labeled, the pixel at centroid should still be that object
                    if lbl[rr, cc] == candidates[i]["cid"]:
                        found = i
                        break
                    # otherwise accept any non-zero label containing that point
                    if lbl[rr, cc] > 0:
                        found = i
                        break
            if found is not None:
                out_indices.append(found)
                continue

        # (3) nearest-centroid within threshold
        if cy is not None and cx is not None:
            best_i, best_d2 = None, None
            for i in fov_idxs:
                yx = cand_centroids[i]
                if np.isnan(yx[0]): 
                    continue
                d2 = (yx[0]-cy)**2 + (yx[1]-cx)**2
                if best_d2 is None or d2 < best_d2:
                    best_i, best_d2 = i, d2
            if best_i is not None and best_d2 <= (nearest_px**2):
                out_indices.append(best_i)
                continue

        # give up for this item
        misses += 1

    return out_indices


def _cell_bbox(lbl, cid, pad=0):
    rr, cc = np.where(lbl == cid)
    if rr.size == 0:
        return None
    r0, r1 = rr.min(), rr.max() + 1
    c0, c1 = cc.min(), cc.max() + 1
    H, W = lbl.shape
    r0 = max(0, r0 - pad); c0 = max(0, c0 - pad)
    r1 = min(H, r1 + pad); c1 = min(W, c1 + pad)
    return (r0, r1, c0, c1)

def _compose_overlap_like(phase_img, green_img, red_img, overlap_img=None):
    """Return an RGB image similar to the Designer OVERLAP."""
    if overlap_img is not None:
        base = overlap_img
        if base.ndim == 2:  # grayscale → RGB
            base = np.stack([base]*3, axis=-1)
        return rescale_for_display(base)
    # synthesize from phase + channels
    ph = rescale_for_display(phase_img)
    g  = rescale_for_display(green_img)
    r  = rescale_for_display(red_img)
    return np.dstack([np.maximum(r, ph*0.25), np.maximum(g, ph*0.25), ph*0.25])

def _save_fov_selection_map(sample_id, fov_key, base_rgb, selections_for_fov, out_dir, cfg):
    """
    Draw red rectangles (and optional index labels) on the OVERLAP-like image for one FOV.
    selections_for_fov: list of dicts with keys:
        - 'rank' (1-based global index)
        - 'cid', 'lbl_cells'
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.imshow(base_rgb)
    ax.axis("off")

    color = cfg.selection_map_box_color
    for sel in selections_for_fov:
        cid = sel["cid"]
        lbl = sel["lbl_cells"]
        bbox = _cell_bbox(lbl, cid, pad=cfg.selection_map_pad_px)
        if bbox is None:
            continue
        r0, r1, c0, c1 = bbox
        h = r1 - r0; w = c1 - c0
        rect = Rectangle((c0, r0), w, h, fill=False,
                         edgecolor=color, linewidth=cfg.selection_map_box_thickness)
        ax.add_patch(rect)
        if cfg.selection_map_label_boxes:
            ax.text(c0+2, r0+12, str(sel["rank"]),
                    color="white", fontsize=9,
                    bbox=dict(facecolor="black", alpha=0.6, edgecolor="none", boxstyle="round,pad=0.15"))

    fname = f"{sample_id}_FOV_{fov_key.well}_site{fov_key.site}_{fov_key.timecode}_selection_map.png"
    out_path = Path(out_dir) / fname
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

def imagej_metrics_from_region(region, mpp):
    area_px = float(region.area)
    perimeter_px = float(region.perimeter)
    major_px = float(region.major_axis_length) if region.major_axis_length else 0.0
    minor_px = float(region.minor_axis_length) if region.minor_axis_length else 0.0
    orient = float(region.orientation) if region.orientation is not None else 0.0
    angle_deg = -orient * 180.0 / math.pi
    area_um2 = area_px * (mpp**2)
    perimeter_um = perimeter_px * mpp
    major_um = major_px * mpp
    minor_um = minor_px * mpp
    circ = (4*math.pi*area_px)/(perimeter_px**2) if perimeter_px>0 else 0.0
    ar = (major_px/minor_px) if minor_px>0 else 0.0
    rnd = (4*area_px)/(math.pi*(major_px**2)) if major_px>0 else 0.0
    convex_area_px = float(region.convex_area) if hasattr(region, "convex_area") else float(morphology.convex_hull_image(region.image).sum())
    sol = (area_px/convex_area_px) if convex_area_px>0 else 0.0
    return {
        "area_px": area_px, "perimeter_px": perimeter_px, "major_px": major_px, "minor_px": minor_px, "angle_deg": angle_deg,
        "area_um2": area_um2, "perimeter_um": perimeter_um, "major_um": major_um, "minor_um": minor_um,
        "circularity": circ, "aspect_ratio": ar, "roundness": rnd, "solidity": sol,
        "centroid_row": float(region.centroid[0]), "centroid_col": float(region.centroid[1]),
    }

def refine_nucleus_inside_cell(green, red, cell_mask, cfg: PipelineConfig):
    """Re-segment nuclei inside one cell ROI to enforce single-nucleus rule."""
    roi = cell_mask.astype(bool)
    combo = np.maximum(_robust_01(green), _robust_01(red))
    combo = exposure.equalize_adapthist(combo, clip_limit=max(0.008, cfg.nuc_clahe_clip))
    try:
        bg = restoration.rolling_ball(combo, radius=max(18, cfg.nuc_top_hat_radius*2))
        combo = np.clip(combo - bg, 0, 1)
    except Exception:
        combo = np.clip(combo - filters.gaussian(combo, sigma=cfg.nuc_top_hat_radius), 0, 1)
    try:
        combo = white_tophat(combo, footprint=disk(max(3, cfg.nuc_top_hat_radius//2)))
    except TypeError:
        combo = white_tophat(combo, selem=disk(max(3, cfg.nuc_top_hat_radius//2)))

    combo *= roi
    bw = _adaptive_nuc_threshold(combo)
    bw &= roi
    bw = remove_small_objects(bw, 20)
    bw = remove_small_holes(bw, 8)

    dist = ndi.distance_transform_edt(bw)
    peaks = feature.peak_local_max(dist, labels=bw, footprint=np.ones((3,3)), exclude_border=False)
    seeds = np.zeros_like(dist, dtype=np.int32)
    for i,(r,c) in enumerate(peaks, start=1):
        seeds[r,c] = i
    lbl = segmentation.watershed(-dist, seeds, mask=bw)
    return lbl.astype(np.int32)

def measure_cell_and_nucleus(lbl_cells, lbl_nuclei, nuc_classes, selected_ids, mpp, cfg: PipelineConfig,
                             green_img=None, red_img=None):
    """Measure cell + nucleus; optionally refine nucleus inside each selected cell using raw GFP/MC."""
    cell_props = {r.label:r for r in measure.regionprops(lbl_cells) if r.label in selected_ids}
    out_cell, out_nuc = [], []
    H,W = lbl_cells.shape
    for cid in selected_ids:
        cprop = cell_props[cid]
        if cfg.exclude_phase_edge_touching:
            minr,minc,maxr,maxc = cprop.bbox
            if minr==0 or minc==0 or maxr==H or maxc==W: 
                continue
        mask = (lbl_cells==cid)

        hits = np.unique(lbl_nuclei[mask]); hits = [int(h) for h in hits if h!=0]

        # optional refinement if not exactly one nucleus overlapped
        if (len(hits) != 1) and cfg.refine_nuclei_within_cell and (green_img is not None) and (red_img is not None):
            rr, cc = np.where(mask)
            r0, r1 = rr.min(), rr.max()+1
            c0, c1 = cc.min(), cc.max()+1
            roi_cell = mask[r0:r1, c0:c1]
            g_roi = green_img[r0:r1, c0:c1]
            r_roi = red_img[r0:r1, c0:c1]
            lbl_ref = refine_nucleus_inside_cell(g_roi, r_roi, roi_cell, cfg)
            tmp = np.zeros_like(lbl_nuclei)
            tmp[r0:r1, c0:c1] = lbl_ref
            hits = np.unique(tmp[mask]); hits = [int(h) for h in hits if h!=0]

        if len(hits)!=1:
            continue
        nid = hits[0]

        if cfg.exclude_nucleus_edge_touching:
            nprop = [r for r in measure.regionprops(lbl_nuclei) if r.label==nid][0]
            minr,minc,maxr,maxc = nprop.bbox
            if minr==0 or minc==0 or maxr==H or maxc==W:
                continue

        cm = imagej_metrics_from_region(cprop, mpp); cm.update({"cell_id": int(cid)})
        nm = imagej_metrics_from_region([r for r in measure.regionprops(lbl_nuclei) if r.label==nid][0], mpp)
        nm.update({"cell_id": int(cid), "nucleus_id": int(nid), "nucleus_class": nuc_classes.get(nid,"NA")})
        out_cell.append(cm); out_nuc.append(nm)
    return pd.DataFrame(out_cell), pd.DataFrame(out_nuc)


# =========================
# Orchestration (auto)
# =========================
# def smooth_label_mask(lbl, radius=1):
#     """Light morphological smoothing per label (close→open)."""
#     if radius <= 0:
#         return lbl
#     out = np.zeros_like(lbl, dtype=lbl.dtype)
#     se1 = disk(radius)
#     se2 = disk(max(1, radius-1))
#     for lab in np.unique(lbl)[1:]:
#         m = (lbl == lab)
#         m = morphology.binary_closing(m, se1)
#         m = morphology.binary_opening(m, se2)
#         out[m] = lab
#     return out

def smooth_label_mask_subpixel(lbl, sigma=0.8, pad=3):
    """Round label edges using signed distance + Gaussian blur (ROI-based)."""
    if sigma <= 0:
        return lbl
    out = np.zeros_like(lbl, dtype=lbl.dtype)
    labs = np.unique(lbl); labs = labs[labs != 0]
    for lab in labs:
        m = (lbl == lab)
        rr, cc = np.where(m)
        if rr.size == 0:
            continue
        r0, r1 = max(0, rr.min()-pad), min(lbl.shape[0], rr.max()+1+pad)
        c0, c1 = max(0, cc.min()-pad), min(lbl.shape[1], cc.max()+1+pad)
        mr = m[r0:r1, c0:c1]
        # signed distance inside/outside on ROI only
        d_in  = ndi.distance_transform_edt(mr)
        d_out = ndi.distance_transform_edt(~mr)
        sdf = d_in - d_out
        sdf = ndi.gaussian_filter(sdf.astype(np.float32), sigma=sigma)
        m2 = sdf > 0
        out[r0:r1, c0:c1][m2] = lab
    return out


def process_sample(sample_id: str, cfg: PipelineConfig):
    np.random.seed(cfg.random_seed)

    parent = Path(cfg.parent_dir)
    out_root = Path(cfg.outputs_dir) / sample_id
    out_qc = out_root / "qc"
    out_masks = out_root / "masks"
    ensure_dir(out_root); ensure_dir(out_qc); ensure_dir(out_masks)

    fovs = group_fovs_by_key(parent, cfg, sample_filter=sample_id, max_fovs=cfg.max_fovs_per_sample)
    if not fovs: 
        raise RuntimeError(f"No FOVs found for sample {sample_id}. Check names and regex.")
    if cfg.verbose: print(f"FOVs for {sample_id}: {len(fovs)}")

    mpp = calibrate_microns_per_pixel(cfg, fovs)

    # Gather candidates
    candidates = []
    for key, chans in fovs.items():
        phase = load_image(chans["PHASE"])
        g = load_image(chans["GFP"])
        r = load_image(chans["MC"])
        overlap_img = load_image(chans["OVERLAP"]) if "OVERLAP" in chans else None

        # --- Inner circle crop (well edge exclusion) ---
        roi_mask = None
        if cfg.use_inner_circle_crop:
            roi_mask = _circle_mask(phase.shape, radius_frac=cfg.crop_radius_frac, center_xy=cfg.crop_center_xy)
            if cfg.hard_mask_images:
                phase, g, r = _apply_roi_mask_to_images(roi_mask, phase, g, r)
        # -----------------------------------------------

        lbl_cells = segment_cells_phase(phase, cfg)
        lbl_g = segment_nuclei_single_channel(g, cfg)
        lbl_r = segment_nuclei_single_channel(r, cfg)
        lbl_nuclei, nuc_classes = merge_nuclei_classes(lbl_g, lbl_r)

        # size filter in µm²
        keep_cells = np.zeros_like(lbl_cells, dtype=bool)
        for rp in measure.regionprops(lbl_cells):
            a = rp.area * (mpp**2)
            if cfg.min_cell_area_um2 <= a <= cfg.max_cell_area_um2:
                keep_cells[lbl_cells==rp.label] = True
        lbl_cells = measure.label(keep_cells, connectivity=2)

        keep_n = np.zeros_like(lbl_nuclei, dtype=bool)
        for rp in measure.regionprops(lbl_nuclei):
            a = rp.area * (mpp**2)
            if cfg.min_nuc_area_um2 <= a <= cfg.max_nuc_area_um2:
                keep_n[lbl_nuclei==rp.label] = True
        lbl_nuclei = measure.label(keep_n, connectivity=2)

        # --- Drop labels outside ROI (centroid rule) ---
        if roi_mask is not None:
            lbl_cells  = _keep_labels_inside_mask(lbl_cells,  roi_mask)
            lbl_nuclei = _keep_labels_inside_mask(lbl_nuclei, roi_mask)
        # -----------------------------------------------

        # --- apply subpixel smoothing to working masks ---
        if not cfg.smooth_for_display_only:
            lbl_cells  = smooth_label_mask_subpixel(lbl_cells,  sigma=0.8)   # try 0.6–1.2
            lbl_nuclei = smooth_label_mask_subpixel(lbl_nuclei, sigma=0.8)
        # --------------------------------------------------

        singles = find_single_cells(lbl_cells, lbl_nuclei, cfg)
        for cid in singles:
            candidates.append({
                "key": key, "cid": cid,
                "lbl_cells": lbl_cells, "lbl_nuclei": lbl_nuclei,
                "nuc_classes": nuc_classes,
                "phase": phase, "green": g, "red": r,
                "overlap": overlap_img,
            })
        if cfg.save_intermediate_masks:
            tiff.imwrite(str(out_masks / f"{key.root}_{sample_id}_{key.well}_site{key.site}_{key.timecode}_cells.tif"), lbl_cells.astype(np.uint16))
            tiff.imwrite(str(out_masks / f"{key.root}_{sample_id}_{key.well}_site{key.site}_{key.timecode}_nuclei.tif"), lbl_nuclei.astype(np.uint16))

    if cfg.verbose: print(f"Single-cell candidates: {len(candidates)}")

    # Spatial selection per FOV, then concat until target
    needed = cfg.target_cells
    selected = []
    for key, _ch in fovs.items():
        subs = [c for c in candidates if c["key"]==key]
        if not subs: continue
        lbl_cells = subs[0]["lbl_cells"]
        ids = [c["cid"] for c in subs]
        per_target = max(1, math.ceil(cfg.target_cells/max(1,len(fovs))))
        sel, _ = sample_spatially_diverse_cells(lbl_cells, ids, min(needed, per_target), cfg.grid_rows, cfg.grid_cols, cfg.random_seed)
        selected.extend([c for c in subs if c["cid"] in sel])
        needed = cfg.target_cells - len(selected)
        if needed <= 0: break

    if len(selected)==0:
        raise RuntimeError("No single-cell selections could be made. Loosen filters or increase max_fovs_per_sample.")
    selected = selected[:cfg.target_cells]
    if cfg.verbose: print(f"Selected {len(selected)} cells for QC/metrics.")

    # Save QC panels + gallery
    if cfg.save_individual_qc_panels:
        for idx, item in enumerate(selected, start=1):
            fig = make_qc_panel(
                item["cid"], item["phase"], item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
                item["green"], item["red"], cfg, overlap_img=item.get("overlap")
            )
            fig.savefig(out_qc / f"{sample_id}_cell_{idx:03d}_panel.png", dpi=150)
            plt.close(fig)
    cols = cfg.qc_panels_per_row
    rows = math.ceil(len(selected) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols*3, rows*3))
    axes = np.atleast_1d(axes).ravel()
    for i, item in enumerate(selected):
        fig_i = make_qc_panel(
            item["cid"], item["phase"], item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
            item["green"], item["red"], cfg, overlap_img=item.get("overlap")
        )
        fig_i.canvas.draw()
        w, h = fig_i.canvas.get_width_height()
        buf = np.frombuffer(fig_i.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
        axes[i].imshow(buf[..., :3]); axes[i].axis("off"); axes[i].set_title(f"{i+1}")
        plt.close(fig_i)
    # turn off unused axes
    for j in range(i+1, len(axes)): axes[j].axis("off")
    fig.tight_layout()
    if cfg.save_combined_qc_panel:
        fig.savefig(out_qc / f"{sample_id}_selected_cells_QC.png", dpi=150)
    plt.close(fig)

    # --- NEW: per-FOV selection maps (red boxes around chosen cells) ---
    if cfg.save_selection_maps:
        # group selected items by FOV (ImageKey)
        by_fov = {}
        for rank, item in enumerate(selected, start=1):
            k = item["key"]
            by_fov.setdefault(k, []).append({
                "rank": rank,
                "cid": item["cid"],
                "lbl_cells": item["lbl_cells"],
                "phase": item["phase"],
                "green": item["green"],
                "red": item["red"],
                "overlap": item.get("overlap"),
            })
        # compose base image once per FOV and draw all boxes
        for k, items in by_fov.items():
            base = _compose_overlap_like(items[0]["phase"], items[0]["green"], items[0]["red"], items[0]["overlap"])
            _save_fov_selection_map(sample_id, k, base, items, out_qc, cfg)
    # -------------------------------------------------------------------

    # Measurements & export
    rows_cell, rows_nuc = [], []
    for idx, item in enumerate(selected, start=1):
        key = item["key"]; cid = item["cid"]
        dfc, dfn = measure_cell_and_nucleus(
            item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
            [cid], mpp, cfg, green_img=item["green"], red_img=item["red"]
        )
        for d in dfc.to_dict("records"):
            d.update({"root": key.root, "sample": key.sample, "well": key.well, "site": key.site, "timecode": key.timecode, "selected_rank": idx})
            rows_cell.append(d)
        for d in dfn.to_dict("records"):
            d.update({"root": key.root, "sample": key.sample, "well": key.well, "site": key.site, "timecode": key.timecode, "selected_rank": idx})
            rows_nuc.append(d)
    df_cell = pd.DataFrame(rows_cell); df_nuc = pd.DataFrame(rows_nuc)
    ts = __import__("datetime").datetime.now().strftime("%Y%m%d_%H%M%S")
    cell_csv = out_root / f"cell_morphology_{sample_id}_{ts}.csv"
    nuc_csv = out_root / f"nuclear_morphology_{sample_id}_{ts}.csv"
    df_cell.to_csv(cell_csv, index=False); df_nuc.to_csv(nuc_csv, index=False)
    if cfg.verbose:
        print("Wrote:", cell_csv); print("Wrote:", nuc_csv); print("QC dir:", out_qc)
    return {"mpp": mpp, "cell_csv": str(cell_csv), "nuc_csv": str(nuc_csv), "qc_dir": str(out_qc)}

# =========================
# Interactive review helpers
# =========================

def _gather_candidates(sample_id: str, cfg: PipelineConfig):
    """Internal: gather FOVs, compute mpp, segment all channels, and build candidate list."""
    parent = Path(cfg.parent_dir)
    out_root = Path(cfg.outputs_dir) / sample_id
    out_masks = out_root / "masks"
    ensure_dir(out_root); ensure_dir(out_masks)

    fovs = group_fovs_by_key(parent, cfg, sample_filter=sample_id, max_fovs=cfg.max_fovs_per_sample)
    if not fovs:
        raise RuntimeError(f"No FOVs found for sample {sample_id}. Check names and regex.")

    mpp = calibrate_microns_per_pixel(cfg, fovs)

    candidates = []
    for key, chans in fovs.items():
        phase = load_image(chans["PHASE"])
        g = load_image(chans["GFP"])
        r = load_image(chans["MC"])
        overlap_img = load_image(chans["OVERLAP"]) if "OVERLAP" in chans else None

        # --- Inner circle crop (well edge exclusion) ---
        roi_mask = None
        if cfg.use_inner_circle_crop:
            roi_mask = _circle_mask(phase.shape, radius_frac=cfg.crop_radius_frac, center_xy=cfg.crop_center_xy)
            if cfg.hard_mask_images:
                phase, g, r = _apply_roi_mask_to_images(roi_mask, phase, g, r)
        # -----------------------------------------------

        lbl_cells = segment_cells_phase(phase, cfg)
        lbl_g = segment_nuclei_single_channel(g, cfg)
        lbl_r = segment_nuclei_single_channel(r, cfg)
        lbl_nuclei, nuc_classes = merge_nuclei_classes(lbl_g, lbl_r)

        # size filters in µm²
        keep_cells = np.zeros_like(lbl_cells, dtype=bool)
        for rp in measure.regionprops(lbl_cells):
            a = rp.area * (mpp**2)
            if cfg.min_cell_area_um2 <= a <= cfg.max_cell_area_um2:
                keep_cells[lbl_cells==rp.label] = True
        lbl_cells = measure.label(keep_cells, connectivity=2)

        keep_n = np.zeros_like(lbl_nuclei, dtype=bool)
        for rp in measure.regionprops(lbl_nuclei):
            a = rp.area * (mpp**2)
            if cfg.min_nuc_area_um2 <= a <= cfg.max_nuc_area_um2:
                keep_n[lbl_nuclei==rp.label] = True
        lbl_nuclei = measure.label(keep_n, connectivity=2)

        # --- Drop labels outside ROI (centroid rule) ---
        if roi_mask is not None:
            lbl_cells  = _keep_labels_inside_mask(lbl_cells,  roi_mask)
            lbl_nuclei = _keep_labels_inside_mask(lbl_nuclei, roi_mask)
        # -----------------------------------------------

        # --- apply subpixel smoothing to working masks ---
        if not cfg.smooth_for_display_only:
            lbl_cells  = smooth_label_mask_subpixel(lbl_cells,  sigma=0.8)   # try 0.6–1.2
            lbl_nuclei = smooth_label_mask_subpixel(lbl_nuclei, sigma=0.8)
        # --------------------------------------------------

        singles = find_single_cells(lbl_cells, lbl_nuclei, cfg)
        for cid in singles:
            candidates.append({
                "key": key, "cid": cid,
                "lbl_cells": lbl_cells, "lbl_nuclei": lbl_nuclei,
                "nuc_classes": nuc_classes,
                "phase": phase, "green": g, "red": r,
                "overlap": overlap_img,
            })

        if cfg.save_intermediate_masks:
            tiff.imwrite(str(out_masks / f"{key.root}_{sample_id}_{key.well}_site{key.site}_{key.timecode}_cells.tif"), lbl_cells.astype(np.uint16))
            tiff.imwrite(str(out_masks / f"{key.root}_{sample_id}_{key.well}_site{key.site}_{key.timecode}_nuclei.tif"), lbl_nuclei.astype(np.uint16))

    return fovs, mpp, candidates


def build_review_session(sample_id: str, cfg: PipelineConfig):
    """
    Build an interactive session:
    - precomputes candidates for the sample
    - picks an initial spatially-diverse selection (len = cfg.target_cells or less)
    - optionally REPLAYS a prior selection from a manifest (cfg.use_selection_manifest_path)
    - stores per-FOV per-bin backups for replacements
    Returns a dict 'session' to be used by show_batch/replace_selected/finalize_review.
    """
    np.random.seed(cfg.random_seed)
    fovs, mpp, candidates = _gather_candidates(sample_id, cfg)
    if cfg.verbose:
        print(f"Single-cell candidates found: {len(candidates)}")

    # bucket candidates per FOV and grid-bin
    rows, cols = cfg.grid_rows, cfg.grid_cols
    rng = np.random.RandomState(cfg.random_seed)

    # Map: (fov_key, bin_rc) -> [indices into candidates]
    binmap = {}
    for idx, c in enumerate(candidates):
        k = c["key"]
        lbl_cells = c["lbl_cells"]
        rp = next(r for r in measure.regionprops(lbl_cells) if r.label == c["cid"])
        y, x = rp.centroid
        H, W = lbl_cells.shape
        r_sz = H / rows
        c_sz = W / cols
        rbin = min(int(y // r_sz), rows - 1)
        cbin = min(int(x // c_sz), cols - 1)
        binmap.setdefault((k, (rbin, cbin)), []).append(idx)

    # shuffle bin lists
    for key_bin, idxs in binmap.items():
        rng.shuffle(idxs)

    # Initial selection: cycle bins across FOVs for diversity
    selected_idxs = []
    need = cfg.target_cells
    key_bins = list(binmap.keys())
    key_bins.sort(key=lambda x: (x[0].well, x[0].site, x[1][0], x[1][1]))  # deterministic
    while need > 0 and any(len(binmap[kb]) > 0 for kb in key_bins):
        for kb in key_bins:
            if need <= 0:
                break
            if binmap[kb]:
                cand_idx = binmap[kb].pop(0)
                if cand_idx not in selected_idxs:
                    selected_idxs.append(cand_idx)
                    need -= 1

    # Prepare backup pools (remaining items in each bin)
    backups = {kb: idxs[:] for kb, idxs in binmap.items()}

    # ---------- NEW: manifest replay (overrides initial selection if present) ----------
    if getattr(cfg, "use_selection_manifest_path", None):
        mp = Path(cfg.use_selection_manifest_path)
        if mp.exists():
            try:
                # uses centroid-aware matching with fallback to nearest within 10 px (default)
                sel_from_manifest = _indices_from_manifest(candidates, str(mp))
                if sel_from_manifest:
                    selected_idxs = sel_from_manifest[:cfg.target_cells]
                    if cfg.verbose:
                        print(f"Replayed {len(selected_idxs)} selections from manifest.")
                else:
                    if cfg.verbose:
                        print("Manifest replay found 0 matches; using diversity-based selection.")
            except Exception as e:
                if cfg.verbose:
                    print(f"Warning: failed to load selection manifest: {e}")
        else:
            if cfg.verbose:
                print(f"Manifest path not found: {mp}")
    # ----------------------------------------------------------------------------------

    session = {
        "sample": sample_id,
        "cfg": cfg,
        "mpp": mpp,
        "fovs": fovs,
        "candidates": candidates,
        "selected_idxs": selected_idxs[:cfg.target_cells],
        "backups": backups,           # (fov_key, (rbin,cbin)) -> [cand idxs]
        "rows": rows,
        "cols": cols,
        "rng": rng,
    }
    if cfg.verbose:
        print(f"Initial selection: {len(session['selected_idxs'])} cells.")
    return session


def _candidate_bin(cand, rows, cols):
    lbl_cells = cand["lbl_cells"]
    rp = next(r for r in measure.regionprops(lbl_cells) if r.label == cand["cid"])
    y, x = rp.centroid
    H, W = lbl_cells.shape
    r_sz = H/rows; c_sz = W/cols
    rbin = min(int(y//r_sz), rows-1); cbin = min(int(x//c_sz), cols-1)
    return (rbin, cbin)

def draw_qc_row(axes_row, cell_id, phase_img, lbl_cells, lbl_nuclei, nuc_classes,
                green_img, red_img, cfg: PipelineConfig, overlap_img=None):
    """Draw the 5-tile QC row (Phase, Phase Seg, Raw Nuclei, Nuc Seg, OVERLAP+outline) into axes_row[0..4]."""
    # crop around the cell on PHASE
    crop_phase, bbox = crop_around_label(phase_img, lbl_cells, cell_id, cfg.qc_panel_crop_pad_px)
    r0,r1,c0,c1 = bbox

    # cropped masks
    crop_cellmask = (lbl_cells[r0:r1, c0:c1]==cell_id).astype(np.uint8)
    crop_nuclbl   = lbl_nuclei[r0:r1, c0:c1]
    crop_nucmask  = (crop_nuclbl > 0).astype(np.uint8)

    # raw fluorescence composite (R/G)
    crop_green = green_img[r0:r1, c0:c1]
    crop_red   = red_img[r0:r1, c0:c1]
    g_disp = rescale_for_display(crop_green)
    r_disp = rescale_for_display(crop_red)
    nuc_raw_rgb = np.dstack([r_disp, g_disp, np.zeros_like(g_disp)])

    # real OVERLAP crop if provided (else synthesize from r/g + phase)
    if overlap_img is not None:
        crop_overlap = overlap_img[r0:r1, c0:c1]
        if crop_overlap.ndim == 2:
            crop_overlap = np.stack([crop_overlap]*3, axis=-1)
        ov_disp = rescale_for_display(crop_overlap)
    else:
        ph = rescale_for_display(crop_phase)
        ov_disp = np.dstack([np.maximum(r_disp, ph*0.2), np.maximum(g_disp, ph*0.2), ph*0.2])

    outline = segmentation.find_boundaries(crop_cellmask, mode="inner")
    ov2 = ov_disp.copy()
    ov2[outline] = [1,1,1]

    # draw into provided axes
    ax0, ax1, ax2, ax3, ax4 = axes_row
    ax0.imshow(rescale_for_display(crop_phase), cmap="gray"); ax0.set_title(f"Phase (raw) #{cell_id}"); ax0.axis("off")
    ax1.imshow(crop_cellmask, cmap="viridis");                 ax1.set_title("Phase Seg");            ax1.axis("off")
    ax2.imshow(nuc_raw_rgb);                                   ax2.set_title("Raw Nuclei (R/G)");     ax2.axis("off")
    ax3.imshow(crop_nucmask,  cmap="magma");                   ax3.set_title("Nuc Seg (merged)");     ax3.axis("off")
    ax4.imshow(ov2);                                           ax4.set_title("OVERLAP + Cell outline"); ax4.axis("off")

def show_batch(session: dict, start: int = 0, count: int = 5):
    sel = session["selected_idxs"][start:start+count]
    if not sel:
        print("Nothing to show for this range.")
        return

    n = len(sel)
    fig, axes = plt.subplots(n, 5, figsize=(12, 3*n))
    if n == 1:
        axes = np.expand_dims(axes, 0)  # ensure 2D

    total = len(session["selected_idxs"])
    for row, cand_idx in enumerate(sel):
        item = session["candidates"][cand_idx]
        draw_qc_row(
            axes[row, :],
            item["cid"], item["phase"], item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
            item["green"], item["red"], session["cfg"], overlap_img=item.get("overlap")
        )
        # >>> overlay "k/total" in the top-left of the FIRST tile for this row
        global_idx = start + row + 1
        axes[row, 0].text(
            0.02, 0.06, f"{global_idx}/{total}",
            transform=axes[row, 0].transAxes,
            ha="left", va="bottom",
            color="white",
            fontsize=10,
            bbox=dict(facecolor="black", alpha=0.5, edgecolor="none", boxstyle="round,pad=0.2")
        )

    plt.suptitle(f"Selected cells {start+1}–{start+n} (of {total})")
    plt.tight_layout()
    plt.show()

def replace_selected(session: dict, bad_indices: List[int]):
    """
    Replace 1-based indices from the current selection with spatially-similar alternatives.
    Priority: same FOV + same grid bin; then same FOV any bin; then any FOV any bin.
    """
    cfg = session["cfg"]
    rows, cols = session["rows"], session["cols"]
    selected = session["selected_idxs"]
    cands = session["candidates"]
    backups = session["backups"]

    # track already used indices to avoid duplicates
    used = set(selected)

    def _pull_replacement(orig_idx):
        orig = cands[orig_idx]
        fov_key = orig["key"]
        bin_rc = _candidate_bin(orig, rows, cols)

        # 1) same FOV + same bin
        key_bin = (fov_key, bin_rc)
        if key_bin in backups:
            while backups[key_bin]:
                j = backups[key_bin].pop(0)
                if j not in used:
                    return j
        # 2) same FOV, any bin
        fov_bins = [kb for kb in backups.keys() if kb[0]==fov_key]
        for kb in fov_bins:
            while backups[kb]:
                j = backups[kb].pop(0)
                if j not in used:
                    return j
        # 3) any FOV, any bin
        for kb in list(backups.keys()):
            while backups[kb]:
                j = backups[kb].pop(0)
                if j not in used:
                    return j
        return None

    replaced = []
    for b in bad_indices:
        idx0 = b - 1
        if idx0 < 0 or idx0 >= len(selected):
            continue
        orig_cand_idx = selected[idx0]
        repl = _pull_replacement(orig_cand_idx)
        if repl is not None:
            selected[idx0] = repl
            used.add(repl)
            replaced.append((b, repl))
        else:
            if cfg.verbose:
                print(f"No replacement available for index {b}")

    if cfg.verbose and replaced:
        print("Replaced:", ", ".join([f"{i}->{cands[j]['cid']}" for (i,j) in replaced]))
    return session


def finalize_review(session: dict):
    """
    Write per-cell panels, a gallery, and the two CSVs using the *current* selection.
    Mirrors process_sample() outputs but uses the reviewed selection.
    """
    sample_id = session["sample"]
    cfg = session["cfg"]
    mpp = session["mpp"]
    out_root = Path(cfg.outputs_dir) / sample_id
    out_qc = out_root / "qc"
    ensure_dir(out_root); ensure_dir(out_qc)

    selected = [session["candidates"][i] for i in session["selected_idxs"]]

    # Panels + gallery
    if cfg.save_individual_qc_panels:
        for idx, item in enumerate(selected, start=1):
            fig = make_qc_panel(
                item["cid"], item["phase"], item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
                item["green"], item["red"], cfg, overlap_img=item.get("overlap")
            )
            fig.savefig(out_qc / f"{sample_id}_cell_{idx:03d}_panel.png", dpi=150)
            plt.close(fig)

    cols = cfg.qc_panels_per_row
    rows = math.ceil(len(selected) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols*3, rows*3))
    axes = np.atleast_1d(axes).ravel()

    for i, item in enumerate(selected):
        fig_i = make_qc_panel(
            item["cid"], item["phase"], item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
            item["green"], item["red"], cfg, overlap_img=item.get("overlap")
        )
        fig_i.canvas.draw()
        w, h = fig_i.canvas.get_width_height()
        buf = np.frombuffer(fig_i.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
        axes[i].imshow(buf[..., :3]); axes[i].axis("off"); axes[i].set_title(f"{i+1}")
        plt.close(fig_i)

    for j in range(i+1, len(axes)): axes[j].axis("off")
    fig.tight_layout()
    if cfg.save_combined_qc_panel:
        fig.savefig(out_qc / f"{sample_id}_selected_cells_QC.png", dpi=150)
    plt.close(fig)


    # --- NEW: per-FOV selection maps (red boxes around chosen cells) ---
    if cfg.save_selection_maps:
        by_fov = {}
        for rank, item in enumerate(selected, start=1):
            k = item["key"]
            by_fov.setdefault(k, []).append({
                "rank": rank,
                "cid": item["cid"],
                "lbl_cells": item["lbl_cells"],
                "phase": item["phase"],
                "green": item["green"],
                "red": item["red"],
                "overlap": item.get("overlap"),
            })
        for k, items in by_fov.items():
            base = _compose_overlap_like(items[0]["phase"], items[0]["green"], items[0]["red"], items[0]["overlap"])
            _save_fov_selection_map(sample_id, k, base, items, out_qc, cfg)
    # -------------------------------------------------------------------

    # --- Selection manifest (reproduce exact picks later) ---
    if getattr(cfg, "save_selection_manifest", False):
        manifest_path = out_qc / f"{sample_id}_selection_manifest.json"
        _write_selection_manifest(selected, manifest_path)
        if cfg.verbose:
            print("Wrote:", manifest_path)
    # --------------------------------------------------------

    # Measurements + CSVs
    rows_cell, rows_nuc = [], []
    for idx, item in enumerate(selected, start=1):
        key = item["key"]; cid = item["cid"]
        dfc, dfn = measure_cell_and_nucleus(
            item["lbl_cells"], item["lbl_nuclei"], item["nuc_classes"],
            [cid], mpp, cfg, green_img=item["green"], red_img=item["red"]
        )
        for d in dfc.to_dict("records"):
            d.update({"root": key.root, "sample": key.sample, "well": key.well, "site": key.site, "timecode": key.timecode, "selected_rank": idx})
            rows_cell.append(d)
        for d in dfn.to_dict("records"):
            d.update({"root": key.root, "sample": key.sample, "well": key.well, "site": key.site, "timecode": key.timecode, "selected_rank": idx})
            rows_nuc.append(d)

    df_cell = pd.DataFrame(rows_cell); df_nuc = pd.DataFrame(rows_nuc)
    ts = __import__("datetime").datetime.now().strftime("%Y%m%d_%H%M%S")
    cell_csv = out_root / f"cell_morphology_{sample_id}_{ts}.csv"
    nuc_csv = out_root / f"nuclear_morphology_{sample_id}_{ts}.csv"
    df_cell.to_csv(cell_csv, index=False); df_nuc.to_csv(nuc_csv, index=False)
    if cfg.verbose:
        print("Wrote:", cell_csv); print("Wrote:", nuc_csv); print("QC dir:", out_qc)

    return {"mpp": mpp, "cell_csv": str(cell_csv), "nuc_csv": str(nuc_csv), "qc_dir": str(out_qc)}