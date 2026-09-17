"""
scale_invariance_subroi_width.py
--------------------------------
For each timepoint in a rep, filters the existing width_profile.csv to only
points whose basal endpoint falls inside the sub-ROI polygon, then:
  1. Saves subroi_width_stats.csv   → ITER_DIR
  2. Saves subroi_width_overlay.png → ITER_DIR
     (tissue image + width lines + sub-ROI polygon boundary labelled on top)

Usage
-----
    python scale_invariance_subroi_width.py          # uses config below
    python scale_invariance_subroi_width.py rep4
"""

import sys
import os
import numpy as np
import pandas as pd
import cv2
import imageio.v3 as iio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================
REP        = sys.argv[1] if len(sys.argv) > 1 else "rep2"
TIMEPOINTS = ["day0", "day3", "day5"]
N_ITER     = 20

try:
    from utils.constant import VISIUM_DATA_ROOT
    DATA_ROOT = Path(VISIUM_DATA_ROOT)
except ImportError:
    DATA_ROOT = Path(
        "/Users/yaelheyman/Library/CloudStorage/"
        "GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/"
        "mouse visium data/versi data set"
    )

MASK_DIR = Path(__file__).parent.parent / "extractedData" / "visium_masks"
RAW_DIR = DATA_ROOT / "GSE303705_RAW"
GSM_BY_REP = {
    "rep1": "GSM9134407",
    "rep2": "GSM9134408",
    "rep3": "GSM9134409",
    "rep4": "GSM9134410",
}
# ============================================================


def poly2mask(vertices_xy, H, W):
    """vertices_xy: (N,2) in (x,y)=(col,row). Returns uint8 mask (H,W)."""
    pts = np.asarray(vertices_xy, dtype=np.float32)
    mask = np.zeros((H, W), dtype=np.uint8)
    cv2.fillPoly(mask, [pts.astype(np.int32)], 1)
    return mask


def process_timepoint(rep, day, n_iter, data_root=None, mask_dir=None):
    _data_root = Path(data_root) if data_root else DATA_ROOT
    _mask_dir  = Path(mask_dir)  if mask_dir  else MASK_DIR
    roi_dir  = _data_root / rep / day
    iter_dir = roi_dir / "infl_out" / f"inflation_num_iterations_{n_iter}"
    sample   = f"{GSM_BY_REP[rep]}_{rep}"
    hires_image = _data_root / "GSE303705_RAW" / f"{sample}_tissue_hires_image.png"

    width_csv  = iter_dir / "width_profile.csv"
    subroi_csv = _mask_dir / f"{rep}_{day}_subroi.csv"

    # ── checks ──────────────────────────────────────────────────────────────
    for p, name in [(width_csv, "width_profile.csv"), (subroi_csv, "subroi polygon CSV")]:
        if not p.exists():
            print(f"  [{day}] SKIP — {name} not found: {p}")
            return

    print(f"  [{day}] loading data…")

    # ── load width profile ───────────────────────────────────────────────────
    width_df = pd.read_csv(width_csv)

    # ── load sub-ROI polygon → mask ──────────────────────────────────────────
    img = iio.imread(str(hires_image))
    if img.ndim == 2:
        img = np.stack([img]*3, axis=-1)
    elif img.shape[2] == 4:
        img = img[:,:,:3]
    H, W = img.shape[:2]

    df_poly = pd.read_csv(subroi_csv)
    vertices_xy = df_poly[["axis-1", "axis-0"]].to_numpy()
    subroi_mask = poly2mask(vertices_xy, H, W)

    # ── filter width points inside sub-ROI (by basal endpoint) ──────────────
    bx = np.clip(width_df["basal_x"].to_numpy().astype(int), 0, W-1)
    by = np.clip(width_df["basal_y"].to_numpy().astype(int), 0, H-1)
    in_subroi = subroi_mask[by, bx] > 0

    width_df["in_subroi"] = in_subroi

    # ── compute stats for included + in-subroi points ────────────────────────
    sub_inc = width_df[in_subroi & width_df["included"]]
    n_total   = int(in_subroi.sum())
    n_included = len(sub_inc)

    if n_included == 0:
        print(f"  [{day}] WARNING: no included points in sub-ROI")
        stats = {"mean_px": float("nan"), "std_px": float("nan"),
                 "median_px": float("nan"), "n_total": n_total, "n_included": 0}
    else:
        stats = {
            "mean_px":   float(sub_inc["width_px"].mean()),
            "std_px":    float(sub_inc["width_px"].std()),
            "median_px": float(sub_inc["width_px"].median()),
            "n_total":   n_total,
            "n_included": n_included,
        }

    print(f"  [{day}] sub-ROI width: {stats['mean_px']:.1f} ± {stats['std_px']:.1f} px "
          f"(n={n_included} points)")

    # ── save stats CSV ────────────────────────────────────────────────────────
    stats_path = iter_dir / "subroi_width_stats.csv"
    pd.DataFrame([stats]).to_csv(stats_path, index=False)
    print(f"  [{day}] saved: {stats_path.name}")

    # ── build overlay image ───────────────────────────────────────────────────
    img_f = img.astype(np.float32) / 255.0

    # color scale based on sub-ROI included widths
    w_px = width_df["width_px"].to_numpy()
    inc  = width_df["included"].to_numpy()
    w_inc_sub = w_px[in_subroi & inc]
    vmin = w_inc_sub.min() if len(w_inc_sub) else w_px.min()
    vmax = w_inc_sub.max() if len(w_inc_sub) else w_px.max()
    if vmin == vmax:
        vmin, vmax = vmin * 0.9, vmax * 1.1 + 1
    norm    = Normalize(vmin=vmin, vmax=vmax)
    cmap_fn = cm.get_cmap("plasma")

    # subsample lines
    n_lines = 400
    idx = np.round(np.linspace(0, len(width_df)-1, min(n_lines, len(width_df)))).astype(int)
    df_sub = width_df.iloc[idx]
    lw = max(1.0, H / 800)

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # --- Left panel: tissue + width lines + sub-ROI polygon ---
    ax = axes[0]
    ax.imshow(img_f, origin="upper")

    for _, row in df_sub.iterrows():
        x0, y0 = row["basal_x"], row["basal_y"]
        x1, y1 = row["lum_x"],   row["lum_y"]
        is_in  = subroi_mask[int(np.clip(y0,0,H-1)), int(np.clip(x0,0,W-1))] > 0
        if row["included"] and is_in:
            color = cmap_fn(norm(row["width_px"]))
            alpha = 0.85
        elif not row["included"] or not is_in:
            color = (0.5, 0.5, 0.5, 0.3)
            alpha = 0.3
        ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw, alpha=alpha)

    # sub-ROI polygon boundary
    poly_closed = np.vstack([vertices_xy, vertices_xy[0]])
    ax.plot(poly_closed[:, 0], poly_closed[:, 1],
            color="yellow", linewidth=max(2, lw*2),
            linestyle="--", label="Sub-ROI", zorder=5)

    # colorbar
    sm = cm.ScalarMappable(cmap="plasma", norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, label="Width (px)")

    ax.set_title(f"{rep} {day} — local width + sub-ROI", fontsize=13, fontweight="bold")
    ax.legend(loc="upper right", fontsize=11, facecolor="black",
              labelcolor="white", framealpha=0.7)
    ax.axis("off")

    # --- Right panel: histogram ---
    ax2 = axes[1]
    all_inc = w_px[inc]
    if len(all_inc):
        ax2.hist(all_inc, bins=30, color="steelblue", alpha=0.5,
                 label=f"Full ROI (n={len(all_inc)})", density=True)
    if len(w_inc_sub):
        ax2.hist(w_inc_sub, bins=30, color="orange", alpha=0.7,
                 label=f"Sub-ROI (n={len(w_inc_sub)})", density=True)
        ax2.axvline(stats["mean_px"], color="darkorange", lw=2, ls="--",
                    label=f"Sub-ROI mean: {stats['mean_px']:.1f} px")

    ax2.set_xlabel("Width (px)", fontsize=11)
    ax2.set_ylabel("Density", fontsize=11)
    ax2.set_title("Width distribution", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    out_img = iter_dir / "subroi_width_overlay.png"
    plt.savefig(out_img, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [{day}] saved: {out_img.name}")


# ============================================================
if __name__ == "__main__":
    print(f"Sub-ROI width analysis — {REP}")
    for day in TIMEPOINTS:
        process_timepoint(REP, day, N_ITER)
    print("\nDone.")
