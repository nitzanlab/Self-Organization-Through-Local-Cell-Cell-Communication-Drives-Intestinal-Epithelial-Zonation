#!/usr/bin/env python3
"""
scale_invariance_figure.py  –  200 mm wide, Illustrator-compatible SVG

Top row (same height, total = 200 mm):
  A – tissue image (same FOV for all), all 3 timepoints' sub-ROI polygons
      (dashed, 2× thick, no lumenal/basal borders)
  B – same FOV as A, rainbow ring-assignment illustration overlaid on
      day5 sub-ROI (shows how the distance-from-basal profiles are computed)
  C – section-width bar chart (slightly boosted width)

Legend strip (Arial 12 pt) between top row and profiles.

Bottom full-width:
  D – raw CPM expression profiles (2-row × 4-gene panel)

Usage: python scale_invariance_figure.py [rep]   (default: rep2)
"""

import base64, io, json, re, sys
import numpy as np
import pandas as pd
import cv2
from pathlib import Path
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as mcm

# ── MODULE-LEVEL CONSTANTS ────────────────────────────────────────────────────
TIMEPOINTS = ["day0", "day3", "day5"]
N_ITER     = 20
GSM_BY_REP = {"rep1": "GSM9134407", "rep2": "GSM9134408",
              "rep3": "GSM9134409", "rep4": "GSM9134410"}
COLORS   = {"day0": "#8C1A6A", "day3": "#E8941A", "day5": "#1B7FA3"}
TP_LABEL = {"day0": "Pre-irradiation",
            "day3": "3 Days post-irradiation",
            "day5": "5 Days post-irradiation"}

# layout constants (mm)
FIG_W       = 205.0   # +5 mm so Panel D x-axis labels have room
L_MARGIN    = 5.0
R_MARGIN    = 5.0
GAP         = 5.0       # gap between top-row panels
LEGEND_H    = 7.0    # one horizontal row at 8pt
PROF_GAP    = 4.0
B_MARGIN    = 5.0
PAD_PX      = 200       # hires-px padding around combined sub-ROI bbox
SCALEBAR_UM = 200
RAINBOW_OPACITY = 0.75  # opacity of solid ring colors over tissue
B_PAD_PX    = 100      # padding around day5 subroi for Panel B

# SVG font sizes (mm = pt / 72 * 25.4)
PT14_MM = 14 / 72 * 25.4
PT12_MM = 12 / 72 * 25.4
PT8_MM  = 8  / 72 * 25.4

# Panel subheadings
SUBHEAD = {
    "A": ("Intestinal sections, 0–5 days", "post-irradiation"),
    "B": ("Inflation analysis",),
    "C": "Width post-irradiation",
}

# ── HELPERS ───────────────────────────────────────────────────────────────────
def f(v): return f"{v:.4f}"

def load_poly(p):
    df = pd.read_csv(p)
    return df[["axis-0", "axis-1"]].to_numpy()   # [row, col]

def poly_to_svgpts(rc, c0, r0, cw, ch, dw, dh):
    xs = (rc[:, 1] - c0) / cw * dw
    ys = (rc[:, 0] - r0) / ch * dh
    return " ".join(f"{x:.3f},{y:.3f}" for x, y in zip(xs, ys))

def rc_to_svg(rc, c0, r0, cw, ch, dw, dh):
    return float((rc[1]-c0)/cw*dw), float((rc[0]-r0)/ch*dh)

def encode_png_b64(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG", compress_level=1)
    return base64.b64encode(buf.getvalue()).decode("ascii")

def extract_svg(src):
    text = Path(src).read_text("utf-8") if (
        isinstance(src, Path) or not str(src).lstrip().startswith("<")
    ) else src
    vb   = re.search(r'viewBox=["\']([^"\']+)["\']', text)
    body = re.search(r'<svg[^>]*>(.*)</svg>', text, re.DOTALL)
    return (vb.group(1) if vb else "0 0 100 100",
            body.group(1) if body else "")

def label_pos_rc(sr_rc, border_rc, frac=0.20):
    t = cKDTree(border_rc)
    d, _ = t.query(sr_rc)
    n = max(1, int(len(sr_rc) * frac))
    return sr_rc[np.argsort(d)[:n]].mean(axis=0)

def panel_subheading(lines_tuple, panel_w, font_mm=None, pad=1.0):
    """White semi-transparent box with 1 or 2 lines of text at top-centre of a panel."""
    if font_mm is None:
        font_mm = PT8_MM
    lh   = font_mm * 1.4
    n    = len(lines_tuple)
    bh   = lh * n + pad * 2
    cx   = panel_w / 2
    out  = (f'<rect x="0" y="0" width="{f(panel_w)}" height="{f(bh)}" '
            f'fill="white" opacity="0.72"/>')
    for i, txt in enumerate(lines_tuple):
        y = pad + lh * (i + 0.82)
        out += (f'<text x="{f(cx)}" y="{f(y)}" text-anchor="middle" '
                f'font-size="{f(font_mm)}" font-family="Arial">{txt}</text>')
    return out

def pill_label(cx, cy, text, col, font_mm, txt_col="white"):
    lw = len(text) * font_mm * 0.58 + 2.0
    return (
        f'<rect x="{f(cx-lw/2)}" y="{f(cy-font_mm*0.72)}" '
        f'width="{f(lw)}" height="{f(font_mm*1.05)}" '
        f'rx="{f(font_mm*0.3)}" fill="{col}" opacity="0.85"/>'
        f'<text x="{f(cx)}" y="{f(cy)}" text-anchor="middle" '
        f'fill="{txt_col}" font-size="{f(font_mm)}" font-weight="bold">{text}</text>'
    )

def pill_label_2l(cx, cy, line1, line2, col, font_mm, txt_col="white"):
    """Two-line pill label centered at (cx, cy)."""
    lw   = max(len(line1), len(line2)) * font_mm * 0.58 + 2.0
    lh   = font_mm * 1.35
    rh   = lh * 2 + font_mm * 0.2
    ry   = cy - lh - font_mm * 0.1
    y1   = cy - font_mm * 0.1
    y2   = cy + lh - font_mm * 0.1
    return (
        f'<rect x="{f(cx-lw/2)}" y="{f(ry)}" '
        f'width="{f(lw)}" height="{f(rh)}" '
        f'rx="{f(font_mm*0.3)}" fill="{col}" opacity="0.85"/>'
        f'<text text-anchor="middle" fill="{txt_col}" '
        f'font-size="{f(font_mm)}" font-weight="bold">'
        f'<tspan x="{f(cx)}" y="{f(y1)}">{line1}</tspan>'
        f'<tspan x="{f(cx)}" y="{f(y2)}">{line2}</tspan>'
        f'</text>'
    )

TP_LABEL_2L = {
    "day0": ("Pre-", "irradiation"),
    "day3": ("3 Days", "post-irradiation"),
    "day5": ("5 Days", "post-irradiation"),
}


# ── MAIN FIGURE ASSEMBLY ──────────────────────────────────────────────────────

def create_figure_2(rep="rep2", data_root=None, out_dir=None):
    """
    Assemble the full Figure 2 SVG (panels A–D).

    Parameters
    ----------
    rep       : replicate identifier (default "rep2")
    data_root : path to Visium dataset root (defaults to VISIUM_DATA_ROOT in constant.py)
    out_dir   : where to save output SVGs (defaults to SCALE_INVARIANCE_PLOTS_FOLDER_PATH)
    """
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["font.family"]  = "Arial"
    matplotlib.rcParams["font.size"]    = 8

    # ── resolve paths ─────────────────────────────────────────────────────────
    if data_root:
        _data_root = Path(data_root)
    else:
        try:
            from utils.constant import VISIUM_DATA_ROOT
            _data_root = Path(VISIUM_DATA_ROOT)
        except ImportError:
            _data_root = Path(
                "/Users/yaelheyman/Library/CloudStorage/"
                "GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/"
                "mouse visium data/versi data set"
            )

    if out_dir:
        _out_dir = Path(out_dir)
    else:
        try:
            from utils.constant import SCALE_INVARIANCE_PLOTS_FOLDER_PATH
            _out_dir = Path(SCALE_INVARIANCE_PLOTS_FOLDER_PATH)
        except ImportError:
            _out_dir = _data_root / rep / "subroi_profiles"

    # polygon CSVs live in the repo's extractedData folder
    mask_dir = Path(__file__).parent.parent / "extractedData" / "visium_masks"

    RAW_DIR       = _data_root / "GSE303705_RAW"
    GSM           = GSM_BY_REP[rep]
    HIRES_IMG     = RAW_DIR / f"{GSM}_{rep}_tissue_hires_image.png"
    DAY_B         = "day5"
    MASK_DIR_B    = _data_root / rep / DAY_B / "infl_out" / "inflation_num_iterations_10"
    PROFILES_SVG  = _out_dir / f"{rep}_panel_c_3row_subroi_ribbon.svg"
    PROFILES_SVG_RAW = _out_dir / f"{rep}_panel_d_raw_2row_subroi_ribbon.svg"
    OUT_SVG       = _out_dir / f"{rep}_zonation_figure.svg"

    # ── SCALE FACTORS ─────────────────────────────────────────────────────────
    with open(RAW_DIR / f"{GSM}_{rep}_scalefactors_json.json") as fh:
        sf = json.load(fh)
    UM_PER_PX = sf["microns_per_pixel"] / sf["tissue_hires_scalef"]
    print(f"µm/hires-px: {UM_PER_PX:.4f}")

    # ── LOAD POLYGON DATA & WIDTH STATS ───────────────────────────────────────
    print("Loading polygons and width stats …")
    poly_data, width_stats = {}, {}
    for day in TIMEPOINTS:
        poly_data[day] = {
            "lum":    load_poly(mask_dir / f"{rep}_{day}.csv"),
            "basal":  load_poly(mask_dir / f"{rep}_{day}_basal.csv"),
            "subroi": load_poly(mask_dir / f"{rep}_{day}_subroi.csv"),
        }
        ws = pd.read_csv(_data_root / rep / day / "infl_out"
                         / f"inflation_num_iterations_{N_ITER}" / "subroi_width_stats.csv")
        width_stats[day] = {
            "mean_um": float(ws["mean_px"].iloc[0]) * UM_PER_PX,
            "std_um":  float(ws["std_px"].iloc[0])  * UM_PER_PX,
        }
        print(f"  {day}: {width_stats[day]['mean_um']:.1f} ± {width_stats[day]['std_um']:.1f} µm")

    # ── SHARED CROP: bounding box of all sub-ROIs ──────────────────────────────
    print("\nCropping …")
    img_arr      = np.array(Image.open(HIRES_IMG))
    IMG_H, IMG_W = img_arr.shape[:2]

    all_sub = np.vstack([poly_data[d]["subroi"] for d in TIMEPOINTS])
    cr0 = max(0,       int(all_sub[:,0].min()) - PAD_PX)
    cr1 = min(IMG_H-1, int(all_sub[:,0].max()) + PAD_PX)
    cc0 = max(0,       int(all_sub[:,1].min()) - PAD_PX)
    cc1 = min(IMG_W-1, int(all_sub[:,1].max()) + PAD_PX)
    crop       = img_arr[cr0:cr1+1, cc0:cc1+1]
    CROP_H, CROP_W = crop.shape[:2]
    A_RATIO    = CROP_W / CROP_H

    # "lumenal"/"basal" annotation positions (day0)
    sr0        = poly_data["day0"]["subroi"]
    lum_lbl_rc = label_pos_rc(sr0, poly_data["day0"]["lum"])
    bas_lbl_rc = label_pos_rc(sr0, poly_data["day0"]["basal"])

    print(f"  Tissue crop: {CROP_W}×{CROP_H} px  (ratio {A_RATIO:.3f})")

    # ── PANEL B: day5 sub-ROI, solid rainbow from per-ring masks ──────────────
    print(f"Building Panel B (solid rainbow, {DAY_B} sub-ROI) …")

    srB_rc  = poly_data[DAY_B]["subroi"]
    b_cr0   = max(0,       int(srB_rc[:,0].min()) - B_PAD_PX)
    b_cr1   = min(IMG_H-1, int(srB_rc[:,0].max()) + B_PAD_PX)
    b_cc0   = max(0,       int(srB_rc[:,1].min()) - B_PAD_PX)
    b_cc1   = min(IMG_W-1, int(srB_rc[:,1].max()) + B_PAD_PX)
    b_crop  = img_arr[b_cr0:b_cr1+1, b_cc0:b_cc1+1]
    B_CROP_H, B_CROP_W = b_crop.shape[:2]
    B_RATIO = B_CROP_W / B_CROP_H

    # subroi mask in crop space
    mask_b_pil = Image.new("L", (B_CROP_W, B_CROP_H), 0)
    draw_b     = ImageDraw.Draw(mask_b_pil)
    sr_px_b    = [(int(c - b_cc0), int(r - b_cr0)) for r, c in srB_rc]
    draw_b.polygon(sr_px_b, fill=255)
    mask_b = np.array(mask_b_pil) > 0

    # paint each inflation ring with a solid rainbow color
    ring_mask_paths = sorted(MASK_DIR_B.glob("inflated_mask_*.png"),
                             key=lambda p: int(p.stem.split("_")[-1]))
    N_RINGS = len(ring_mask_paths)
    cmap_b  = mcm.get_cmap("rainbow")
    panel_b = b_crop[:, :, :3].copy()
    prev_full = np.zeros((IMG_H, IMG_W), dtype=bool)
    for i, mp in enumerate(ring_mask_paths):
        full_mask  = np.array(Image.open(mp)) > 128
        ring_only  = full_mask & ~prev_full
        prev_full  = full_mask
        ring_crop  = ring_only[b_cr0:b_cr1+1, b_cc0:b_cc1+1]
        apply_mask = ring_crop & mask_b
        if apply_mask.sum() > 0:
            rv, gv, bv, _ = cmap_b(i / max(N_RINGS - 1, 1))
            solid = np.array([rv*255, gv*255, bv*255], dtype=np.float32)
            tissue_px = b_crop[:, :, :3].astype(np.float32)[apply_mask]
            panel_b[apply_mask] = np.clip(
                RAINBOW_OPACITY * solid + (1 - RAINBOW_OPACITY) * tissue_px, 0, 255
            ).astype(np.uint8)

    # lumenal/basal label positions in Panel B
    srB_in_b = srB_rc[
        (srB_rc[:,0] >= b_cr0) & (srB_rc[:,0] <= b_cr1) &
        (srB_rc[:,1] >= b_cc0) & (srB_rc[:,1] <= b_cc1)
    ]
    if len(srB_in_b) >= 3:
        lum_lbl_b = label_pos_rc(srB_in_b, poly_data[DAY_B]["lum"])
        bas_lbl_b = label_pos_rc(srB_in_b, poly_data[DAY_B]["basal"])
    else:
        lum_lbl_b = np.array([(b_cr0+b_cr1)//2, b_cc0+100])
        bas_lbl_b = np.array([(b_cr0+b_cr1)//2, b_cc1-100])

    print(f"  {DAY_B} crop: {B_CROP_W}×{B_CROP_H} px  (ratio {B_RATIO:.3f}), {N_RINGS} rings")

    # ── TOP-ROW GEOMETRY ──────────────────────────────────────────────────────
    AVAIL   = FIG_W - L_MARGIN - R_MARGIN - 2 * GAP
    H_top   = AVAIL / (2.5 * A_RATIO + B_RATIO)
    A_W     = H_top * A_RATIO
    B_W     = H_top * B_RATIO
    C_W     = 1.5 * A_W
    C_H     = 1.5 * H_top * 0.7
    A_X   = L_MARGIN
    B_X   = A_X + A_W + GAP
    C_X   = B_X + B_W + GAP
    print(f"Top row: H_AB={H_top:.1f} mm  A={A_W:.1f}  B={B_W:.1f}  C={C_W:.1f}×{C_H:.1f} mm")
    print(f"  Tissue DPI: {CROP_W/A_W*25.4:.0f}")

    # ── PANEL C: bar chart at C_W × C_H ───────────────────────────────────────
    print("Building bar chart …")
    fig_bar, ax_bar = plt.subplots(figsize=(C_W / 25.4, C_H / 25.4))
    xs    = np.arange(len(TIMEPOINTS))
    means = [width_stats[d]["mean_um"] for d in TIMEPOINTS]
    stds  = [width_stats[d]["std_um"]  for d in TIMEPOINTS]
    ax_bar.bar(xs, means, yerr=stds, color=[COLORS[d] for d in TIMEPOINTS],
               width=0.7, capsize=5,
               error_kw={"linewidth": 1.2, "capthick": 1.2}, edgecolor="none")
    ax_bar.set_xticks(xs)
    ax_bar.set_xticklabels([TP_LABEL[d] for d in TIMEPOINTS], rotation=40, ha="right")
    ax_bar.set_title(SUBHEAD["C"])
    ax_bar.set_ylabel("Section width (µm)")
    ax_bar.set_ylim(0, max(means) * 1.45)
    ax_bar.set_xlim(-0.65, len(TIMEPOINTS) - 0.35)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    fig_bar.subplots_adjust(bottom=0.32, top=0.95, left=0.18, right=0.97)
    buf = io.BytesIO()
    fig_bar.savefig(buf, format="svg")
    plt.close(fig_bar)
    bar_vb, bar_body = extract_svg(buf.getvalue().decode("utf-8"))

    # ── PROFILES SVG ──────────────────────────────────────────────────────────
    print("Loading profiles SVG …")
    prof_vb, prof_body = extract_svg(PROFILES_SVG)
    prof_body = re.sub(
        r'<text\b[^>]*>[^<]*(top:\s*mean CPM|ring\s*width)[^<]*</text>',
        '', prof_body, flags=re.IGNORECASE
    )
    pv = [float(x) for x in prof_vb.split()]
    PROF_W_MM = FIG_W - L_MARGIN - R_MARGIN
    PROF_H_MM = PROF_W_MM * (pv[3] / pv[2])

    # ── DOCUMENT DIMENSIONS ───────────────────────────────────────────────────
    T_MARGIN = L_MARGIN
    LEGEND_Y = T_MARGIN + max(H_top, C_H)
    PROF_Y   = LEGEND_Y + LEGEND_H + PROF_GAP
    DOC_H    = PROF_Y + PROF_H_MM + B_MARGIN
    DOC_W    = FIG_W
    print(f"Document: {DOC_W:.0f} × {DOC_H:.0f} mm")

    # ── SVG ASSEMBLY ──────────────────────────────────────────────────────────
    lines = [
        '<?xml version="1.0" encoding="utf-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{f(DOC_W)}mm" height="{f(DOC_H)}mm" '
        f'viewBox="0 0 {f(DOC_W)} {f(DOC_H)}" font-family="Arial">',
        f'<rect x="0" y="0" width="{f(DOC_W)}" height="{f(DOC_H)}" fill="white"/>',
        "<defs>",
    ]

    # clip-path defs (used in Panel A only for the sub-ROI fills)
    for day in TIMEPOINTS:
        pts = poly_to_svgpts(poly_data[day]["subroi"],
                             cc0, cr0, CROP_W, CROP_H, A_W, H_top)
        lines.append(f'  <clipPath id="clip_{day}"><polygon points="{pts}"/></clipPath>')
    lines.append("</defs>")

    # ── Panel A ───────────────────────────────────────────────────────────────
    lines.append(f'<g id="panel_A" transform="translate({f(A_X)},{f(T_MARGIN)})">')
    lines.append(f'<text x="-1" y="-1.5" font-size="{f(PT14_MM)}" font-weight="bold">A</text>')

    b64A = encode_png_b64(crop)
    lines.append(
        f'<image x="0" y="0" width="{f(A_W)}" height="{f(H_top)}" '
        f'preserveAspectRatio="none" xlink:href="data:image/png;base64,{b64A}"/>'
    )

    # sub-ROI polygons: dashed, 2× stroke, no lumenal/basal
    for day in TIMEPOINTS:
        col = COLORS[day]
        pts = poly_to_svgpts(poly_data[day]["subroi"],
                             cc0, cr0, CROP_W, CROP_H, A_W, H_top)
        lines.append(
            f'<polygon points="{pts}" fill="{col}" fill-opacity="0.07" '
            f'stroke="{col}" stroke-width="1.2" stroke-dasharray="2.5,1.2" '
            f'stroke-linejoin="round" stroke-opacity="0.9"/>'
        )
        # timepoint label at sub-ROI centroid
        cen = poly_data[day]["subroi"].mean(axis=0)
        cx, cy = rc_to_svg(cen, cc0, cr0, CROP_W, CROP_H, A_W, H_top)
        lines.append(pill_label_2l(cx, cy, *TP_LABEL_2L[day], col, PT8_MM))

    # subheading overlay (Panel A)
    lines.append(panel_subheading(SUBHEAD["A"], A_W))

    # scale bar
    sb_mm = (SCALEBAR_UM / UM_PER_PX) / CROP_W * A_W
    sb_x0, sb_y = A_W - sb_mm - 2.0, H_top - 3.2
    lines.append(
        f'<line x1="{f(sb_x0)}" y1="{f(sb_y)}" x2="{f(sb_x0+sb_mm)}" y2="{f(sb_y)}" '
        f'stroke="white" stroke-width="0.55" stroke-linecap="round"/>'
        f'<text x="{f(sb_x0+sb_mm/2)}" y="{f(sb_y+2.3)}" text-anchor="middle" '
        f'fill="white" font-size="{f(PT8_MM)}">{SCALEBAR_UM} µm</text>'
    )

    lines.append("</g>")  # end Panel A

    # ── Panel B: same FOV, rainbow overlay ────────────────────────────────────
    lines.append(f'<g id="panel_B" transform="translate({f(B_X)},{f(T_MARGIN)})">')
    lines.append(f'<text x="-1" y="-1.5" font-size="{f(PT14_MM)}" font-weight="bold">B</text>')

    b64B = encode_png_b64(panel_b)
    lines.append(
        f'<image x="0" y="0" width="{f(B_W)}" height="{f(H_top)}" '
        f'preserveAspectRatio="none" xlink:href="data:image/png;base64,{b64B}"/>'
    )

    # day5 sub-ROI outline on top
    pts_dB = poly_to_svgpts(poly_data[DAY_B]["subroi"],
                             b_cc0, b_cr0, B_CROP_W, B_CROP_H, B_W, H_top)
    lines.append(
        f'<polygon points="{pts_dB}" fill="none" '
        f'stroke="{COLORS[DAY_B]}" stroke-width="1.2" stroke-dasharray="2.5,1.2" '
        f'stroke-linejoin="round" stroke-opacity="0.9"/>'
    )

    # lumenal / basal annotation labels in Panel B
    for lbl_rc, label in [(lum_lbl_b, "lumenal"), (bas_lbl_b, "basal")]:
        tx, ty = rc_to_svg(lbl_rc, b_cc0, b_cr0, B_CROP_W, B_CROP_H, B_W, H_top)
        lw = len(label) * PT8_MM * 0.58 + 1.5
        lines.append(
            f'<rect x="{f(tx-lw/2)}" y="{f(ty-PT8_MM*0.75)}" '
            f'width="{f(lw)}" height="{f(PT8_MM*1.1)}" '
            f'rx="{f(PT8_MM*0.3)}" fill="black" opacity="0.65"/>'
            f'<text x="{f(tx)}" y="{f(ty)}" text-anchor="middle" fill="white" '
            f'font-size="{f(PT8_MM)}">{label}</text>'
        )

    # subheading overlay (Panel B)
    lines.append(panel_subheading(SUBHEAD["B"], B_W))

    # scale bar
    sb_mm_b = (SCALEBAR_UM / UM_PER_PX) / B_CROP_W * B_W
    sb_x0b, sb_yb = B_W - sb_mm_b - 2.0, H_top - 3.2
    lines.append(
        f'<line x1="{f(sb_x0b)}" y1="{f(sb_yb)}" x2="{f(sb_x0b+sb_mm_b)}" y2="{f(sb_yb)}" '
        f'stroke="white" stroke-width="0.55" stroke-linecap="round"/>'
        f'<text x="{f(sb_x0b+sb_mm_b/2)}" y="{f(sb_yb+2.3)}" text-anchor="middle" '
        f'fill="white" font-size="{f(PT8_MM)}">{SCALEBAR_UM} µm</text>'
    )

    lines.append("</g>")  # end Panel B

    # ── Panel C: bar chart ────────────────────────────────────────────────────
    lines.append(
        f'<text x="{f(C_X-1)}" y="{f(T_MARGIN-1.5)}" '
        f'font-size="{f(PT14_MM)}" font-weight="bold">C</text>'
    )
    lines.append(
        f'<svg x="{f(C_X)}" y="{f(T_MARGIN)}" '
        f'width="{f(C_W)}" height="{f(C_H)}" '
        f'viewBox="{bar_vb}" preserveAspectRatio="none">'
        + bar_body + '</svg>'
    )

    # ── Legend strip (horizontal) ─────────────────────────────────────────────
    lines.append(f'<g id="legend" transform="translate({f(L_MARGIN)},{f(LEGEND_Y+2)})">')
    tx  = 0.0
    my  = LEGEND_H * 0.5
    for day in TIMEPOINTS:
        col = COLORS[day]
        txt = TP_LABEL[day]
        lines.append(
            f'<line x1="{f(tx)}" y1="{f(my)}" x2="{f(tx+9)}" y2="{f(my)}" '
            f'stroke="{col}" stroke-width="1.5"/>'
            f'<text x="{f(tx+11)}" y="{f(my + PT8_MM*0.38)}" '
            f'font-size="{f(PT8_MM)}">{txt}</text>'
        )
        tx += 9 + 2 + len(txt) * PT8_MM * 0.60 + 7
    lines.append("</g>")

    # ── WRITE — two variants (y-normalised Panel D and raw CPM Panel D) ───────
    _out_dir.mkdir(parents=True, exist_ok=True)

    _variants = [
        (PROFILES_SVG,     OUT_SVG,                                          "y-norm"),
        (PROFILES_SVG_RAW, _out_dir / "figure_2_scale_invariance.svg",         "raw CPM"),
    ]

    for _prof_path, _out_path, _label in _variants:
        if not _prof_path.exists():
            print(f"\nSkipping {_label} variant — {_prof_path.name} not found")
            continue
        _vb, _body = extract_svg(_prof_path)
        _pv = [float(x) for x in _vb.split()]
        _prof_h = PROF_W_MM * (_pv[3] / _pv[2])
        _doc_h  = PROF_Y + _prof_h + B_MARGIN

        _lines = (
            ['<?xml version="1.0" encoding="utf-8"?>',
             f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
             f'width="{f(DOC_W)}mm" height="{f(_doc_h)}mm" '
             f'viewBox="0 0 {f(DOC_W)} {f(_doc_h)}" font-family="Arial">',
             f'<rect x="0" y="0" width="{f(DOC_W)}" height="{f(_doc_h)}" fill="white"/>']
            + lines[3:]   # <defs> … legend </g>
            + [f'<text x="{f(L_MARGIN-1)}" y="{f(PROF_Y-1.5)}" '
               f'font-size="{f(PT14_MM)}" font-weight="bold">D</text>',
               f'<svg x="{f(L_MARGIN)}" y="{f(PROF_Y)}" '
               f'width="{f(PROF_W_MM)}" height="{f(_prof_h)}" '
               f'viewBox="{_vb}" preserveAspectRatio="xMidYMid meet">'
               + _body + '</svg>',
               '</svg>']
        )
        _out_path.write_text("\n".join(_lines), encoding="utf-8")
        sz = _out_path.stat().st_size / 1e6
        print(f"\nSaved ({_label}): {_out_path.name}  ({sz:.1f} MB)")
        print(f"Figure: {DOC_W:.0f} × {_doc_h:.0f} mm")


if __name__ == "__main__":
    _rep = sys.argv[1] if len(sys.argv) > 1 else "rep2"
    create_figure_2(rep=_rep)
