"""
transplant_neighborhoods.py
----------------------------
Plots transplanted-cell neighborhoods for Figure 3 panels B and C.

Format: DAPI (blue) + GFP (green) composite background, Aldob transcript
dots in red, 10 µm scale bar. No cell-outline overlays.

Output per cell: <name>_COMBO.png  (raster background)
                 <name>_COMBO.svg  (vector transcript dots + scale bar,
                                    Illustrator-editable, links the PNG)

Call plot_transplant_neighborhoods() to generate all panels.
"""

import os
import io
import numpy as np
import pandas as pd
import tifffile as tiff
import scipy.spatial.distance as sci_dist
from PIL import Image

try:
    from utils.constant import (
        ZONATION_PLASTICITY_PLOTS_FOLDER_PATH,
        SPRINKLING_NOV23_BASE_PATH,
        SPRINKLING_BG_72,
        SPRINKLING_BG_12,
    )
    _DEFAULT_OUT_DIR = ZONATION_PLASTICITY_PLOTS_FOLDER_PATH
    _DEFAULT_BASE    = SPRINKLING_NOV23_BASE_PATH
    _DEFAULT_BG_72   = SPRINKLING_BG_72
    _DEFAULT_BG_12   = SPRINKLING_BG_12
except ImportError:
    _BASE = '/Users/yaelheyman/Library/CloudStorage/GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/sprinkling_nov_23'
    _DEFAULT_BASE    = _BASE
    _DEFAULT_BG_72   = '/Users/yaelheyman/RajLab Dropbox/Yael Heyman/shared_yael/sg/sprinkling nov 23/72 hr/roi_1/hyb_background_aligned.tiff'
    _DEFAULT_BG_12   = '/Users/yaelheyman/RajLab Dropbox/Yael Heyman/shared_yael/sg/sprinkling nov 23/12 hr/roi_2/hyb_background_aligned.tiff'
    _DEFAULT_OUT_DIR = os.path.join(os.getcwd(), 'paper', 'graphs', 'zonation_plasticity_plots')

# ── chosen GFP cell IDs ────────────────────────────────────────────────────
_GFP_CELLS_72 = [5259, 7872, 9502]   # panel C (72 hr)
_GFP_CELLS_12 = [7328, 8172, 1501]   # panel B (12 hr)

# ── imaging parameters ─────────────────────────────────────────────────────
_DAPI_CH       = 3
_GFP_CH        = 2
_GENE          = 'Aldob'
_NN_RADIUS     = 350
_PIXEL_UM      = 0.10711
_SCALEBAR_UM   = 10
_SCALEBAR_PX   = int(round(_SCALEBAR_UM / _PIXEL_UM))
_Q_LO, _Q_HI  = 1.0, 98.5
_GAMMA         = 0.85


# ── internal helpers ───────────────────────────────────────────────────────
def _norm(gray):
    g = gray.astype(np.float32)
    lo, hi = np.percentile(g, [_Q_LO, _Q_HI])
    if hi <= lo:
        return np.zeros_like(gray, dtype=np.uint8)
    g = np.clip((g - lo) / (hi - lo), 0, 1)
    return (np.power(g, 1.0 / _GAMMA) * 255).astype(np.uint8)


def _find_neighbors(coords, gfp_idx, radius):
    dist = sci_dist.squareform(sci_dist.pdist(coords))
    pairs = np.argwhere(dist[:, gfp_idx] < radius)
    if pairs.size == 0:
        return [[] for _ in gfp_idx]
    isort = np.argsort(pairs[:, 1])
    s = pairs[isort]
    _, starts, _ = np.unique(s[:, 1], return_index=True, return_counts=True)
    return [g.tolist() for g in np.split(s[:, 0], starts[1:])]


def _get_transcripts(transcripts, gene):
    if gene in transcripts.index:
        return transcripts.loc[[gene], ['x', 'y']].to_numpy()
    for col in [c for c in transcripts.columns if c.lower() in ('name', 'gene', 'target')]:
        mask = transcripts[col].astype(str).str.lower() == gene.lower()
        if mask.any():
            return transcripts.loc[mask, ['x', 'y']].to_numpy()
    raise KeyError(f"Gene '{gene}' not found in transcripts table")


def _png_bytes(rgb):
    buf = io.BytesIO()
    Image.fromarray(rgb).save(buf, format='PNG', compress_level=6)
    return buf.getvalue()


def _make_svg(png_path, W, H, spots, x0, y0, add_gene_label=False):
    href = os.path.basename(png_path)
    margin = int(min(W, H) * 0.04)
    bar_h  = max(2, int(H * 0.008))
    bx, by = margin, H - margin - bar_h
    fs     = max(8, int(min(W, H) * 0.04))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        f'<image x="0" y="0" width="{W}" height="{H}" '
        f'preserveAspectRatio="none" href="{href}" xlink:href="{href}"/>',
    ]
    for x, y in spots:
        px, py = x - x0, y - y0
        if 0 <= px < W and 0 <= py < H:
            parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="2.5" '
                         f'fill="#FF0000" fill-opacity="1.0"/>')
    parts.append(
        f'<rect x="{bx}" y="{by}" width="{_SCALEBAR_PX}" height="{bar_h}" fill="#FFFFFF"/>'
    )
    parts.append(
        f'<text x="{bx + _SCALEBAR_PX / 2:.1f}" y="{by - 4}" '
        f'fill="#FFFFFF" font-family="Arial" font-size="{fs}" text-anchor="middle">'
        f'{_SCALEBAR_UM} µm</text>'
    )
    if add_gene_label:
        parts.append(
            f'<text x="{margin}" y="{margin + fs}" '
            f'fill="#FF0000" font-family="Arial" font-style="italic" font-size="{fs}">'
            f'{_GENE}</text>'
        )
    parts.append('</svg>')
    return '\n'.join(parts)


def _plot_one(gfp_cell_id, image, coords, transcripts, out_base, add_gene_label=False):
    neighbor_list = _find_neighbors(coords, [gfp_cell_id], _NN_RADIUS)
    neigh_ids = np.array(neighbor_list[0]) if neighbor_list[0] else np.array([gfp_cell_id])
    neigh_coords = coords[neigh_ids]

    min_x, max_x = neigh_coords[:, 0].min(), neigh_coords[:, 0].max()
    min_y, max_y = neigh_coords[:, 1].min(), neigh_coords[:, 1].max()
    pad_x = (max_x - min_x) / 4.0
    pad_y = (max_y - min_y) / 4.0
    x0 = int(max(0, np.floor(min_x - pad_x)))
    x1 = int(min(image.shape[2], np.ceil(max_x + pad_x)))
    y0 = int(max(0, np.floor(min_y - pad_y)))
    y1 = int(min(image.shape[1], np.ceil(max_y + pad_y)))

    dapi_8 = _norm(image[_DAPI_CH, y0:y1, x0:x1])
    gfp_8  = _norm(image[_GFP_CH,  y0:y1, x0:x1])

    H, W = dapi_8.shape
    combo = np.zeros((H, W, 3), dtype=np.uint8)
    combo[..., 1] = gfp_8
    combo[..., 2] = dapi_8

    spots = _get_transcripts(transcripts, _GENE)
    in_crop = ((spots[:, 0] > x0) & (spots[:, 0] < x1) &
               (spots[:, 1] > y0) & (spots[:, 1] < y1))

    png_path = f'{out_base}_COMBO.png'
    Image.fromarray(combo).save(png_path)

    svg = _make_svg(png_path, W, H, spots[in_crop], x0, y0, add_gene_label=add_gene_label)
    with open(f'{out_base}_COMBO.svg', 'w', encoding='utf-8') as f:
        f.write(svg)

    print(f'  gfp{gfp_cell_id}: {W}×{H} px  |  {in_crop.sum()} {_GENE} dots')


# ── public entry point ─────────────────────────────────────────────────────
def plot_transplant_neighborhoods(output_dir=None, base_path=None,
                                  bg_72=None, bg_12=None):
    """Generate neighborhood images for Figure 3 panels B (12 hr) and C (72 hr)."""
    _out  = output_dir or _DEFAULT_OUT_DIR
    _base = base_path  or _DEFAULT_BASE
    _bg72 = bg_72      or _DEFAULT_BG_72
    _bg12 = bg_12      or _DEFAULT_BG_12
    os.makedirs(_out, exist_ok=True)

    print('  Loading 72hr data …')
    img72 = tiff.imread(_bg72)
    coords72 = pd.read_csv(
        f'{_base}/72hr/roi1/output/attributes/cell_attributes.csv', index_col=0
    )[['center_x', 'center_y']].to_numpy()
    tx72 = pd.read_csv(f'{_base}/72hr/roi1/output/transcripts/transcripts.csv', index_col=0)

    print('  Plotting 72hr neighborhoods (panel C) …')
    for i, cell_id in enumerate(_GFP_CELLS_72):
        out_base = os.path.join(_out, f'neighborhood_72hr_gfp{cell_id}')
        _plot_one(cell_id, img72, coords72, tx72, out_base, add_gene_label=(i == 0))

    print('  Loading 12hr data …')
    img12 = tiff.imread(_bg12)
    coords12 = pd.read_csv(
        f'{_base}/12hr/roi2/output/attributes/cell_attributes.csv', index_col=0
    )[['center_x', 'center_y']].to_numpy()
    tx12 = pd.read_csv(f'{_base}/12hr/roi2/output/transcripts/transcripts.csv', index_col=0)

    print('  Plotting 12hr neighborhoods (panel B) …')
    for i, cell_id in enumerate(_GFP_CELLS_12):
        out_base = os.path.join(_out, f'neighborhood_12hr_gfp{cell_id}')
        _plot_one(cell_id, img12, coords12, tx12, out_base, add_gene_label=(i == 0))


if __name__ == '__main__':
    plot_transplant_neighborhoods()
