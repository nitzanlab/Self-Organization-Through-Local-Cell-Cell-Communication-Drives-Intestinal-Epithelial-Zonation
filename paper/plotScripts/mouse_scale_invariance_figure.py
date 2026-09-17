"""
Figure 2 panels A-C: mouse irradiation (rep2, day0 / day3 / day5).

Panel A  three timepoint sub-ROIs cropped from the rep2 hires H&E, oriented
         luminal (cyan, top) -> basal (red, bottom) with the border arcs.
Panel B  mean section width per timepoint (+/- STD).
Panel C  per-gene profiles, one line per timepoint, vs distance from basal border.

Converted from mouse_final_figure.ipynb. Paths come from utils.constant; the ring
helpers are vendored in utils.ring_profiles so this has no external dependency.
Panel D is produced separately by scale_invariance_profiles.plot_scale_invariance_figures.
"""
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
from PIL import Image

from utils.ring_profiles import ring_profiles, profiles_to_dataframe, smooth1d

Image.MAX_IMAGE_PIXELS = None

try:
    from utils.constant import VISIUM_DATA_ROOT, SCALE_INVARIANCE_PLOTS_FOLDER_PATH
    _DEFAULT_DATA_ROOT = VISIUM_DATA_ROOT
    _DEFAULT_OUT_DIR = SCALE_INVARIANCE_PLOTS_FOLDER_PATH
except ImportError:                                           # standalone use
    _DEFAULT_DATA_ROOT = ('/Users/yaelheyman/Library/CloudStorage/'
                          'GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/'
                          'mouse visium data/versi data set')
    _DEFAULT_OUT_DIR = os.path.join(os.getcwd(), 'paper', 'graphs', 'scale_invariance_plots')

_MASKS_DIR = Path(__file__).resolve().parents[1] / 'extractedData' / 'visium_masks'

REP = 'rep2'
N_ITER = 20
KERNEL_SIZE = 20
UM_PER_PX = 0.3443494972902308 / 0.114027254      # microns_per_pixel / tissue_hires_scalef
RING_WIDTH_UM = (KERNEL_SIZE / 2) * UM_PER_PX
TIMEPOINTS = ['day0', 'day3', 'day5']
TP_LABEL = {'day0': 'day0 (pre-irradiation)', 'day3': 'day3 post-irradiation',
            'day5': 'day5 post-irradiation'}
TCOL = {'day0': '#8C1A6A', 'day3': '#E8941A', 'day5': '#1B7FA3'}
SMOOTH_SIGMA = 0
_HIRES_REL = os.path.join('GSE303705_RAW', 'GSM9134408_rep2_tissue_hires_image.png')

GENE_CATEGORY = {
    'Ada': 'villus top', 'Apoa4': 'villus top', 'Alpi': 'villus mid',
    'Sis': 'villus bottom', 'Slc5a1': 'villus bottom', 'Plac8': 'villus bottom',
    'Olfm4': 'stem cell', 'Lgr5': 'stem cell', 'Mki67': 'progenitor',
    'Lyz1': 'paneth', 'Muc2': 'goblet', 'Col1a1': 'fibrosis', 'Acta2': 'fibrosis',
    'Epha2': 'ephrin', 'Epha1': 'ephrin', 'Efna1': 'ephrin',
}
CANONICAL_GENE_ORDER = ['ADA', 'APOA4', 'ALPI', 'SI', 'SLC5A1', 'PLAC8', 'OLFM4',
                        'LGR5', 'MKI67', 'LYZ', 'DEFA5', 'MUC2', 'COL1A1', 'ACTA2',
                        'EPHA2', 'EPHA1', 'EFNA1']


def _gene_rank(g):
    """Shared gene order across Crohn's / Celiac / mouse; orthologs Sis->SI, Lyz1->LYZ."""
    k = {'SIS': 'SI', 'LYZ1': 'LYZ'}.get(g.upper(), g.upper())
    return CANONICAL_GENE_ORDER.index(k) if k in CANONICAL_GENE_ORDER else 99


def _line(p):
    return pd.read_csv(p)[['axis-1', 'axis-0']].to_numpy()        # (x, y)


def _load_tp(tp, data_root):
    it = (Path(data_root) / REP / tp / 'infl_out' /
          f'inflation_num_iterations_{N_ITER}' / 'ring_subset_bin_by_gene_subroi')
    X = sp.load_npz(it / 'counts_genes_by_bins.npz')
    bc = pd.read_csv(it / 'barcodes.tsv', header=None, sep='\t')[0].astype(str).tolist()
    feat = pd.read_csv(it / 'features.tsv', sep='\t')
    ring = pd.read_csv(it / 'barcode_to_ring_aligned.csv')
    return X, bc, feat, ring


def _rotate_bound(im, ang):
    import cv2
    h, w = im.shape[:2]
    M = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), ang, 1.0)
    cos, sin = abs(M[0, 0]), abs(M[0, 1])
    nW, nH = int(h * sin + w * cos), int(h * cos + w * sin)
    M[0, 2] += nW / 2.0 - w / 2.0
    M[1, 2] += nH / 2.0 - h / 2.0
    return cv2.warpAffine(im, M, (nW, nH), flags=cv2.INTER_LINEAR,
                          borderValue=(255, 255, 255)), M


def _panelA_crop(tp, img_full, pad=70):
    """Crop rep2 hires to a timepoint sub-ROI and orient luminal(top) -> basal(bottom)."""
    import cv2
    H, W = img_full.shape[:2]
    sub = _line(_MASKS_DIR / f'{REP}_{tp}_subroi.csv')
    lum = _line(_MASKS_DIR / f'{REP}_{tp}.csv')            # luminal (inner)
    bas = _line(_MASKS_DIR / f'{REP}_{tp}_basal.csv')      # basal (outer)
    x0, x1 = int(sub[:, 0].min() - pad), int(sub[:, 0].max() + pad)
    y0, y1 = int(sub[:, 1].min() - pad), int(sub[:, 1].max() + pad)
    x0, y0 = max(0, x0), max(0, y0)
    x1, y1 = min(W, x1), min(H, y1)
    crop = img_full[y0:y1, x0:x1].copy()

    def arc(P):
        m = (P[:, 0] >= x0) & (P[:, 0] <= x1) & (P[:, 1] >= y0) & (P[:, 1] <= y1)
        return P[m]

    la, ba = arc(lum), arc(bas)
    d = ba.mean(0) - la.mean(0)
    cx, cy = crop.shape[1] / 2.0, crop.shape[0] / 2.0
    best = (-1e9, 0)
    for phi in range(0, 360):
        R = cv2.getRotationMatrix2D((cx, cy), phi, 1.0)[:, :2]
        dv = R @ d
        down = dv[1] / (np.hypot(*dv) + 1e-9)
        if down > best[0]:
            best = (down, phi)
    raw_o, M = _rotate_bound(crop, best[1])

    def tf(P):
        q = P - [x0, y0]
        return (M @ np.c_[q, np.ones(len(q))].T).T

    return raw_o, tf(la), tf(ba)


def plot_mouse_scale_invariance_figure(output_dir=None, data_root=None):
    """Build Figure 2 panels A-C. Returns the list of files written (empty if none)."""
    out_dir = Path(output_dir or _DEFAULT_OUT_DIR)
    root = Path(data_root or _DEFAULT_DATA_ROOT)
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = sorted(GENE_CATEGORY, key=_gene_rank)
    profiles, width = {}, {}
    for tp in TIMEPOINTS:
        X, bc, feat, ring = _load_tp(tp, root)
        ur, m, s, n, found, na = ring_profiles(X, bc, feat, ring, genes,
                                               exclude_last_ring=True, min_bin_counts=1)
        profiles[tp] = profiles_to_dataframe(ur, m, s, n, found, RING_WIDTH_UM, ns_all=na)
        w = pd.read_csv(root / REP / tp / 'infl_out' /
                        f'inflation_num_iterations_{N_ITER}' / 'subroi_width_stats.csv')
        width[tp] = (float(w['mean_px'][0]) * UM_PER_PX, float(w['std_px'][0]) * UM_PER_PX)
        print(f'  {tp}: width {width[tp][0]:.0f} +/- {width[tp][1]:.0f} um')

    genes = [g for g in genes if g in set(profiles[TIMEPOINTS[0]]['gene'])]
    if not genes:
        warnings.warn('mouse scale-invariance figure: no genes resolved in the ring '
                      'profiles - panels A-C not generated.', RuntimeWarning)
        return []
    print('  genes:', genes)

    img_full = np.array(Image.open(root / _HIRES_REL))[:, :, :3]

    ncolC = 5
    nrowC = int(np.ceil(len(genes) / ncolC))
    fig = plt.figure(figsize=(18, 11 + 2.7 * nrowC))
    outer = fig.add_gridspec(2, 1, height_ratios=[11.0, 2.7 * nrowC], hspace=0.14)

    top = outer[0].subgridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.14)
    gsA = top[0].subgridspec(3, 1, hspace=0.15)
    for i, tp in enumerate(TIMEPOINTS):
        ax = fig.add_subplot(gsA[i])
        raw_o, la, ba = _panelA_crop(tp, img_full)
        ax.imshow(raw_o)
        ax.plot(la[:, 0], la[:, 1], color='#00c8ff', lw=2, solid_capstyle='round')
        ax.plot(ba[:, 0], ba[:, 1], color='#ff2828', lw=2, solid_capstyle='round')
        lm, bm = la.mean(0), ba.mean(0)
        ax.annotate('luminal', (lm[0], lm[1]), color='white', fontsize=7, fontweight='bold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.2', fc='#00c8ff', ec='none', alpha=0.85))
        ax.annotate('basal', (bm[0], bm[1]), color='white', fontsize=7, fontweight='bold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.2', fc='#ff2828', ec='none', alpha=0.85))
        ax.set_title(TP_LABEL[tp], fontsize=10, fontweight='bold')
        ax.axis('off')
        ax.set_xlim(0, raw_o.shape[1])
        ax.set_ylim(raw_o.shape[0], 0)
        if i == 0:
            ax.text(-0.03, 1.08, 'A', transform=ax.transAxes, fontsize=17, fontweight='bold')

    axB = fig.add_subplot(top[1])
    xs = np.arange(len(TIMEPOINTS))
    axB.bar(xs, [width[t][0] for t in TIMEPOINTS], yerr=[width[t][1] for t in TIMEPOINTS],
            capsize=7, color=[TCOL[t] for t in TIMEPOINTS], width=0.6,
            error_kw=dict(lw=1.4, capthick=1.4))
    axB.set_xticks(xs)
    axB.set_xticklabels(TIMEPOINTS, fontsize=11)
    axB.set_ylabel('section width (um)\n(lumen -> basal)', fontsize=10)
    axB.set_title('Mean section width (+/- STD)', fontsize=11, fontweight='bold')
    axB.spines[['top', 'right']].set_visible(False)
    axB.text(-0.16, 1.04, 'B', transform=axB.transAxes, fontsize=17, fontweight='bold')

    gsC = outer[1].subgridspec(nrowC, ncolC, hspace=0.85, wspace=0.35)
    for k, gene in enumerate(genes):
        ax = fig.add_subplot(gsC[k // ncolC, k % ncolC])
        for tp in TIMEPOINTS:
            sub = profiles[tp][profiles[tp]['gene'] == gene].sort_values('ring_order_lumen_to_basal')
            if sub.empty:
                continue
            x = sub['dist_from_basal_um'].to_numpy()
            y = smooth1d(sub['mean_cpm'].to_numpy(), SMOOTH_SIGMA)
            e = smooth1d((sub['std_cpm'] / np.sqrt(sub['n_bins'].clip(lower=1))).to_numpy(),
                         SMOOTH_SIGMA)
            ax.plot(x, y, color=TCOL[tp], lw=2, marker='o', ms=3, label=tp)
            ax.fill_between(x, y - e, y + e, color=TCOL[tp], alpha=0.15, lw=0)
        ax.set_title(gene, fontsize=9, fontweight='bold')
        ax.text(0.5, 1.28, GENE_CATEGORY[gene], transform=ax.transAxes, ha='center',
                va='bottom', fontsize=8, style='italic', color='#444')
        ax.spines[['top', 'right']].set_visible(False)
        if k == 0:
            ax.legend(fontsize=7, frameon=False)
            ax.text(-0.42, 1.55, 'C', transform=ax.transAxes, fontsize=17, fontweight='bold')
        if k // ncolC == nrowC - 1:
            ax.set_xlabel('distance from basal border (um)', fontsize=7)
        if k % ncolC == 0:
            ax.set_ylabel('mean expression', fontsize=8)
    for k in range(len(genes), nrowC * ncolC):
        fig.add_subplot(gsC[k // ncolC, k % ncolC]).axis('off')

    fig.suptitle('Mouse irradiation (rep2): crypt-villus zonation across timepoints. '
                 'Absolute counts.', fontsize=13, fontweight='bold', y=0.995)

    written = []
    for ext in ('png', 'pdf', 'svg'):
        p = out_dir / f'mouse_scale_invariance_panels_ABC.{ext}'
        fig.savefig(p, dpi=200, bbox_inches='tight')
        written.append(str(p))
    plt.close(fig)
    print('  saved:', out_dir / 'mouse_scale_invariance_panels_ABC.{png,pdf,svg}')
    return written


if __name__ == '__main__':
    plot_mouse_scale_invariance_figure()
