"""
neighbor_correlation_scatter.py
--------------------------------
Plots scatter of sprinkled-cell expression vs. nearest-neighbor mean expression
for the Aldob gene at 12 hr and 72 hr post-transplantation (Figure 3 panels D/E).

Call plot_neighbor_correlation_scatter() to generate both PDFs.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy.spatial.distance as sci_dist
import numpy.ma as ma
from itertools import compress
from numpy.core.numeric import argwhere

try:
    from utils.constant import ZONATION_PLASTICITY_PLOTS_FOLDER_PATH, SPRINKLING_NOV23_BASE_PATH
    _DEFAULT_OUT_DIR  = ZONATION_PLASTICITY_PLOTS_FOLDER_PATH
    _DEFAULT_BASE     = SPRINKLING_NOV23_BASE_PATH
except ImportError:
    _DEFAULT_BASE = '/Users/yaelheyman/Library/CloudStorage/GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/sprinkling_nov_23'
    _DEFAULT_OUT_DIR  = os.path.join(os.getcwd(), 'paper', 'graphs', 'zonation_plasticity_plots')

_GFP_THRESH = 2
_RADIUS     = 400
_GENE       = 'Aldob'
_FIGURE_SIZE = [2, 2]
_DPI        = 300
_FONT_SIZE  = 8


def _load_roi(base, timepoint, roi):
    cbg   = os.path.join(base, timepoint, roi, 'output', 'cell_by_gene', 'cell_by_gene.csv')
    attrs = os.path.join(base, timepoint, roi, 'output', 'attributes', 'cell_attributes.csv')
    metadata = pd.read_csv(cbg, index_col=0)
    coords   = pd.read_csv(attrs, index_col=0).reindex(columns=['center_x', 'center_y']).to_numpy()
    gfp_idx  = list(compress(range(len(metadata['GFP'])), metadata['GFP'] > _GFP_THRESH))
    return metadata, coords, gfp_idx


def _find_neighbors(coords, gfp_idx, radius):
    dist_mat = sci_dist.squareform(sci_dist.pdist(coords))
    neighbors_idx = argwhere(dist_mat[:, gfp_idx] < radius)
    idx_sort = np.argsort(neighbors_idx[:, 1])
    s = neighbors_idx[idx_sort]
    _, idx_start, _ = np.unique(s[:, 1], return_counts=True, return_index=True)
    return [grp.tolist() for grp in np.split(s[:, 0], idx_start[1:])]


def _plot_scatter(gene, metadata, gfp_id, coords, out_path, radius=_RADIUS):
    self_exp     = metadata.iloc[gfp_id][gene].to_numpy()
    neighbor_list = _find_neighbors(coords, gfp_id, radius)
    nn_mean      = [np.mean(metadata.iloc[np.setdiff1d(nb, gfp_id)][gene]) for nb in neighbor_list]
    corr         = ma.corrcoef(ma.masked_invalid(self_exp), ma.masked_invalid(np.array(nn_mean)), rowvar=False)[0, 1]

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42})
        fig, ax = plt.subplots(figsize=_FIGURE_SIZE)
        ax.scatter(nn_mean, self_exp, c='g', s=0.8)
        ax.set_xlabel('nearest neighbors mean expression of ' + gene, fontsize=_FONT_SIZE)
        ax.set_ylabel('inserted cell expression of  ' + gene, fontsize=_FONT_SIZE)
        ax.text(0.1, 0.9, 'correlation with neighbors  ' + str(np.round(corr, decimals=3)),
                transform=ax.transAxes, color='black', fontsize=_FONT_SIZE, ha='left', va='center')
        ax.tick_params(labelsize=_FONT_SIZE)
        ax.set_xlim(0, 150)
        ax.set_ylim(0, 150)
        ax.set_box_aspect(1)          # square plotting box
        ax.grid(True, which='major')
        plt.tight_layout()
        plt.savefig(out_path, dpi=_DPI, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        print(f"  Saved: {out_path}")
    finally:
        matplotlib.rcParams.update(_saved_rc)


def plot_neighbor_correlation_scatter(output_dir=None, base_path=None):
    """Generate neighbor correlation scatter plots for 12 hr and 72 hr (panels D and E)."""
    _out  = output_dir or _DEFAULT_OUT_DIR
    _base = base_path  or _DEFAULT_BASE
    os.makedirs(_out, exist_ok=True)

    print("  Loading 12hr data…")
    meta12, c12, gfp12 = _load_roi(_base, '12hr', 'roi2')
    _plot_scatter(_GENE, meta12, gfp12, c12,
                  os.path.join(_out, 'neighbor_corr_scatter_12hr.pdf'))

    print("  Loading 72hr data…")
    meta72, c72, gfp72 = _load_roi(_base, '72hr', 'roi1')
    _plot_scatter(_GENE, meta72, gfp72, c72,
                  os.path.join(_out, 'neighbor_corr_scatter_72hr.pdf'))


if __name__ == '__main__':
    plot_neighbor_correlation_scatter()
