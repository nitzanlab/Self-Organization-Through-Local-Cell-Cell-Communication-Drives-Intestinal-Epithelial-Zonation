"""
cluster_monolayer.py
--------------------
Plots monolayer cells colored by unsupervised cluster identity (Leiden).
Uses pre-computed cluster annotations from the full-image clustering run.

Output: monolayer_cluster_identity.pdf → AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH

Call plot_monolayer_cluster_identity() to generate the figure.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    from utils.constant import AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, MONOLAYER_RAW_DATA_ROOT
    _DEFAULT_OUT_DIR  = AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH
    _DEFAULT_RAW_ROOT = MONOLAYER_RAW_DATA_ROOT
except ImportError:
    _DEFAULT_RAW_ROOT = (
        '/Users/yaelheyman/Library/CloudStorage/'
        'GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/pasadena_run_no_gel'
    )
    _DEFAULT_OUT_DIR = os.path.join(
        os.getcwd(), 'paper', 'graphs', 'autonomous_zonation_figure_plots'
    )

_PIXEL_SIZE_NM = 107.11
_CROP          = (9219, 17780, 14269, 20986)   # min_x, max_x, min_y, max_y

_CLUSTER_COLOR_DICT = {
    '0': (.4, .4, .4),
    '1': (.4, .4, .4),
    '2': (.4, .4, .4),
    '4': (.4, .4, .4),
    '3': 'lime',
    '5': 'magenta',
    '6': 'cyan',
}

_LEGEND_PATCHES = [
    mpatches.Patch(color=(.4, .4, .4), label='Regenerative'),
    mpatches.Patch(color='lime',       label='Enterocyte'),
    mpatches.Patch(color='magenta',    label='Secretory'),
    mpatches.Patch(color='cyan',       label='Progenitor'),
]


def plot_monolayer_cluster_identity(output_dir=None, raw_data_root=None):
    """Generate monolayer_cluster_identity.pdf for the autonomous zonation panel."""
    from SGanalysis.SGobject import SGobject

    _out_dir  = output_dir    or _DEFAULT_OUT_DIR
    _raw_root = raw_data_root or _DEFAULT_RAW_ROOT

    nuclei_seg  = os.path.join(_raw_root, 'attributes', 'segmentations', 'segmentation_20230809_nuclei.tiff')
    cluster_csv = os.path.join(_raw_root, 'umap_figures_8', 'cell_by_gene_cluster_annotations.csv')
    out_path    = os.path.join(_out_dir, 'monolayer_cluster_identity.pdf')

    os.makedirs(_out_dir, exist_ok=True)

    # load nucleus geometries (no transcript loading needed)
    sg = SGobject()
    print("  Loading segmentation…")
    sg.mask_to_objects(nuclei_seg)
    print(f"  Total nuclei: {len(sg.gdf)}")

    # map cluster IDs from pre-computed CSV
    cl = pd.read_csv(cluster_csv)
    cl['object_id']  = pd.to_numeric(cl['object_id'],  errors='coerce').astype('Int64')
    cl['cluster_id'] = cl['cluster_id'].astype(str)
    cl_map = dict(zip(cl['object_id'], cl['cluster_id']))
    sg.gdf['object_id'] = pd.to_numeric(sg.gdf['object_id'], errors='coerce').astype('Int64')
    sg.gdf['leiden']    = sg.gdf['object_id'].map(cl_map)

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({
            'figure.titlesize': 8, 'axes.titlesize': 8, 'axes.labelsize': 8,
            'ytick.labelsize': 8,  'xtick.labelsize': 8, 'legend.fontsize': 6,
            'font.family': 'DejaVu Sans', 'savefig.dpi': 300, 'pdf.fonttype': 42,
        })

        fig, ax = plt.subplots(figsize=(3, 3))

        # plot all cells (full monolayer), crop axes afterward
        for cid, color in _CLUSTER_COLOR_DICT.items():
            subset = sg.gdf[sg.gdf['leiden'] == cid]
            if not subset.empty:
                subset.plot(ax=ax, color=color, linewidth=0, edgecolor='none', alpha=0.9)

        min_x, max_x, min_y, max_y = _CROP
        ax.set_xlim([min_x, max_x])
        ax.set_ylim([min_y, max_y])
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(''); ax.set_ylabel('')

        # scale bar
        scale_bar_px = int(100 * 1000.0 / _PIXEL_SIZE_NM)
        bar_height   = scale_bar_px / 4.665
        ax.add_patch(plt.Rectangle(
            (min_x + 200, max_y - 200 - bar_height),
            scale_bar_px, bar_height, color='black', lw=0,
        ))

        plt.legend(handles=_LEGEND_PATCHES, title="Cell Types",
                   loc='upper left', prop={'size': 6}, title_fontsize=6)
        for sp in ax.spines.values():
            sp.set_visible(False)

        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out_path}")
    finally:
        matplotlib.rcParams.update(_saved_rc)


if __name__ == '__main__':
    plot_monolayer_cluster_identity()
