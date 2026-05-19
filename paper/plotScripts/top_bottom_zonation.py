"""
top_bottom_zonation.py
----------------------
Plots the relative expression of villus-top genes spatially over the monolayer,
using segmented nuclei (dilated to cells), transcript assignments, and cluster
annotations to retain only enterocyte-like cells (cluster 3).

Output: monolayer_top_bottom_villus_expression.pdf → AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH

Call plot_top_bottom_villus_expression() to generate the figure.
"""

import os
import sys
import threading
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

_PIXEL_SIZE_NM      = 107.11
_SCALE_BAR_UM       = 100
_SCALE_BAR_MARGIN   = 200
_CROP               = (9219, 17780, 14269, 20986)   # min_x, max_x, min_y, max_y
_KEEP_CLUSTERS      = {'3'}                          # only enterocyte-like cells
_DILATION_RADIUS    = 30


def _heartbeat(msg="Working", interval_s=10):
    stop = threading.Event()
    def _run():
        sys.stdout.write(msg + "\n"); sys.stdout.flush()
        i = 0
        while not stop.wait(interval_s):
            i += 1; sys.stdout.write(f".  ({i*interval_s}s)\n"); sys.stdout.flush()
    threading.Thread(target=_run, daemon=True).start()
    return stop


def _build_cell_gene_table(nuclei_seg, transcripts, crop, cbg_cache):
    from SGanalysis.SGobject import SGobject
    min_x, max_x, min_y, max_y = crop
    sg = SGobject()
    print(f"  Loading segmentation…")
    sg.mask_to_objects(nuclei_seg)
    print(f"  Loading transcripts…")
    sg.load_points(transcripts)
    if hasattr(sg, 'points_gdf') and sg.points_gdf is not None:
        n_before = len(sg.points_gdf)
        sg.points_gdf = sg.points_gdf[
            (sg.points_gdf['x'] >= min_x) & (sg.points_gdf['x'] <= max_x) &
            (sg.points_gdf['y'] >= min_y) & (sg.points_gdf['y'] <= max_y)
        ].copy()
        print(f"  Points filtered to ROI: {n_before} → {len(sg.points_gdf)}")
    print("  Dilating nuclei to approximate cells…")
    sg.dilate_objects(radius=_DILATION_RADIUS, identifier='nucleus',
                      output_name='dilated_nucleus', set_geometry=True)
    if os.path.exists(cbg_cache):
        print(f"  Loading cached cell-by-gene: {cbg_cache}")
        cell_gene_df = pd.read_csv(cbg_cache)
    else:
        print("  Building cell-by-gene table (slow)…")
        tick = _heartbeat("  Assigning transcripts to cells…")
        try:
            sg.create_cell_gene_table(index_col='object_id')
        finally:
            tick.set()
        cell_gene_df = sg.get_cell_gene_table_df().reset_index()
        cell_gene_df.to_csv(cbg_cache, index=False)
        print(f"  Cached to {cbg_cache}")
    return sg, cell_gene_df


def _merge_annotations(cell_gene_df, annotated_cbg):
    print(f"  Merging top/bottom annotations…")
    tp_bottom = pd.read_csv(annotated_cbg)
    if 'Unnamed: 0' in tp_bottom.columns and 'name' not in tp_bottom.columns:
        tp_bottom = tp_bottom.rename(columns={'Unnamed: 0': 'name'})
    cell_gene_df = cell_gene_df.copy()
    cell_gene_df['object_id'] = pd.to_numeric(cell_gene_df['object_id'], errors='coerce')
    cell_gene_df['name_adjusted'] = cell_gene_df['object_id'] - 1
    tp_bottom['name'] = pd.to_numeric(tp_bottom['name'], errors='coerce')
    merged = cell_gene_df.merge(
        tp_bottom, left_on='name_adjusted', right_on='name',
        how='left', suffixes=('_gene', '_tp')
    )
    if 'signal' not in merged.columns:
        raise KeyError("Expected a 'signal' column in the annotations file")
    return merged


def _add_scale_bar(ax, min_x, max_y, length_um, margin_px, px_size_nm):
    length_px = int(length_um * 1000.0 / px_size_nm)
    height_px = max(2, int(length_px / 20))
    bar_left   = min_x + margin_px
    bar_bottom = max_y - margin_px - height_px
    ax.add_patch(plt.Rectangle((bar_left, bar_bottom), length_px, height_px, color='black', lw=0))
    ax.text(bar_left + length_px / 2, bar_bottom - 0.75 * height_px,
            f"{length_um} µm", ha='center', va='bottom', fontsize=7)


def _plot_and_save(sg, merged, out_path, cluster_csv, crop):
    from shapely.geometry import box as shapely_box
    from matplotlib import cm, colors as mcolors

    print("  Preparing spatial plot…")
    sg.gdf['object_id'] = pd.to_numeric(sg.gdf['object_id'], errors='coerce').astype('Int64')
    signal_map = dict(zip(merged['object_id'], merged['signal']))
    sg.gdf['signal'] = sg.gdf['object_id'].map(signal_map)

    gdf_plot = sg.gdf.copy()
    if cluster_csv and os.path.exists(cluster_csv):
        cl = pd.read_csv(cluster_csv)
        if {'object_id', 'cluster_id'}.issubset(cl.columns):
            cl['object_id'] = pd.to_numeric(cl['object_id'], errors='coerce').astype('Int64')
            cl['cluster_id'] = cl['cluster_id'].astype(str)
            gdf_plot['leiden'] = gdf_plot['object_id'].map(dict(zip(cl['object_id'], cl['cluster_id']))).astype('category')
            before = len(gdf_plot)
            gdf_plot = gdf_plot[gdf_plot['leiden'].astype(str).isin(_KEEP_CLUSTERS)].copy()
            print(f"  Kept cluster(s) {_KEEP_CLUSTERS}: {len(gdf_plot)} cells (removed {before - len(gdf_plot)})")
    else:
        print("  [warn] cluster CSV not found; no cluster filtering applied.")

    min_x, max_x, min_y, max_y = crop
    roi_bbox = shapely_box(min_x, min_y, max_x, max_y)

    vmin = float(np.nanmin(gdf_plot['signal'])) if gdf_plot['signal'].notna().any() else 0.0
    vmax = float(np.nanmax(gdf_plot['signal'])) if gdf_plot['signal'].notna().any() else 1.0

    fig, ax = plt.subplots(figsize=(3, 3))

    # gray outlines for all cells in ROI
    sg.gdf[sg.gdf.geometry.intersects(roi_bbox)].boundary.plot(
        ax=ax, color='0.8', linewidth=0.2, alpha=0.9
    )
    # colored overlay for kept cells
    gdf_plot[gdf_plot.geometry.intersects(roi_bbox)].plot(
        column='signal', ax=ax, cmap='viridis',
        vmin=vmin, vmax=vmax, linewidth=0, edgecolor='none', legend=False
    )

    ax.set_xlim([min_x, max_x]); ax.set_ylim([min_y, max_y])
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.invert_yaxis()

    try:
        xmin, _ = ax.get_xlim()
        y_top, _ = ax.get_ylim()
        _add_scale_bar(ax, min_x=xmin, max_y=y_top,
                       length_um=_SCALE_BAR_UM, margin_px=_SCALE_BAR_MARGIN,
                       px_size_nm=_PIXEL_SIZE_NM)
    except Exception as e:
        print(f"  [warn] scale bar failed: {e}")

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm   = cm.ScalarMappable(norm=norm, cmap='viridis')
    cax  = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label('relative expression of villus top genes')
    for sp in ax.spines.values():
        sp.set_visible(False)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_top_bottom_villus_expression(output_dir=None, raw_data_root=None):
    """Generate monolayer_top_bottom_villus_expression.pdf for the autonomous zonation panel."""
    _out_dir  = output_dir    or _DEFAULT_OUT_DIR
    _raw_root = raw_data_root or _DEFAULT_RAW_ROOT

    nuclei_seg    = os.path.join(_raw_root, 'attributes', 'segmentations', 'segmentation_20230809_nuclei.tiff')
    transcripts   = os.path.join(_raw_root, 'all transcripts', 'transcrips_20240925_refid_dummy.csv')
    annotated_cbg = os.path.join(_raw_root, 'cell by gene', 'top_bottom_villus_exp.csv')
    cluster_csv   = os.path.join(_raw_root, 'umap_figures_8', 'cell_by_gene_cluster_annotations.csv')
    cbg_cache     = os.path.join(_raw_root, 'umap_figures_8', 'cell_by_gene_cached.csv')
    out_path      = os.path.join(_out_dir, 'monolayer_top_bottom_villus_expression.pdf')

    os.makedirs(_out_dir, exist_ok=True)

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({
            'figure.titlesize': 8, 'figure.titleweight': 'bold',
            'axes.titlesize': 8,   'axes.titleweight': 'bold',
            'axes.labelsize': 8,   'axes.labelweight': 'bold',
            'ytick.labelsize': 8,  'xtick.labelsize': 8,
            'legend.fontsize': 8,  'font.family': 'Arial',
            'figure.figsize': (4, 3), 'savefig.dpi': 300,
            'pdf.fonttype': 42,    'ps.fonttype': 42,
            'text.usetex': False,
        })
        sg, cell_gene_df = _build_cell_gene_table(nuclei_seg, transcripts, _CROP, cbg_cache)
        merged = _merge_annotations(cell_gene_df, annotated_cbg)
        _plot_and_save(sg, merged, out_path, cluster_csv, _CROP)
    finally:
        matplotlib.rcParams.update(_saved_rc)


if __name__ == '__main__':
    plot_top_bottom_villus_expression()
