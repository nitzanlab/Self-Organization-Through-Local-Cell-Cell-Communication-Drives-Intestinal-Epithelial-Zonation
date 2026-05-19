"""
crop_monolayer.py
-----------------
Creates two raw monolayer figures for the autonomous zonation panel:
  1. monolayer_cell_type_markers.pdf — 8 cell-type-marker genes, wide crop
  2. monolayer_zonation.pdf          — 2 zonation genes (Alpi, Apoa4), zoomed crop

Call plot_monolayer_raw_figures() to generate both figures.
"""

import os
import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib
import matplotlib.pyplot as plt

try:
    from utils.constant import (
        AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH,
        MONOLAYER_RAW_DATA_ROOT,
        MONOLAYER_BACKGROUND_IMAGE,
    )
    _DEFAULT_OUT_DIR      = AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH
    _DEFAULT_RAW_ROOT     = MONOLAYER_RAW_DATA_ROOT
    _DEFAULT_BG_IMAGE     = MONOLAYER_BACKGROUND_IMAGE
except ImportError:
    _DEFAULT_RAW_ROOT = (
        '/Users/yaelheyman/Library/CloudStorage/'
        'GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/pasadena_run_no_gel'
    )
    _DEFAULT_BG_IMAGE = (
        '/Users/yaelheyman/RajLab Dropbox/Yael Heyman/shared_yael/'
        'sg/pasadena run no gel/roi_1/hyb_background_aligned.tiff'
    )
    _DEFAULT_OUT_DIR = os.path.join(
        os.getcwd(), 'paper', 'graphs', 'autonomous_zonation_figure_plots'
    )

_GENE_COLOR_DICT = {
    'Ada':   [0.58, 0,    0.83],
    'Lyz1':  [1,    0,    1   ],
    'Muc2':  [1,    0.45, 0   ],
    'Chga':  [0,    0,    1   ],
    'Dclk1': [1,    1,    1   ],
    'Mki67': [0,    1,    1   ],
    'Lgr5':  [1,    0,    0   ],
    'Alpi':  [0.20, 1.00, 0.30],
    'Apoa4': [0.85, 0.30, 0.95],
}

_GENE_SIZE_DICT = {
    'Alpi':  0.5,
    'Apoa4': 0.5,
    'Ada':   0.3,
    'Lyz1':  0.7,
    'Muc2':  0.7,
    'Chga':  0.7,
    'Dclk1': 0.7,
    'Mki67': 0.5,
    'Lgr5':  0.9,
}


def _show_neighborhood_with_adjusted_scalebar(
    min_x, max_x, min_y, max_y,
    image_path, transcripts_file, genes,
    output_file_path,
    figure_size=(3, 3),
    dpi=600,
    brightness_factor=5.0,
    legend_marker_size=5,
    legend_font_size=6,
    legend_box_scale=0.5,
    gene_size_dict=None,
    scale_bar_length_um=100,
    pixel_size_nm=107.11,
    scale_bar_color='white',
    scale_bar_height=70,
):
    transcripts = pd.read_csv(transcripts_file, index_col=0)
    if isinstance(genes, str):
        genes = [genes]

    spot_coords_list = []
    for gene in genes:
        if gene in transcripts.index:
            gene_coords = transcripts.loc[[gene], ['x', 'y']].copy()
            gene_coords['gene'] = gene
            spot_coords_list.append(gene_coords)
        else:
            print(f"  Gene '{gene}' not found in transcripts.")
    if not spot_coords_list:
        print("  No valid genes found — skipping.")
        return
    spot_coords = pd.concat(spot_coords_list).reset_index(drop=True)

    image = tiff.imread(image_path)   # shape (C, Y, X)
    pad_y = (max_y - min_y) / 4
    pad_x = (max_x - min_x) / 4
    ry0 = max(0,              int(min_y - pad_y))
    ry1 = min(image.shape[1], int(max_y + pad_y))
    rx0 = max(0,              int(min_x - pad_x))
    rx1 = min(image.shape[2], int(max_x + pad_x))

    in_range = (
        (spot_coords['y'] > ry0) & (spot_coords['y'] < ry1) &
        (spot_coords['x'] > rx0) & (spot_coords['x'] < rx1)
    )
    spot_coords_in_range = spot_coords[in_range]

    gray = np.array(image[3, ry0:ry1, rx0:rx1])
    norm = np.clip(
        (gray - gray.min()) / (gray.max() - gray.min()) * 255 * brightness_factor,
        0, 255,
    ).astype(np.uint8)

    blue_img = np.zeros((*gray.shape, 3), dtype=np.uint8)
    blue_img[..., 2] = norm

    fig, ax = plt.subplots()
    ax.imshow(blue_img)
    ax.axis('off')

    default_spot_size = 50
    for gene in genes:
        gdf = spot_coords_in_range[spot_coords_in_range['gene'] == gene]
        if gdf.empty:
            continue
        spot_size = (gene_size_dict or {}).get(gene, default_spot_size)
        ax.scatter(
            gdf['x'] - rx0, gdf['y'] - ry0,
            color=_GENE_COLOR_DICT.get(gene, [1, 1, 1]),
            s=spot_size,
            edgecolor='none',
            linewidths=0,
            alpha=0.7,
            label=gene,
        )

    fig.set_size_inches(figure_size[0], figure_size[1])

    scale_bar_px = int(scale_bar_length_um * 1000 / pixel_size_nm)
    ax.add_patch(plt.Rectangle(
        (200, gray.shape[0] - scale_bar_height - 300),
        scale_bar_px, scale_bar_height,
        color=scale_bar_color, lw=0,
    ))

    legend = ax.legend(
        markerscale=legend_marker_size,
        fontsize=legend_font_size,
        handletextpad=0.4,
        labelspacing=0.2,
        ncol=1,
        loc='lower right',
        borderpad=0.8,
    )
    legend.get_frame().set_linewidth(0.5)
    legend.get_frame().set_facecolor('lightgrey')
    legend.get_frame().set_alpha(legend_box_scale)

    os.makedirs(os.path.dirname(os.path.abspath(output_file_path)), exist_ok=True)
    plt.savefig(output_file_path, dpi=dpi, bbox_inches='tight', pad_inches=0)
    plt.close(fig)


def plot_monolayer_raw_figures(output_dir=None, raw_data_root=None, background_image_path=None):
    """Generate the two raw monolayer crop figures for the autonomous zonation panel."""
    _out_dir     = output_dir            or _DEFAULT_OUT_DIR
    _raw_root    = raw_data_root         or _DEFAULT_RAW_ROOT
    _bg_image    = background_image_path or _DEFAULT_BG_IMAGE
    _transcripts = os.path.join(_raw_root, 'all transcripts', 'transcrips_20240925.csv')

    os.makedirs(_out_dir, exist_ok=True)

    _saved_rc = dict(matplotlib.rcParams)
    try:
        print("  Generating intestinal_monolayers_express_canonical_cell_type_markers.pdf …")
        _show_neighborhood_with_adjusted_scalebar(
            min_x=100, max_x=8000, min_y=47118, max_y=50469,
            image_path=_bg_image,
            transcripts_file=_transcripts,
            genes=['Alpi', 'Apoa4', 'Muc2', 'Chga', 'Lgr5', 'Dclk1', 'Lyz1', 'Mki67'],
            output_file_path=os.path.join(_out_dir, 'intestinal_monolayers_express_canonical_cell_type_markers.pdf'),
            figure_size=(3, 3),
            dpi=600,
            brightness_factor=5.0,
            legend_marker_size=5,
            legend_font_size=6,
            legend_box_scale=0.5,
            gene_size_dict=_GENE_SIZE_DICT,
            scale_bar_length_um=100,
            pixel_size_nm=107.11,
            scale_bar_color='white',
            scale_bar_height=200,
        )

        print("  Generating intestinal_monolayers_are_zonated.pdf …")
        _show_neighborhood_with_adjusted_scalebar(
            min_x=9219, max_x=17780, min_y=14269, max_y=20986,
            image_path=_bg_image,
            transcripts_file=_transcripts,
            genes=['Alpi', 'Apoa4'],
            output_file_path=os.path.join(_out_dir, 'intestinal_monolayers_are_zonated.pdf'),
            figure_size=(3, 3),
            dpi=600,
            brightness_factor=5.0,
            legend_marker_size=5,
            legend_font_size=6,
            legend_box_scale=0.5,
            gene_size_dict=_GENE_SIZE_DICT,
            scale_bar_length_um=100,
            pixel_size_nm=107.11,
            scale_bar_color='white',
            scale_bar_height=200,
        )
        print(f"  Monolayer figures saved to: {_out_dir}")
    finally:
        matplotlib.rcParams.update(_saved_rc)


if __name__ == "__main__":
    plot_monolayer_raw_figures()
