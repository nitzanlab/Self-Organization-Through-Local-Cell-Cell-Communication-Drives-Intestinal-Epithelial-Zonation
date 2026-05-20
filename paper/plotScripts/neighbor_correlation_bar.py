"""
neighbor_correlation_bar.py
----------------------------
Bar plot of neighbor-correlation difference (72hr − 12hr) per gene,
colored by Moran's I in the unperturbed monolayer (Figure 3 panel F).

Output: neighbor_correlation_morans_i_bar.pdf → ZONATION_PLASTICITY_PLOTS_FOLDER_PATH

Call plot_neighbor_correlation_bar() to generate the figure.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

try:
    from utils.constant import ZONATION_PLASTICITY_PLOTS_FOLDER_PATH, SPRINKLING_NOV23_BASE_PATH
    _DEFAULT_OUT_DIR = ZONATION_PLASTICITY_PLOTS_FOLDER_PATH
    _DEFAULT_BASE    = SPRINKLING_NOV23_BASE_PATH
except ImportError:
    _DEFAULT_BASE    = '/Users/yaelheyman/Library/CloudStorage/GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/sprinkling_nov_23'
    _DEFAULT_OUT_DIR = os.path.join(os.getcwd(), 'paper', 'graphs', 'zonation_plasticity_plots')


def plot_neighbor_correlation_bar(output_dir=None, base_path=None):
    """Bar plot of 72hr−12hr neighbor-correlation difference per gene, colored by Moran's I."""
    _out  = output_dir or _DEFAULT_OUT_DIR
    _base = base_path  or _DEFAULT_BASE
    os.makedirs(_out, exist_ok=True)

    corr_path    = os.path.join(_base, '72hr', 'roi1', 'output', 'cell_by_gene', 'correlation_by_gene_100_.csv')
    morans_path  = os.path.join(_base, '72hr', 'roi1', 'output', 'morans_i_wt_monolayer.csv')

    corr_df      = pd.read_csv(corr_path)
    morans_df    = pd.read_csv(morans_path)
    morans_i     = morans_df.set_index('Unnamed: 0')['morans_i']

    time_point_comp = pd.DataFrame({
        "12": corr_df.loc[corr_df["id"] == 12]["corr_coeff"].to_numpy(),
        "72": corr_df.loc[corr_df["id"] == 72]["corr_coeff"].to_numpy(),
    }, index=corr_df.loc[corr_df["id"] == 12]["gene"])

    diff        = (time_point_comp['72'] - time_point_comp['12']).sort_values()
    x_labels    = diff.index.tolist()
    common      = [g for g in x_labels if g in morans_i.index]
    sorted_mi   = morans_i.loc[common].values

    norm   = plt.Normalize(vmin=sorted_mi.min(), vmax=sorted_mi.max())
    colors = plt.cm.viridis(norm(sorted_mi))

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({
            'figure.titlesize': 8, 'figure.titleweight': 'bold',
            'axes.titlesize': 8,   'axes.titleweight': 'bold',
            'axes.labelsize': 8,   'axes.labelweight': 'bold',
            'ytick.labelsize': 8,  'xtick.labelsize': 8,
            'legend.fontsize': 8,  'font.family': 'Arial',
            'savefig.dpi': 300,    'pdf.fonttype': 42, 'ps.fonttype': 42,
            'text.usetex': False,
        })

        fig, ax = plt.subplots(figsize=(4.1, 2))
        ax.bar(range(len(common)), diff.loc[common].values, color=colors)

        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label="Moran's I")

        offset = 4
        idxs   = np.arange(offset, len(common), 5)
        ax.set_xticks(idxs)
        ax.set_xticklabels([common[i] for i in idxs], rotation=90, fontsize=7)
        ax.set_ylabel('Neighbor correlation difference\n(72hr - 12hr)', fontsize=8)

        out_path = os.path.join(_out, 'neighbor_correlation_morans_i_bar.pdf')
        plt.savefig(out_path, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        print(f'  Saved: {out_path}')
    finally:
        matplotlib.rcParams.update(_saved_rc)


if __name__ == '__main__':
    plot_neighbor_correlation_bar()
