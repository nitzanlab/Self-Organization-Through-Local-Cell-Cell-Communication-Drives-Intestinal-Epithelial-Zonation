"""
pharmacological_perturbations.py
---------------------------------
Plotting functions for the pharmacological perturbations figure (Figure 4).

Panel E: Per-condition grouped bar chart — mean ± SD cell counts (Sis, Msln, Apoa4)
         across 3 experiments, one sub-panel per condition.
Panel F: Ring-colored polygon image — cell segmentations colored by distance from
         monolayer edge (ENR condition, frame 3).
Panel G: All-conditions %Double-Positive line plot — Apoa4+ & Sis+ / Apoa4+ [%]
         vs. distance from edge, all conditions overlaid (same color per condition,
         different markers per experiment).
Panel H: In vivo EphA2 heatmap (produced by autonomous_zonation_plots.py, saved here).

Entry points
------------
  plot_perturbation_ring_image()          → Panel F
  plot_perturbation_cell_count_bars()     → Panel E
  plot_perturbation_double_positive()     → Panel G
  plot_all_pharmacological_perturbation_plots() → run all three
"""

import os
import re
import glob
import json
import colorsys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, Patch
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

try:
    from utils.constant import (
        PHARMACOLOGICAL_PERTURBATIONS_PLOTS_FOLDER_PATH,
        PHARMACOLOGICAL_PERTURBATION_EXPERIMENT_BASES,
        PHARMACOLOGICAL_PERTURBATION_PANEL_F_POLYGONS_PATH,
        PHARMACOLOGICAL_PERTURBATION_PANEL_F_CSV_PATH,
    )
    _DEFAULT_OUT_DIR       = PHARMACOLOGICAL_PERTURBATIONS_PLOTS_FOLDER_PATH
    _DEFAULT_EXP_BASES     = PHARMACOLOGICAL_PERTURBATION_EXPERIMENT_BASES
    _DEFAULT_PANEL_F_POLY  = PHARMACOLOGICAL_PERTURBATION_PANEL_F_POLYGONS_PATH
    _DEFAULT_PANEL_F_CSV   = PHARMACOLOGICAL_PERTURBATION_PANEL_F_CSV_PATH
except ImportError:
    _SPRINKLING = "/Users/yaelheyman/RajLab Dropbox/Yael Heyman/shared_yael/Zonation"
    _DEFAULT_OUT_DIR = os.path.join(os.getcwd(), "paper", "graphs", "pharmacological_perturbations")
    _DEFAULT_EXP_BASES = [
        os.path.join(_SPRINKLING, "20250529_monolayer_conditions_re", "different_conditions"),
        os.path.join(_SPRINKLING, "20250718_monolayer_conditions", "different_conditions"),
        os.path.join(_SPRINKLING, "20250719_monolayer_conditions", "different_conditions"),
    ]
    _DEFAULT_PANEL_F_POLY = os.path.join(
        _SPRINKLING, "20250529_monolayer_conditions_re", "different_conditions",
        "ENR", "polygons_by_frame", "frame_2_polygons.json"
    )
    _DEFAULT_PANEL_F_CSV = os.path.join(
        _SPRINKLING, "20250529_monolayer_conditions_re", "different_conditions",
        "ENR", "processedData", "erosion_analysis",
        "frame_3", "num_iterations_11", "cell_by_gene_with_ring.csv"
    )

# ── shared style ──────────────────────────────────────────────────────────────
_GENES      = ["Sis", "Msln", "Apoa4"]
_THRESHOLDS = {"Sis": 11000, "Msln": 8000, "Apoa4": 15000}

_EXP_LABELS = {
    "20250529_monolayer_conditions_re": "Experiment 1",
    "20250718_monolayer_conditions":    "Experiment 2",
    "20250719_monolayer_conditions":    "Experiment 3",
}
_EXP_COLORS = {
    "20250529_monolayer_conditions_re": "#9467bd",
    "20250718_monolayer_conditions":    "#ff7f0e",
    "20250719_monolayer_conditions":    "#8c564b",
}
_MARKERS = {
    "20250529_monolayer_conditions_re": "o",
    "20250718_monolayer_conditions":    "^",
    "20250719_monolayer_conditions":    "s",
}
_KEEP_CONDITIONS = {"ENR", "IWP LOW", "LDN LOW", "ALW LOW", "R0 HIGH"}
_COND_COLORS = {
    "ENR":      "#9467bd",
    "IWP LOW":  "#8c564b",
    "LDN LOW":  "#ff7f0e",
    "ALW LOW":  "#17becf",
    "R0 HIGH":  "#2ca02c",
}
_COND_DISPLAY = {
    "ENR":      "Standard media",
    "LDN LOW":  "BMP inhibitor",
    "IWP LOW":  "WNT inhibitor",
    "ALW LOW":  "Ephrin inhibitor",
    "R0 HIGH":  "Notch inhibitor",
}
_COND_ORDER = ["ENR", "LDN LOW", "IWP LOW", "ALW LOW", "R0 HIGH"]
_AVG_RING_WIDTH_UM = 8.6

_ITER_DIR_GLOB    = os.path.join("processedData", "erosion_analysis", "overview_plots_iterations_*")
_PERCENT_CSV_NAME = "percent_apoa4_and_sis_out_of_apoa4.csv"
_COUNTS_CSV_NAME  = "threshold_histogram_counts_per_frame.csv"


# ── helpers ───────────────────────────────────────────────────────────────────

def _shades_of(base_color, n, lightest=0.65, darkest=0.0):
    """Return n shades of base_color going from light to dark.

    lightest/darkest control how much whiteness is added (0=none, 1=fully white).
    """
    r, g, b = mcolors.to_rgb(base_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    shades = []
    for i in range(n):
        t = i / max(n - 1, 1)         # 0 → lightest, 1 → darkest
        blend = lightest + t * (darkest - lightest)
        new_l = min(1.0, l + blend * (1.0 - l))
        r2, g2, b2 = colorsys.hls_to_rgb(h, new_l, s)
        shades.append((r2, g2, b2))
    return shades


def _exp_key(base_path):
    tail = os.path.basename(base_path.rstrip(os.sep))
    if tail == "different_conditions":
        return os.path.basename(os.path.dirname(base_path.rstrip(os.sep)))
    return tail


def _canonical_condition(name):
    s = name.strip().upper().replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\bB2\b", "BMP2", s)
    s = re.sub(r"\bB4\b", "BMP4", s)
    s = re.sub(r"\bH\b",  "HIGH", s)
    s = re.sub(r"\bL\b",  "LOW",  s)
    s = re.sub(r"\b(BMP2|BMP4|LDN|ALW|WNT|IWP|R0)(HIGH|LOW)\b", r"\1 \2", s)
    return re.sub(r"\s+", " ", s).strip()


def _find_iteration_dirs(base, raw_cond):
    cond_dir = os.path.join(base, raw_cond)
    out = {}
    for iterdir in glob.glob(os.path.join(cond_dir, _ITER_DIR_GLOB)):
        m = re.search(r"overview_plots_iterations_(\d+)$", os.path.basename(iterdir))
        if m:
            out[int(m.group(1))] = iterdir
    return out


def _build_canonical_map(experiment_bases):
    per_exp = {}
    for base in experiment_bases:
        key = _exp_key(base)
        mapping = {}
        if not os.path.isdir(base):
            per_exp[key] = mapping
            continue
        for raw in os.listdir(base):
            if os.path.isdir(os.path.join(base, raw)):
                mapping.setdefault(_canonical_condition(raw), []).append(raw)
        per_exp[key] = mapping
    return per_exp


# ── Panel F: ring-colored polygon image ──────────────────────────────────────

def _load_polygons(json_path):
    with open(json_path) as f:
        data = json.load(f)
    polys = []
    for item in data:
        arr = np.asarray(item, dtype=float)
        if arr.ndim == 2 and arr.shape[1] == 2 and len(arr) >= 3:
            polys.append(arr)
    if not polys:
        raise ValueError(f"No valid polygons in {json_path}")
    return polys


def _poly_centroid(P):
    x, y   = P[:, 0], P[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    cross  = x * y2 - x2 * y
    A      = cross.sum() / 2.0
    if abs(A) < 1e-9:
        return np.array([x.mean(), y.mean()], dtype=float)
    cx = ((x + x2) * cross).sum() / (6 * A)
    cy = ((y + y2) * cross).sum() / (6 * A)
    return np.array([cx, cy], dtype=float)


def _compute_bounds(polys):
    all_pts = np.vstack(polys)
    xmin, ymin = all_pts.min(axis=0)
    xmax, ymax = all_pts.max(axis=0)
    return xmin, xmax, ymin, ymax


def _unique_nn_match(src_points, dst_points):
    tree   = cKDTree(src_points)
    dists, idx = tree.query(dst_points, k=1)
    order  = np.argsort(dists)
    taken  = set()
    dst2src = np.full(len(dst_points), -1, dtype=int)
    for j in order:
        s = idx[j]
        if s not in taken:
            taken.add(s)
            dst2src[j] = s
    return dst2src, dists


def _add_scale_bar(ax, scale_bar_length_um=100, pixel_size_um=543.2 / 100,
                   where="tr", pad_frac=0.04, height_px_frac=0.2):
    scale_bar_px = scale_bar_length_um * pixel_size_um
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    pad_x  = pad_frac * abs(x1 - x0)
    pad_y  = pad_frac * abs(y1 - y0)
    bar_h  = scale_bar_px * height_px_frac
    if where == "tr":
        bl = max(x0, x1) - pad_x - scale_bar_px
        bb = max(y0, y1) - pad_y - bar_h
    else:
        bl = min(x0, x1) + pad_x
        bb = min(y0, y1) + pad_y
    ax.add_patch(plt.Rectangle((bl, bb), scale_bar_px, bar_h,
                                color="white", lw=0, zorder=10))


def plot_perturbation_ring_image(output_dir=None,
                                  polygons_path=None,
                                  csv_path=None,
                                  no_ring_alpha=0.5):
    """Panel F: cell polygons colored by erosion ring (distance from edge), black background."""
    _out   = output_dir   or _DEFAULT_OUT_DIR
    _poly  = polygons_path or _DEFAULT_PANEL_F_POLY
    _csv   = csv_path      or _DEFAULT_PANEL_F_CSV
    os.makedirs(_out, exist_ok=True)

    polys = _load_polygons(_poly)

    df    = pd.read_csv(_csv)
    xcol  = next(c for c in df.columns if "centroid" in c.lower() and c.strip().lower().replace(" ", "").endswith("/x"))
    ycol  = next(c for c in df.columns if "centroid" in c.lower() and c.strip().lower().replace(" ", "").endswith("/y"))
    cells_xy = df[[xcol, ycol]].to_numpy(float)
    rings    = df["Ring"].to_numpy()

    poly_cents  = np.vstack([_poly_centroid(P) for P in polys])
    dst2src, d  = _unique_nn_match(poly_cents, cells_xy)

    order = np.argsort(d)
    poly_to_ring = {}
    for j in order:
        pidx = dst2src[j]
        if pidx >= 0 and pidx not in poly_to_ring:
            poly_to_ring[pidx] = int(rings[j])

    pos_rings   = sorted(r for r in set(poly_to_ring.values()) if r > 0)
    tab         = plt.get_cmap("tab10")
    ring_colors = {r: tab(i % tab.N) for i, r in enumerate(pos_rings)}
    ring_colors[0] = (0.7, 0.7, 0.7, float(no_ring_alpha))

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({
            "pdf.fonttype": 42, "ps.fonttype": 42, "text.usetex": False,
        })
        fig, ax = plt.subplots(figsize=(5, 5))
        fig.patch.set_facecolor("black")
        ax.set_facecolor("black")

        for pidx, r in poly_to_ring.items():
            if r == 0:
                ax.add_patch(MplPolygon(polys[pidx], closed=True,
                                        facecolor=ring_colors[0], edgecolor="none",
                                        linewidth=0, zorder=1))
        for pidx, r in poly_to_ring.items():
            if r > 0:
                ax.add_patch(MplPolygon(polys[pidx], closed=True,
                                        facecolor=ring_colors[r], edgecolor="none",
                                        linewidth=0, zorder=2))

        xmin, xmax, ymin, ymax = _compute_bounds(polys)
        margin = 5
        ax.set_xlim(xmin - margin, xmax + margin)
        ax.set_ylim(ymax + margin, ymin - margin)
        ax.set_aspect("equal")
        ax.axis("off")
        _add_scale_bar(ax)

        handles = [Patch(color=ring_colors[0], label="Ring 0")] + \
                  [Patch(color=ring_colors[r], label=f"Ring {r}") for r in pos_rings]
        leg = ax.legend(handles=handles, loc="lower left", framealpha=0.6, fontsize=7)
        leg.get_frame().set_facecolor("black")
        leg.get_frame().set_edgecolor("white")
        for txt in leg.get_texts():
            txt.set_color("white")

        out_path = os.path.join(_out, "ENR_ring_colored_polygons.pdf")
        fig.savefig(out_path, dpi=300, bbox_inches="tight",
                    facecolor="black", edgecolor="none")
        plt.close(fig)
        print(f"  Saved: {out_path}")
    finally:
        matplotlib.rcParams.update(_saved_rc)


# ── Panel E: per-condition grouped bar chart ─────────────────────────────────

def plot_perturbation_cell_count_bars(output_dir=None, experiment_bases=None):
    """Panel E: all conditions in one combined row figure, shaded by experiment."""
    _out   = output_dir        or _DEFAULT_OUT_DIR
    _bases = experiment_bases  or _DEFAULT_EXP_BASES
    os.makedirs(_out, exist_ok=True)

    per_exp_map   = _build_canonical_map(_bases)
    available     = {c for m in per_exp_map.values() for c in m} & _KEEP_CONDITIONS
    conds_ordered = [c for c in _COND_ORDER if c in available]
    exp_order     = [_exp_key(b) for b in _bases]

    # Collect all data: iter_num -> cond -> exp_key -> DataFrame
    all_data = {}
    for cond in conds_ordered:
        exp_iters = {}
        for base in _bases:
            key      = _exp_key(base)
            raw_list = per_exp_map.get(key, {}).get(cond, [])
            merged   = {}
            for raw in raw_list:
                for k, v in _find_iteration_dirs(base, raw).items():
                    merged.setdefault(k, v)
            if merged:
                exp_iters[key] = merged

        for key, imap in exp_iters.items():
            for iter_num, iter_dir in imap.items():
                counts_path = os.path.join(iter_dir, _COUNTS_CSV_NAME)
                if not os.path.isfile(counts_path):
                    continue
                try:
                    dfc = pd.read_csv(counts_path)
                    dfc = dfc[[g for g in _GENES if g in dfc.columns]]
                    all_data.setdefault(iter_num, {}).setdefault(cond, {})[key] = dfc
                except Exception:
                    pass

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({
            "pdf.fonttype": 42, "ps.fonttype": 42, "text.usetex": False,
            "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
            "font.family": "Arial",
        })

        for iter_num, cond_data in sorted(all_data.items()):
            conds_in_plot = [c for c in conds_ordered if c in cond_data]
            if not conds_in_plot:
                continue

            n    = len(conds_in_plot)
            fig, axes = plt.subplots(1, n, figsize=(1.6 * n, 2.3), sharey=True)
            if n == 1:
                axes = [axes]

            fig.suptitle("Mean cell counts under different perturbations",
                         fontsize=8, y=1.03)

            width = 0.22
            x     = np.arange(len(_GENES))

            for ax_idx, (ax, cond) in enumerate(zip(axes, conds_in_plot)):
                shades         = _shades_of(_COND_COLORS.get(cond, "gray"), len(_bases))
                per_exp_counts = cond_data[cond]

                for i, (key, shade) in enumerate(zip(exp_order, shades)):
                    dfc = per_exp_counts.get(key)
                    means_, stds_ = [], []
                    for g in _GENES:
                        if dfc is None or g not in dfc.columns:
                            means_.append(np.nan); stds_.append(np.nan)
                        else:
                            means_.append(dfc[g].mean()); stds_.append(dfc[g].std())

                    ax.bar(
                        x + (i - (len(_bases) - 1) / 2) * width,
                        means_, width=width,
                        yerr=stds_, capsize=3,
                        label=_EXP_LABELS.get(key, key),
                        color=shade,
                        error_kw=dict(elinewidth=0.8, capthick=0.8),
                    )

                ax.set_xticks(x)
                ax.set_xticklabels(_GENES, rotation=30, ha="right", fontsize=7)
                ax.set_title(_COND_DISPLAY.get(cond, cond), fontsize=7)
                ax.tick_params(labelsize=7)
                ax.yaxis.set_major_formatter(lambda val, _: f"{val:.0e}")
                ax.set_ylim(0, 1e4)
                ax.spines[["top", "right"]].set_visible(False)

                if ax_idx == 0:
                    ax.set_ylabel("Number of cells", fontsize=7)

                legend_handles = [
                    Patch(color=s, label=_EXP_LABELS.get(k, k))
                    for k, s in zip(exp_order, shades)
                ]
                ax.legend(handles=legend_handles, fontsize=5, frameon=False,
                          loc="upper right", handlelength=1, handletextpad=0.3,
                          borderpad=0.2)

            fig.tight_layout()
            out_path = os.path.join(_out, f"cell_counts_all_conditions_iter{iter_num}.pdf")
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {out_path}")
    finally:
        matplotlib.rcParams.update(_saved_rc)


# ── Panel G: all-conditions %Double-Positive line plot ───────────────────────

def plot_perturbation_double_positive(output_dir=None, experiment_bases=None):
    """Panel G: %Double-Positive (Apoa4+ & Sis+)/Apoa4+ vs distance from edge,
    all conditions overlaid (color=condition, marker=experiment)."""
    _out   = output_dir        or _DEFAULT_OUT_DIR
    _bases = experiment_bases  or _DEFAULT_EXP_BASES
    os.makedirs(_out, exist_ok=True)

    per_exp_map   = _build_canonical_map(_bases)
    available     = {c for m in per_exp_map.values() for c in m} & _KEEP_CONDITIONS
    conds_sorted  = [c for c in _COND_ORDER if c in available]

    combined_by_iter = {}
    for cond in conds_sorted:
        exp_iters = {}
        for base in _bases:
            key      = _exp_key(base)
            raw_list = per_exp_map.get(key, {}).get(cond, [])
            merged   = {}
            for raw in raw_list:
                for k, v in _find_iteration_dirs(base, raw).items():
                    merged.setdefault(k, v)
            if merged:
                exp_iters[key] = merged

        for key, imap in exp_iters.items():
            for iter_num, iter_dir in imap.items():
                pct_path = os.path.join(iter_dir, _PERCENT_CSV_NAME)
                if not os.path.isfile(pct_path):
                    continue
                try:
                    dfp = pd.read_csv(pct_path)
                except Exception:
                    continue
                if not {"Ring", "%Double_Positive"}.issubset(dfp.columns):
                    continue
                dfp = dfp[["Ring", "%Double_Positive"]].dropna().sort_values("Ring")
                if not dfp.empty:
                    combined_by_iter \
                        .setdefault(iter_num, {}) \
                        .setdefault(cond, {})[key] = dfp

    _saved_rc = dict(matplotlib.rcParams)
    try:
        matplotlib.rcParams.update({
            "pdf.fonttype": 42, "ps.fonttype": 42, "text.usetex": False,
            "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
            "font.family": "Arial",
        })

        for iter_num, cond_dict in sorted(combined_by_iter.items()):
            if not cond_dict:
                continue

            fig, ax = plt.subplots(figsize=(6, 4))
            ax.tick_params(labelsize=8)
            ax.set_title("Percent Double-Positive", fontsize=9)

            # Pre-compute shades per condition (same scheme as Panel E)
            cond_shades = {
                c: _shades_of(_COND_COLORS.get(c, "gray"), len(_bases))
                for c in conds_sorted
            }
            _exp_order = [_exp_key(b) for b in _bases]

            for cond in conds_sorted:
                exp_map = cond_dict.get(cond, {})
                shades  = cond_shades[cond]
                for i, base in enumerate(_bases):
                    key = _exp_key(base)
                    dfp = exp_map.get(key)
                    if dfp is None or dfp.empty:
                        continue
                    ax.plot(
                        dfp["Ring"].values * _AVG_RING_WIDTH_UM,
                        dfp["%Double_Positive"].values,
                        marker=_MARKERS.get(key, "o"),
                        linestyle="-",
                        color=shades[i],
                        linewidth=1.2,
                        markersize=3.5,
                    )

            ax.set_xlabel("Distance from monolayer edge (µm)", fontsize=8)
            ax.set_ylabel("Apoa4+ & Sis+ / Apoa4+ [%]", fontsize=8)
            ax.set_ylim(0, 100)
            ax.grid(True, alpha=0.3)
            ax.spines[["top", "right"]].set_visible(False)

            # Conditions legend: show darkest shade (Exp 3) as representative
            cond_handles = [
                Line2D([0], [0], color=cond_shades[c][-1], lw=2,
                       label=_COND_DISPLAY.get(c, c))
                for c in conds_sorted if c in cond_dict
            ]
            exp_handles = [
                Line2D([0], [0], marker=_MARKERS.get(k, "o"), color="black",
                       linestyle="None", label=_EXP_LABELS.get(k, k), markersize=5)
                for k in _exp_order
            ]
            leg1 = ax.legend(handles=cond_handles, title="Conditions", fontsize=6,
                             title_fontsize=7, frameon=False,
                             loc="upper left", bbox_to_anchor=(1.02, 1.0))
            ax.add_artist(leg1)
            ax.legend(handles=exp_handles, title="Experiments", fontsize=6,
                      title_fontsize=7, frameon=False,
                      loc="lower left", bbox_to_anchor=(1.02, 0.0))

            fig.tight_layout()
            out_path = os.path.join(_out, f"double_positive_all_conditions_iter{iter_num}.pdf")
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {out_path}")
    finally:
        matplotlib.rcParams.update(_saved_rc)


# ── main entry point ──────────────────────────────────────────────────────────

def plot_all_pharmacological_perturbation_plots(output_dir=None, experiment_bases=None):
    """Generate all pharmacological perturbation panels (E, F, G)."""
    _out   = output_dir       or _DEFAULT_OUT_DIR
    _bases = experiment_bases or _DEFAULT_EXP_BASES

    print("\nGenerating Panel F: ring-colored polygon image …")
    plot_perturbation_ring_image(output_dir=_out)

    print("\nGenerating Panel E: per-condition cell count bar charts …")
    plot_perturbation_cell_count_bars(output_dir=_out, experiment_bases=_bases)

    print("\nGenerating Panel G: %Double-Positive line plot …")
    plot_perturbation_double_positive(output_dir=_out, experiment_bases=_bases)


if __name__ == "__main__":
    plot_all_pharmacological_perturbation_plots()
