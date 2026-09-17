"""
scale_invariance_profiles.py
----------------------------
Plots CPM-normalised mean expression per ring vs distance from the basal border
for the Visium scale-invariance analysis (Figure 2).

Can be called two ways:

  1. From plotALL.py (default):
       from paper.plotScripts.scale_invariance_profiles import plot_scale_invariance_figures
       plot_scale_invariance_figures()          # uses VISIUM_DATA_ROOT from utils/constant.py

  2. From the command line (standalone):
       python paper/plotScripts/scale_invariance_profiles.py               # full ROI, rep2
       python paper/plotScripts/scale_invariance_profiles.py rep4          # full ROI, rep4
       python paper/plotScripts/scale_invariance_profiles.py rep2 subroi   # sub-ROI, rep2

Output
------
    Figures → paper/graphs/scale_invariance_plots/
    Scale invariance scores (CSV) → same folder
"""

import os
import sys
import json
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
from scipy.sparse import diags
from scipy.ndimage import gaussian_filter1d
from scipy.stats import pearsonr

# ── resolve data root and output dir from utils/constant.py (with CLI fallback) ──
try:
    from utils.constant import (
        VISIUM_DATA_ROOT, SCALE_INVARIANCE_PLOTS_FOLDER_PATH,
        SCALE_INVARIANCE_PANEL_D_GENES, SCALE_INVARIANCE_GENES_OF_INTEREST,
    )
    _DATA_ROOT_DEFAULT = Path(VISIUM_DATA_ROOT)
    _OUT_DIR_DEFAULT   = Path(SCALE_INVARIANCE_PLOTS_FOLDER_PATH)
except ImportError:
    _DATA_ROOT_DEFAULT = Path(
        "/Users/yaelheyman/Library/CloudStorage/"
        "GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/"
        "mouse visium data/versi data set"
    )
    _OUT_DIR_DEFAULT = None   # set per-call below if still None
    SCALE_INVARIANCE_PANEL_D_GENES = ["Ada", "Apoa4", "Plac8", "Mki67"]
    SCALE_INVARIANCE_GENES_OF_INTEREST = [
        "Ada", "Apoa4", "Plac8", "Mki67", "Enpep", "Apoa1", "Aldob", "Sis",
    ]

# ============================================================
# CONFIG — shared by both CLI and programmatic calls
# ============================================================
TIMEPOINTS    = ["day0", "day3", "day5"]
N_ITER        = 20
KERNEL_SIZE   = 20

# Gene lists now live in utils/constant.py (SCALE_INVARIANCE_* constants).
GENES_OF_INTEREST = SCALE_INVARIANCE_GENES_OF_INTEREST

ERRORBAR_MODE = "sem"   # "sd" = ±1 SD  |  "sem" = ±1 SEM

TIMEPOINT_COLORS = ["#8C1A6A", "#E8941A", "#1B7FA3"]
TP_LABEL = {
    "day0": "Pre-irradiation",
    "day3": "3 Days post-irradiation",
    "day5": "5 Days post-irradiation",
}

EXCLUDE_LAST_RING = True
SMOOTH_SIGMA      = 1.0
# ============================================================


def load_um_per_hires_px(rep, data_root):
    """
    Load scalefactors JSON for the given rep and return µm per hires-image pixel.
    Formula: microns_per_pixel / tissue_hires_scalef
    """
    gsm = {"rep1": "GSM9134407", "rep2": "GSM9134408",
           "rep3": "GSM9134409", "rep4": "GSM9134410"}[rep]
    raw = data_root / "GSE303705_RAW"
    for fname in [f"{gsm}_{rep}_scalefactors_json.json",
                  f"{gsm}_{rep}_scalefactors_json.json.gz"]:
        p = raw / fname
        if p.exists():
            opener = gzip.open if fname.endswith(".gz") else open
            with opener(p, "rt") as f:
                sf = json.load(f)
            um_per_fullres = sf["microns_per_pixel"]
            hires_scalef   = sf["tissue_hires_scalef"]
            um_per_hires   = um_per_fullres / hires_scalef
            print(f"  [{rep}] microns_per_pixel={um_per_fullres:.4f}, "
                  f"tissue_hires_scalef={hires_scalef:.4f} "
                  f"→ {um_per_hires:.2f} µm/hires-px, "
                  f"ring width ≈ {um_per_hires * KERNEL_SIZE / 2:.1f} µm")
            return um_per_hires
    raise FileNotFoundError(f"scalefactors JSON not found for {rep}")


def load_ring_profiles(rep, day, genes_of_interest, n_iter, exclude_last, subset, data_root):
    """
    Returns, for each gene, per-ring arrays (ordered lumen→basal):
        means[gene]   : (n_rings,)  mean CPM per ring
        stds[gene]    : (n_rings,)  std CPM per ring
        ns[gene]      : (n_rings,)  number of bins per ring
        ring_ids      : (n_rings,)  ring labels (1-based, lumen=1)
    """
    out = (data_root / rep / day
           / "infl_out" / f"inflation_num_iterations_{n_iter}"
           / subset)

    X        = sp.load_npz(out / "counts_genes_by_bins.npz")   # genes × bins
    feat_df  = pd.read_csv(out / "features.tsv", sep="\t", header=None, skiprows=1)
    all_genes = feat_df[1].astype(str).tolist()
    ring_map  = pd.read_csv(out / "barcode_to_ring_aligned.csv")

    # CPM-normalise every bin (in-place CSR row scaling)
    X_bg = X.T.tocsr().astype(np.float64)                      # bins × genes
    totals = np.array(X_bg.sum(axis=1)).flatten()
    totals[totals == 0] = 1
    row_idx = np.repeat(np.arange(X_bg.shape[0]), np.diff(X_bg.indptr))
    X_bg.data *= (1e6 / totals[row_idx])

    # subset to requested genes
    g2i     = {g: i for i, g in enumerate(all_genes)}
    found   = [g for g in genes_of_interest if g in g2i]
    missing = [g for g in genes_of_interest if g not in g2i]
    if missing:
        print(f"  [{day}] not in dataset: {missing}")

    sel_idx = [g2i[g] for g in found]
    X_sel   = X_bg[:, sel_idx].toarray()                       # bins × n_found (dense)

    rings        = ring_map["ring"].to_numpy()
    unique_rings = sorted(np.unique(rings))
    if exclude_last:
        unique_rings = unique_rings[:-1]

    means, stds, ns = {}, {}, {}
    for gi, gene in enumerate(found):
        expr = X_sel[:, gi]
        m, s, n = [], [], []
        for r in unique_rings:
            vals = expr[rings == r]
            m.append(float(vals.mean()))
            s.append(float(vals.std(ddof=1)) if len(vals) > 1 else 0.0)
            n.append(int(len(vals)))
        means[gene] = np.array(m)
        stds[gene]  = np.array(s)
        ns[gene]    = np.array(n)

    return unique_rings, means, stds, ns, found


def plot_gene(gene, all_data, timepoints, palette, errorbar_mode, um_per_hires_px,
              ax=None, normalize_x=False):
    """
    Plot one gene onto ax.
    X-axis:
      normalize_x=False → physical distance from basal border in µm
      normalize_x=True  → fractional position (0=basal, 1=lumen)
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
    else:
        fig = ax.get_figure()

    err_label     = "SEM" if errorbar_mode == "sem" else "SD"
    ring_width_um = (KERNEL_SIZE / 2) * um_per_hires_px

    for i, day in enumerate(timepoints):
        d = all_data[day]
        if gene not in d["means"]:
            continue

        mean_raw = d["means"][gene][::-1]   # reverse: index 0 = basal
        std_raw  = d["stds"][gene][::-1]
        n_raw    = d["ns"][gene][::-1]
        n_rings  = len(mean_raw)

        err = std_raw / np.sqrt(n_raw) if errorbar_mode == "sem" else std_raw

        if normalize_x:
            x = np.linspace(0, 1, n_rings)
        else:
            ring_idx = np.arange(1, n_rings + 1)
            x = (ring_idx - 0.5) * ring_width_um

        color = palette[i]
        ax.errorbar(x, mean_raw,
                    yerr=err,
                    color=color,
                    linewidth=1.8,
                    marker="o",
                    markersize=5,
                    capsize=3,
                    capthick=1.2,
                    elinewidth=1.0,
                    label=f"{TP_LABEL[day]}  (n={n_rings} rings)",
                    zorder=3)

    if normalize_x:
        ax.axvline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ylims = ax.get_ylim()
        ax.text(0.02, ylims[1], "basal", fontsize=8, color="gray", va="top")
        ax.set_xlabel("Fractional position (0=basal, 1=lumen)", fontsize=8)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(-0.05, 1.05)
    else:
        ax.axvline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ylims = ax.get_ylim()
        ax.text(ring_width_um * 0.05, ylims[1], "basal", fontsize=8, color="gray", va="top")
        ax.set_xlabel("Distance from basal border (µm)", fontsize=8)
        all_n = [len(all_data[d]["means"].get(gene, [])) for d in timepoints]
        max_rings = max(all_n) if all_n else 1
        ax.set_xlim(-ring_width_um * 0.3, (max_rings + 0.2) * ring_width_um)

    ax.set_ylabel(f"Mean CPM  (±1 {err_label})", fontsize=8)
    ax.set_title(gene, fontsize=8, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if standalone:
        ax.set_box_aspect(1)

    return fig


def plot_gene_ribbon(gene, all_data, timepoints, palette, errorbar_mode, um_per_hires_px,
                     ax=None, normalize_x=False, normalize=False, smooth_sigma=None):
    """
    Ribbon-style plot (Moor et al. 2018 style): solid line + shaded ±SEM band.
    normalize=True → divide each trace by its own max CPM (y = fraction of max).
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
    else:
        fig = ax.get_figure()

    err_label     = "SEM" if errorbar_mode == "sem" else "SD"
    ring_width_um = (KERNEL_SIZE / 2) * um_per_hires_px

    for i, day in enumerate(timepoints):
        d = all_data[day]
        if gene not in d["means"]:
            continue

        mean_raw = d["means"][gene][::-1]
        std_raw  = d["stds"][gene][::-1]
        n_raw    = d["ns"][gene][::-1]
        n_rings  = len(mean_raw)

        err = std_raw / np.sqrt(n_raw) if errorbar_mode == "sem" else std_raw

        sigma = SMOOTH_SIGMA if smooth_sigma is None else smooth_sigma
        if sigma > 0 and n_rings >= 4:
            mean_raw = gaussian_filter1d(mean_raw, sigma=sigma)
            err      = gaussian_filter1d(err,      sigma=sigma)

        if normalize_x:
            x = np.linspace(0, 1, n_rings)
        else:
            ring_idx = np.arange(1, n_rings + 1, dtype=float)
            x = (ring_idx - 0.5) * ring_width_um

        if normalize:
            trace_max = float(mean_raw.max()) or 1.0
            mean_raw  = mean_raw / trace_max
            err       = err      / trace_max

        xp, yp, ylo, yhi = x, mean_raw, mean_raw - err, mean_raw + err

        color = palette[i]
        ax.fill_between(xp, ylo, yhi, color=color, alpha=0.20, linewidth=0)
        ax.plot(xp, yp, color=color, linewidth=2.2,
                label=f"{TP_LABEL[day]}  (n={n_rings} rings)", zorder=3)

    if normalize_x:
        ax.axvline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ylims = ax.get_ylim()
        ax.text(0.02, ylims[1], "basal", fontsize=8, color="gray", va="top")
        ax.set_xlabel("Fractional position (0=basal, 1=lumen)", fontsize=8)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(-0.05, 1.05)
    else:
        ax.axvline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ylims = ax.get_ylim()
        ax.text(ring_width_um * 0.05, ylims[1], "basal", fontsize=8, color="gray", va="top")
        ax.set_xlabel("Distance from basal border (µm)", fontsize=8)
        all_n = [len(all_data[d]["means"].get(gene, [])) for d in timepoints]
        max_rings = max(all_n) if all_n else 1
        ax.set_xlim(-ring_width_um * 0.3, (max_rings + 0.2) * ring_width_um)

    if normalize:
        ax.set_ylabel(f"Fraction of max CPM  (±1 {err_label})", fontsize=8)
        ylo_cur, yhi_cur = ax.get_ylim()
        ax.set_ylim(min(ylo_cur, -0.02), max(yhi_cur, 1.05) * 1.05)
    else:
        ax.set_ylabel(f"Mean CPM  (±1 {err_label})", fontsize=8)

    ax.set_title(gene, fontsize=8, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if standalone:
        ax.set_box_aspect(1)

    return fig


# ============================================================
# Scale invariance scoring
# ============================================================

def _curve_correlations(means_by_day, days, pairs, ring_width_um, norm_fn, n_grid):
    raw_c  = {}
    frac_c = {}
    for day in days:
        mean = means_by_day[day]
        n    = len(mean)
        x_um   = (np.arange(1, n + 1) - 0.5) * ring_width_um
        x_frac = np.linspace(0, 1, n)
        yn     = norm_fn(mean)
        raw_c[day]  = (x_um,   yn)
        frac_c[day] = (x_frac, yn)

    x_max = max(raw_c[d][0][-1] for d in days)
    x_grid_raw  = np.linspace(0, x_max,  n_grid)
    x_grid_frac = np.linspace(0, 1,      n_grid)
    raw_i  = {d: np.interp(x_grid_raw,  *raw_c[d])  for d in days}
    frac_i = {d: np.interp(x_grid_frac, *frac_c[d]) for d in days}

    out = {}
    for di, dj in pairs:
        if di not in days or dj not in days:
            continue
        r_raw,  _ = pearsonr(raw_i[di],  raw_i[dj])
        r_norm, _ = pearsonr(frac_i[di], frac_i[dj])
        out[f"{di}_vs_{dj}"] = (float(r_raw), float(r_norm))
    return out


def score_scale_invariance(all_data, genes, timepoints, um_per_hires_px,
                           apply_ynorm=True, n_boot=1000, seed=42):
    """
    For each gene compute mean pairwise Pearson r in physical µm space (r_raw)
    and fractional position space (r_norm), and Δr = r_norm − r_raw.
    Bootstrap 95% CI via parametric resampling from Normal(mean, SEM).
    """
    ring_width_um = (KERNEL_SIZE / 2) * um_per_hires_px
    rng   = np.random.default_rng(seed)
    pairs = [(timepoints[i], timepoints[j])
             for i in range(len(timepoints))
             for j in range(i + 1, len(timepoints))]

    def ynorm(y):
        m = float(y.max())
        return y / m if m > 0 else y

    norm_fn = ynorm if apply_ynorm else (lambda y: y)

    summary_rows = []
    pair_rows    = []

    for gene in genes:
        obs_means   = {}
        n_rings_day = {}
        sem_day     = {}

        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m  = d["means"][gene][::-1]
            s  = d["stds"][gene][::-1]
            n  = d["ns"][gene][::-1]
            obs_means[day]   = m
            n_rings_day[day] = len(m)
            sem_day[day]     = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        r_raw_obs, r_norm_obs = [], []
        for di, dj in pairs:
            if di not in days or dj not in days:
                continue
            n_eff  = min(n_rings_day[di], n_rings_day[dj])
            corrs  = _curve_correlations(obs_means, [di, dj],
                                         [(di, dj)], ring_width_um, norm_fn, n_eff)
            key    = f"{di}_vs_{dj}"
            r_r, r_n = corrs[key]
            r_raw_obs.append(r_r)
            r_norm_obs.append(r_n)
            pair_rows.append({
                "gene":    gene,
                "pair":    key,
                "n_rings": n_eff,
                "r_raw":   round(r_r,       3),
                "r_norm":  round(r_n,       3),
                "delta_r": round(r_n - r_r, 3),
            })

        obs_delta = float(np.mean(r_norm_obs)) - float(np.mean(r_raw_obs))

        boot_deltas = []
        for _ in range(n_boot):
            boot_means = {day: np.maximum(
                              rng.normal(obs_means[day], sem_day[day]), 0.0)
                          for day in days}
            br, bn = [], []
            for di, dj in pairs:
                if di not in days or dj not in days:
                    continue
                n_eff = min(n_rings_day[di], n_rings_day[dj])
                try:
                    c = _curve_correlations(boot_means, [di, dj],
                                            [(di, dj)], ring_width_um, norm_fn, n_eff)
                    r_r, r_n = c[f"{di}_vs_{dj}"]
                    br.append(r_r)
                    bn.append(r_n)
                except Exception:
                    pass
            if br:
                boot_deltas.append(float(np.mean(bn)) - float(np.mean(br)))

        ci_lo, ci_hi = (np.percentile(boot_deltas, [2.5, 97.5])
                        if boot_deltas else (np.nan, np.nan))

        summary_rows.append({
            "gene":        gene,
            "mean_r_raw":  round(float(np.mean(r_raw_obs)),  3),
            "mean_r_norm": round(float(np.mean(r_norm_obs)), 3),
            "delta_r":     round(obs_delta,                  3),
            "CI95_lo":     round(float(ci_lo),               3),
            "CI95_hi":     round(float(ci_hi),               3),
        })

    summary_df = (pd.DataFrame(summary_rows)
                  .sort_values("delta_r", ascending=False)
                  .reset_index(drop=True))
    pair_df = pd.DataFrame(pair_rows)
    return summary_df, pair_df


# ============================================================
# Main entry point (called by plotALL.py or CLI)
# ============================================================


def plot_scale_invariance_figures(rep="rep2", subroi=True,
                                   data_root=None, out_dir=None):
    """
    Produce Figure 2 (scale invariance, Visium).
    Only output: figure_2_scale_invariance.svg in out_dir.
    """
    _data_root = Path(data_root) if data_root else _DATA_ROOT_DEFAULT
    if out_dir:
        _out_dir = Path(out_dir)
    elif _OUT_DIR_DEFAULT is not None:
        _out_dir = _OUT_DIR_DEFAULT
    else:
        _out_dir = _data_root / rep / ("subroi_profiles" if subroi else "expression_profiles")

    subset = "ring_subset_bin_by_gene_subroi" if subroi else "ring_subset_bin_by_gene"
    suffix = "_subroi" if subroi else ""
    palette = TIMEPOINT_COLORS[:len(TIMEPOINTS)]

    # scope matplotlib style changes to this function only
    _saved_rc = dict(matplotlib.rcParams)
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["font.family"]  = "Arial"
    matplotlib.rcParams["font.size"]    = 8

    try:
        print(f"Loading data for {rep} (N_ITER={N_ITER}, subset={subset}): {TIMEPOINTS}")
        um_per_hires_px = load_um_per_hires_px(rep, _data_root)

        GENES_D = SCALE_INVARIANCE_PANEL_D_GENES
        all_data = {}
        for day in TIMEPOINTS:
            print(f"  {day}…")
            rings, means, stds, ns, found = load_ring_profiles(
                rep, day, GENES_D, N_ITER, EXCLUDE_LAST_RING, subset, _data_root)
            all_data[day] = {"rings": rings, "means": means, "stds": stds, "ns": ns}
        print("Done.\n")

        os.makedirs(_out_dir, exist_ok=True)

        # ── Panel D: raw CPM, 2 rows (µm x-axis / fractional x-axis) ─────────
        n_c       = len(GENES_D)
        _fig_w_in = 4.5 * n_c
        _fig_h_in = 4.5 * 2
        _comp_fs  = 8.0 * (_fig_w_in * 25.4) / 195.0

        fig_r, ax_r = plt.subplots(2, n_c, figsize=(_fig_w_in, _fig_h_in), squeeze=False)
        for gi, gene in enumerate(GENES_D):
            for ri, (norm_x, norm_y) in enumerate([(False, False), (True, False)]):
                plot_gene_ribbon(gene, all_data, TIMEPOINTS, palette, ERRORBAR_MODE,
                                 um_per_hires_px, ax=ax_r[ri][gi],
                                 normalize_x=norm_x, normalize=norm_y, smooth_sigma=0.5)
                ax = ax_r[ri][gi]
                if ax.get_legend():
                    ax.get_legend().remove()
                ax.set_title(ax.get_title(), fontsize=_comp_fs, fontweight="bold")
                ax.set_xlabel(ax.get_xlabel(), fontsize=_comp_fs)
                ax.set_ylabel(ax.get_ylabel() if gi == 0 else "", fontsize=_comp_fs)
                ax.tick_params(labelsize=_comp_fs)
                for txt in ax.texts:
                    txt.set_fontsize(_comp_fs)

        for ri, lbl in enumerate(["Mean CPM  (µm x-axis)", "Mean CPM  (fractional x-axis)"]):
            ax_r[ri][0].set_ylabel(lbl, fontsize=_comp_fs)

        fig_r.tight_layout()
        panel_d_path = _out_dir / f"{rep}_panel_d_raw_2row{suffix}_ribbon.svg"
        fig_r.savefig(panel_d_path, bbox_inches="tight", format="svg")
        plt.close(fig_r)
        print(f"Panel D saved: {panel_d_path.name}")

        # ── Panel C: same genes, y-normalised (fraction of max CPM) ──────────
        # Ported from monolayer_analysis/plot_rep_profiles.py, which produced this
        # panel before the July refactor. The paper's Figure 2 uses the raw-CPM
        # Panel D above, so this only feeds the y-normalised composite variant:
        # set EMIT_PANEL_C = True here and EMIT_YNORM_VARIANT = True in
        # scale_invariance_figure.py to build that one.
        EMIT_PANEL_C = False
        if EMIT_PANEL_C:
            fig_c, ax_c = plt.subplots(2, n_c, figsize=(_fig_w_in, _fig_h_in), squeeze=False)
            for gi, gene in enumerate(GENES_D):
                for ri, (norm_x, norm_y) in enumerate([(False, True), (True, True)]):
                    plot_gene_ribbon(gene, all_data, TIMEPOINTS, palette, ERRORBAR_MODE,
                                     um_per_hires_px, ax=ax_c[ri][gi],
                                     normalize_x=norm_x, normalize=norm_y, smooth_sigma=0.5)
                    ax = ax_c[ri][gi]
                    if ax.get_legend():
                        ax.get_legend().remove()
                    ax.set_title(ax.get_title(), fontsize=_comp_fs, fontweight="bold")
                    ax.set_xlabel(ax.get_xlabel(), fontsize=_comp_fs)
                    ax.set_ylabel(ax.get_ylabel() if gi == 0 else "", fontsize=_comp_fs)
                    ax.tick_params(labelsize=_comp_fs)
                    for txt in ax.texts:
                        txt.set_fontsize(_comp_fs)

            for ri, lbl in enumerate(["Fraction of max CPM  (µm x-axis)",
                                      "Fraction of max CPM  (fractional x-axis)"]):
                ax_c[ri][0].set_ylabel(lbl, fontsize=_comp_fs)

            fig_c.tight_layout()
            panel_c_path = _out_dir / f"{rep}_panel_c_3row{suffix}_ribbon.svg"
            fig_c.savefig(panel_c_path, bbox_inches="tight", format="svg")
            plt.close(fig_c)
            print(f"Panel C saved: {panel_c_path.name}")

        # ── Sub-ROI width stats (Panel C data) ───────────────────────────────
        if subroi:
            from paper.plotScripts.scale_invariance_subroi_width import process_timepoint
            _mask_dir = Path(__file__).parent.parent / "extractedData" / "visium_masks"
            print("\nComputing sub-ROI width stats …")
            for day in TIMEPOINTS:
                process_timepoint(rep, day, N_ITER, data_root=_data_root, mask_dir=_mask_dir)

        # ── Assemble Figure 2 ─────────────────────────────────────────────────
        if subroi:
            from paper.plotScripts.scale_invariance_figure import create_figure_2
            print("\nAssembling Figure 2 …")
            create_figure_2(rep=rep, data_root=_data_root, out_dir=_out_dir)

        print(f"\nFigure 2 saved to: {_out_dir / 'figure_2_scale_invariance.svg'}")

    finally:
        matplotlib.rcParams.update(_saved_rc)


# ============================================================
# CLI entry point
# ============================================================
if __name__ == "__main__":
    _rep    = sys.argv[1] if len(sys.argv) > 1 else "rep2"
    _subroi = len(sys.argv) > 2 and sys.argv[2].lower() == "subroi"
    plot_scale_invariance_figures(rep=_rep, subroi=_subroi)
