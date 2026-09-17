#!/usr/bin/env python3
"""
validate_scale_invariance_scoring.py
--------------------------------------
Tests two scale-invariance scoring methods against synthetic NEGATIVE CONTROLS
(day0 profiles truncated to day3/day5 widths — same shape, just cut short).

SCORING METHODS
---------------
1. Δr      (original)  = mean_r_norm − mean_r_raw
   Pairwise Pearson r in fractional space minus Pearson r in µm space.
   Positive = profiles align better in fractional space = scale-invariant.
   LIMITATION: r is insensitive to shape distortions; a truncated monotone
   profile correlates strongly in fractional space even though it isn't SI.

2. Δ(RMSE) (new)       = mean_RMSE_µm − mean_RMSE_frac
   Pairwise RMSE in µm space minus RMSE in fractional space (both y-normalised).
   Positive = profiles are CLOSER in fractional space = scale-invariant.
   MORE ROBUST: RMSE captures pointwise shape differences that r misses.
   A truncated profile stretched to [0,1] has a different shape than the full
   profile → higher RMSE_frac relative to RMSE_µm.

SYNTHETIC NEGATIVE CONTROL
---------------------------
  synthetic_day0      = day0 full profile  (N_d0 rings)
  synthetic_day3_trunc = day0 truncated to N_d3 rings  ← basal side kept
  synthetic_day5_trunc = day0 truncated to N_d5 rings  ← basal side kept

A good method should give score ≤ 0 for ALL synthetic genes.
Any gene with CI entirely > 0 is a FALSE POSITIVE.

Usage:  python validate_scale_invariance_scoring.py [rep]   (default: rep2)
"""

import sys
import json
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"]  = "Arial"
matplotlib.rcParams["font.size"]    = 8
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

# ── CONFIG ────────────────────────────────────────────────────────────────────
REP        = sys.argv[1] if len(sys.argv) > 1 else "rep2"
_tp_arg    = sys.argv[2] if len(sys.argv) > 2 else "day0,day3,day5"
TIMEPOINTS = _tp_arg.split(",")
TP_SUFFIX  = "_".join(TIMEPOINTS)   # used in output filenames
# SYNTH_TPS built later once we know which timepoints are in use
N_ITER      = 20
KERNEL_SIZE = 20
N_BOOT      = 1000

GENES_10 = [
    "Ada", "Apoa4", "Sis", "Pigr", "Aldob",
    "Mki67", "Enpep", "Apoa1", "Fabp1", "Lgr5",
]

COLORS_REAL  = {"day0": "#8C1A6A", "day3": "#E8941A", "day5": "#1B7FA3"}
COLORS_SYNTH = {"day0": "#8C1A6A", "day3_trunc": "#E8941A", "day5_trunc": "#1B7FA3"}
TP_LABEL_REAL  = {
    "day0": "Pre-irradiation",
    "day3": "3d post-irradiation",
    "day5": "5d post-irradiation",
}
TP_LABEL_SYNTH = {
    "day0":       "day0 (full)",
    "day3_trunc": "day0 → day3 width",
    "day5_trunc": "day0 → day5 width",
}

GSM_BY_REP = {
    "rep1": "GSM9134407", "rep2": "GSM9134408",
    "rep3": "GSM9134409", "rep4": "GSM9134410",
}
try:
    from utils.constant import VISIUM_DATA_ROOT, SCALE_INVARIANCE_PLOTS_FOLDER_PATH
    DATA_ROOT = Path(VISIUM_DATA_ROOT)
except ImportError:
    DATA_ROOT = Path(
        "/Users/yaelheyman/Library/CloudStorage/"
        "GoogleDrive-yaelhei@gmail.com/My Drive/SPRINKLING/SG/"
        "mouse visium data/versi data set"
    )
    SCALE_INVARIANCE_PLOTS_FOLDER_PATH = str(DATA_ROOT / REP / "scale_invariance_validation")
RAW_DIR = DATA_ROOT / "GSE303705_RAW"
OUT_DIR = Path(SCALE_INVARIANCE_PLOTS_FOLDER_PATH) / "scale_invariance_validation"
SUBSET  = "ring_subset_bin_by_gene_subroi"

# ── DATA LOADING ──────────────────────────────────────────────────────────────

def load_um_per_hires_px():
    gsm = GSM_BY_REP[REP]
    for fname in [f"{gsm}_{REP}_scalefactors_json.json",
                  f"{gsm}_{REP}_scalefactors_json.json.gz"]:
        p = RAW_DIR / fname
        if p.exists():
            opener = gzip.open if fname.endswith(".gz") else open
            with opener(p, "rt") as fh:
                sf = json.load(fh)
            return sf["microns_per_pixel"] / sf["tissue_hires_scalef"]
    raise FileNotFoundError(f"scalefactors JSON not found for {REP}")


def load_ring_profiles(day, genes_of_interest):
    out = (DATA_ROOT / REP / day
           / "infl_out" / f"inflation_num_iterations_{N_ITER}" / SUBSET)

    X         = sp.load_npz(out / "counts_genes_by_bins.npz")
    feat_df   = pd.read_csv(out / "features.tsv", sep="\t", header=None, skiprows=1)
    all_genes = feat_df[1].astype(str).tolist()
    ring_map  = pd.read_csv(out / "barcode_to_ring_aligned.csv")

    X_bg   = X.T.tocsr().astype(np.float64)
    totals = np.array(X_bg.sum(axis=1)).flatten()
    totals[totals == 0] = 1
    row_idx = np.repeat(np.arange(X_bg.shape[0]), np.diff(X_bg.indptr))
    X_bg.data *= (1e6 / totals[row_idx])

    g2i     = {g: i for i, g in enumerate(all_genes)}
    found   = [g for g in genes_of_interest if g in g2i]
    missing = [g for g in genes_of_interest if g not in g2i]
    if missing:
        print(f"  [{day}] not in dataset: {missing}")

    X_sel  = X_bg[:, [g2i[g] for g in found]].toarray()
    rings  = ring_map["ring"].to_numpy()
    uring  = sorted(np.unique(rings))[:-1]   # exclude catch-all last ring

    means, stds, ns = {}, {}, {}
    for gi, gene in enumerate(found):
        expr = X_sel[:, gi]
        m, s, n = [], [], []
        for r in uring:
            vals = expr[rings == r]
            m.append(float(vals.mean()))
            s.append(float(vals.std(ddof=1)) if len(vals) > 1 else 0.0)
            n.append(int(len(vals)))
        means[gene] = np.array(m)
        stds[gene]  = np.array(s)
        ns[gene]    = np.array(n)

    return uring, means, stds, ns, found


def load_mean_width_rings(day, ring_width_um):
    ws = pd.read_csv(DATA_ROOT / REP / day / "infl_out"
                     / f"inflation_num_iterations_{N_ITER}" / "subroi_width_stats.csv")
    mean_um = float(ws["mean_px"].iloc[0]) * UM_PER_PX
    return max(1, int(round(mean_um / ring_width_um))), mean_um


def make_truncated_data(real_day0, n_rings_target):
    """
    Truncate day0 to n_rings_target rings, keeping the BASAL side.
    Storage is lumenal-first → keep last n_rings_target elements.
    """
    means_t, stds_t, ns_t = {}, {}, {}
    for gene in real_day0["means"]:
        m, s, n = (real_day0[k][gene] for k in ("means", "stds", "ns"))
        k = min(n_rings_target, len(m))
        means_t[gene] = m[len(m) - k:]
        stds_t[gene]  = s[len(s) - k:]
        ns_t[gene]    = n[len(n) - k:]
    return {"means": means_t, "stds": stds_t, "ns": ns_t}


# ── SCORING ───────────────────────────────────────────────────────────────────

N_GRID = 100   # interpolation grid points for Euclidean distance scoring

def _build_interp_grids(means_by_day, days, ring_width_um, n_grid,
                        zero_pad=False):
    """
    Interpolate raw CPM profiles for each day onto common µm and
    fractional grids of length n_grid. Returns (raw_i, frac_i) dicts.
    zero_pad: if True, shorter profiles are zero-padded beyond their extent
              in the raw µm grid (right=0.0); otherwise held at last value.
    """
    raw_c, frac_c = {}, {}
    for day in days:
        mean   = means_by_day[day]   # basal-first
        n      = len(mean)
        x_um   = (np.arange(1, n + 1) - 0.5) * ring_width_um
        x_frac = np.linspace(0, 1, n)
        raw_c[day]  = (x_um,   mean)
        frac_c[day] = (x_frac, mean)

    x_max       = max(raw_c[d][0][-1] for d in days)
    x_grid_raw  = np.linspace(0, x_max, n_grid)
    x_grid_frac = np.linspace(0, 1,     n_grid)
    right_val   = 0.0 if zero_pad else None
    raw_i  = {d: np.interp(x_grid_raw,  *raw_c[d], right=right_val)
              for d in days}
    frac_i = {d: np.interp(x_grid_frac, *frac_c[d]) for d in days}
    return raw_i, frac_i


def score_euclidean_dist(all_data, genes, timepoints, ring_width_um,
                         n_grid=N_GRID, n_boot=N_BOOT, seed=42):
    """
    For each gene:
      1. Y-normalise each timepoint profile (divide by its own max).
      2. Interpolate onto a common fractional grid of n_grid points.
      3. At each grid point x: compute mean pairwise |y_i(x) - y_j(x)| across all pairs.
      4. Mean over x → one score per gene.

    Lower score = profiles cluster tightly in fractional space = scale-invariant.
    Higher score = profiles diverge in fractional space = NOT scale-invariant.

    Bootstrap: parametric, sample ring means from Normal(observed, SEM).
    95% CI on the score.
    """
    rng   = np.random.default_rng(seed)
    pairs = [(timepoints[i], timepoints[j])
             for i in range(len(timepoints))
             for j in range(i + 1, len(timepoints))]
    x_grid = np.linspace(0, 1, n_grid)

    rows = []
    for gene in genes:
        obs_means, n_rings_day, sem_day = {}, {}, {}
        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]   # basal-first
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day]   = m
            n_rings_day[day] = len(m)
            sem_day[day]     = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        def _dist(means_dict):
            # interpolate raw CPM onto fractional grid
            interp = {}
            for day in days:
                x_day = np.linspace(0, 1, len(means_dict[day]))
                interp[day] = np.interp(x_grid, x_day, means_dict[day])
            # mean pairwise |diff| per grid point, then mean over grid
            pair_dists = []
            for di, dj in pairs:
                if di not in days or dj not in days:
                    continue
                pair_dists.append(np.abs(interp[di] - interp[dj]))
            if not pair_dists:
                return np.nan
            # shape: (n_pairs, n_grid) → mean over pairs at each x → mean over x
            return float(np.mean(np.stack(pair_dists)))

        obs_score = _dist(obs_means)

        # bootstrap
        boot_scores = []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            try:
                boot_scores.append(_dist(bm))
            except Exception:
                pass

        ci_lo, ci_hi = (np.percentile(boot_scores, [2.5, 97.5])
                        if boot_scores else (np.nan, np.nan))

        rows.append({
            "gene":     gene,
            "euc_dist": round(obs_score, 5),
            "CI95_lo":  round(float(ci_lo), 5),
            "CI95_hi":  round(float(ci_hi), 5),
        })

    return (pd.DataFrame(rows)
            .sort_values("euc_dist")   # ascending: lower = more scale-invariant
            .reset_index(drop=True))


def score_euclidean_zeropad(all_data, genes, timepoints, ring_width_um,
                            n_grid=N_GRID, n_boot=N_BOOT, seed=42):
    """
    Euclidean distance in ABSOLUTE ring space, with zero-padding.

    All profiles share a common grid spanning the longest timepoint's ring range.
    Profiles shorter than the longest are zero-padded at the lumenal end
    (np.interp right=0), NOT held at their last value.

    For each gene:
      1. Y-normalise each profile (divide by its own max).
      2. Place each profile at absolute positions [0 … n_rings] on the common grid.
         Lumenal positions beyond the profile's actual extent = 0.
      3. Mean pairwise |y_i(x) - y_j(x)| across pairs, then mean over grid points.

    Interpretation:
      LOW  score → profiles agree at the same absolute positions (non-scale-invariant)
      HIGH score → profiles differ in absolute space (consistent with scale invariance
                   or biological changes across timepoints)

    For the synthetic NEGATIVE CONTROL, day3_trunc/day5_trunc are EXACTLY day0
    within their extent, then zero → they should have LOW scores.
    Truly scale-invariant real genes compress their gradient → distance within
    the overlapping region → should have HIGHER scores than the synthetic controls.
    """
    rng   = np.random.default_rng(seed)
    pairs = [(timepoints[i], timepoints[j])
             for i in range(len(timepoints))
             for j in range(i + 1, len(timepoints))]

    rows = []
    for gene in genes:
        obs_means, n_rings_day, sem_day = {}, {}, {}
        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]   # basal-first
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day]   = m
            n_rings_day[day] = len(m)
            sem_day[day]     = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        n_max = max(n_rings_day[d] for d in days)
        # common absolute grid: 0 to n_max in ring units
        x_grid = np.linspace(0, n_max, n_grid)

        def _interp_zeropad(means_dict):
            interp = {}
            for day in days:
                raw  = means_dict[day]
                n    = len(raw)
                # ring centres at 0.5, 1.5, ..., n-0.5  (in ring units)
                x_day = np.arange(0.5, n, 1.0)
                # right=0: zero-pad beyond the last ring centre
                interp[day] = np.interp(x_grid, x_day, raw, right=0.0)
            return interp

        def _mean_dist(interp):
            dists = []
            for di, dj in pairs:
                if di not in days or dj not in days:
                    continue
                dists.append(np.abs(interp[di] - interp[dj]))
            if not dists:
                return np.nan
            return float(np.mean(np.stack(dists)))   # mean over pairs then over grid

        obs_interp = _interp_zeropad(obs_means)
        obs_score  = _mean_dist(obs_interp)

        boot_scores = []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            try:
                boot_scores.append(_mean_dist(_interp_zeropad(bm)))
            except Exception:
                pass

        ci_lo, ci_hi = (np.percentile(boot_scores, [2.5, 97.5])
                        if boot_scores else (np.nan, np.nan))

        rows.append({
            "gene":     gene,
            "euc_zpad": round(obs_score, 5),
            "CI95_lo":  round(float(ci_lo), 5),
            "CI95_hi":  round(float(ci_hi), 5),
        })

    # higher = more different in absolute space (more scale-invariant)
    return (pd.DataFrame(rows)
            .sort_values("euc_zpad", ascending=False)
            .reset_index(drop=True))


def score_delta_I(all_data, genes, timepoints, ring_width_um,
                  n_grid=N_GRID, n_boot=N_BOOT, seed=42):
    """
    Information-theoretic ΔI score adapted from Nikolic et al. 2024.

    Each profile is mapped to its own fractional [0,1] axis (raw CPM, no y-norm).
    At each fractional position x_s:

      σ²_between(x_s) = sample variance of mean CPM across timepoints
      σ²_meas(x_s)    = mean of SEM² across timepoints (measurement noise floor)

      ΔI = ½ ⟨ log₂( 1 + σ²_between / σ²_meas ) ⟩_{x_s}

    ΔI ≈ 0  →  profiles match in fractional space  →  scale-invariant.
    ΔI >> 0 →  profiles differ in fractional space  →  NOT scale-invariant.

    Works for monotone genes because a truncated profile stretched to [0,1]
    has different absolute CPM at each fractional position; without
    y-normalisation this drives σ²_between > 0.

    Bootstrap: parametric, sample ring means from Normal(observed, SEM).
    """
    rng    = np.random.default_rng(seed)
    x_grid = np.linspace(0, 1, n_grid)

    rows = []
    for gene in genes:
        obs_means, sem_day = {}, {}
        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day] = m
            sem_day[day]   = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        def _interp_frac(means_dict):
            return {day: np.interp(x_grid,
                                   np.linspace(0, 1, len(means_dict[day])),
                                   means_dict[day])
                    for day in days}

        def _interp_sem(sem_dict):
            return {day: np.interp(x_grid,
                                   np.linspace(0, 1, len(sem_dict[day])),
                                   sem_dict[day])
                    for day in days}

        def _dI(means_i, sem_i):
            M = np.stack([means_i[d] for d in days])    # (n_tp, n_grid)
            S = np.stack([sem_i[d]   for d in days])    # (n_tp, n_grid)
            sigma2_between = np.var(M,    axis=0, ddof=1)   # between-timepoint variance
            sigma2_meas    = np.mean(S**2, axis=0)           # mean measurement noise
            mask = sigma2_meas > 0
            if mask.sum() == 0:
                return np.nan
            ratio = 1.0 + sigma2_between[mask] / sigma2_meas[mask]
            return float(0.5 * np.mean(np.log2(ratio)))

        obs_mi    = _interp_frac(obs_means)
        obs_si    = _interp_sem(sem_day)
        obs_score = _dI(obs_mi, obs_si)

        boot_scores = []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            try:
                boot_scores.append(_dI(_interp_frac(bm), obs_si))
            except Exception:
                pass

        ci_lo, ci_hi = (np.percentile(boot_scores, [2.5, 97.5])
                        if boot_scores else (np.nan, np.nan))

        rows.append({
            "gene":    gene,
            "delta_I": round(obs_score, 5),
            "CI95_lo": round(float(ci_lo), 5),
            "CI95_hi": round(float(ci_hi), 5),
        })

    return (pd.DataFrame(rows)
            .sort_values("delta_I")   # lower = more scale-invariant
            .reset_index(drop=True))


def score_delta_I_alt(all_data, genes, timepoints, ring_width_um,
                      n_grid=N_GRID, n_boot=N_BOOT, seed=42):
    """
    ΔI with mean-expression normalization (fixes noise-floor bias in score_delta_I).

    At each fractional position x_s:

      σ²_between(x_s) = sample variance of mean CPM across timepoints
      μ(x_s)          = mean of mean CPM across timepoints

      ΔI_alt = ½ ⟨ log₂( 1 + σ²_between / μ² ) ⟩_{x_s}

    Normalising by μ² (instead of SEM²) makes the score scale-free and
    avoids conflating low-expression / noisy genes with scale-invariant ones.
    Points where μ ≈ 0 are skipped (no expression → not informative).

    Lower ΔI_alt ≈ 0 → profiles match in fractional space → scale-invariant.
    Higher ΔI_alt   → profiles differ in fractional space → NOT scale-invariant.
    """
    rng    = np.random.default_rng(seed)
    x_grid = np.linspace(0, 1, n_grid)

    rows = []
    for gene in genes:
        obs_means, sem_day = {}, {}
        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day] = m
            sem_day[day]   = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        def _interp(means_dict):
            return {day: np.interp(x_grid,
                                   np.linspace(0, 1, len(means_dict[day])),
                                   means_dict[day])
                    for day in days}

        def _dI_alt(means_i):
            M  = np.stack([means_i[d] for d in days])   # (n_tp, n_grid)
            mu = np.mean(M, axis=0)                       # mean across timepoints
            sigma2_between = np.var(M, axis=0, ddof=1)
            mask = mu > 0
            if mask.sum() == 0:
                return np.nan
            ratio = 1.0 + sigma2_between[mask] / mu[mask]**2
            return float(0.5 * np.mean(np.log2(ratio)))

        obs_score = _dI_alt(_interp(obs_means))

        boot_scores = []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            try:
                boot_scores.append(_dI_alt(_interp(bm)))
            except Exception:
                pass

        ci_lo, ci_hi = (np.percentile(boot_scores, [2.5, 97.5])
                        if boot_scores else (np.nan, np.nan))

        rows.append({
            "gene":      gene,
            "dI_alt":    round(obs_score, 5),
            "CI95_lo":   round(float(ci_lo), 5),
            "CI95_hi":   round(float(ci_hi), 5),
        })

    return (pd.DataFrame(rows)
            .sort_values("dI_alt")
            .reset_index(drop=True))


def score_normdist_fractional(all_data, genes, timepoints, ring_width_um,
                              n_grid=N_GRID, n_boot=N_BOOT, seed=42):
    """
    Pointwise normalized squared difference in FRACTIONAL x space:

      score = mean over pairs { mean over x∈[0,1] of (f(x)−g(x))² / ((f(x)+g(x))/2)² }

    Each profile is independently mapped to its own [0,1] fractional axis.
    No y-normalization; raw CPM. No zero-padding (all profiles fill [0,1]).
    Grid points where both profiles are zero are skipped.

    Lower score = profiles more similar in fractional space = scale-invariant.
    """
    rng   = np.random.default_rng(seed)
    pairs = [(timepoints[i], timepoints[j])
             for i in range(len(timepoints))
             for j in range(i + 1, len(timepoints))]
    x_grid = np.linspace(0, 1, n_grid)

    rows = []
    for gene in genes:
        obs_means, n_rings_day, sem_day = {}, {}, {}
        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day]   = m
            n_rings_day[day] = len(m)
            sem_day[day]     = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        def _interp(means_dict):
            out = {}
            for day in days:
                raw   = means_dict[day]
                x_day = np.linspace(0, 1, len(raw))
                out[day] = np.interp(x_grid, x_day, raw)
            return out

        def _norm_dist(interp):
            dists = []
            for di, dj in pairs:
                if di not in days or dj not in days:
                    continue
                f, g  = interp[di], interp[dj]
                denom = ((f + g) / 2.0) ** 2
                numer = (f - g) ** 2
                mask  = denom > 0
                if mask.sum() == 0:
                    continue
                dists.append(float(np.mean(numer[mask] / denom[mask])))
            return float(np.mean(dists)) if dists else np.nan

        obs_score = _norm_dist(_interp(obs_means))

        boot_scores = []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            try:
                boot_scores.append(_norm_dist(_interp(bm)))
            except Exception:
                pass

        ci_lo, ci_hi = (np.percentile(boot_scores, [2.5, 97.5])
                        if boot_scores else (np.nan, np.nan))

        rows.append({
            "gene":        gene,
            "normdist_frac": round(obs_score, 5),
            "CI95_lo":       round(float(ci_lo), 5),
            "CI95_hi":       round(float(ci_hi), 5),
        })

    return (pd.DataFrame(rows)
            .sort_values("normdist_frac")   # lower = more scale-invariant
            .reset_index(drop=True))


def score_euclidean_normalized(all_data, genes, timepoints, ring_width_um,
                               n_grid=N_GRID, n_boot=N_BOOT, seed=42):
    """
    Pointwise normalized squared difference (zero-padded absolute space):

      score = mean over pairs { mean over x of (f(x)−g(x))² / ((f(x)+g(x))/2)² }

    The denominator is the local mean squared at each grid point, making the
    metric independent of absolute CPM level. Grid points where both profiles
    are zero (zero-padded overhang region) are skipped.

    Higher score = profiles more different at same absolute positions
                 = consistent with scale invariance.
    Lower  score = profiles similar at same absolute positions (e.g. truncation).
    """
    rng   = np.random.default_rng(seed)
    pairs = [(timepoints[i], timepoints[j])
             for i in range(len(timepoints))
             for j in range(i + 1, len(timepoints))]

    rows = []
    for gene in genes:
        obs_means, n_rings_day, sem_day = {}, {}, {}
        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day]   = m
            n_rings_day[day] = len(m)
            sem_day[day]     = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        n_max  = max(n_rings_day[d] for d in days)
        x_grid = np.linspace(0, n_max, n_grid)

        def _interp(means_dict):
            out = {}
            for day in days:
                raw   = means_dict[day]
                x_day = np.arange(0.5, len(raw), 1.0)
                out[day] = np.interp(x_grid, x_day, raw, right=0.0)
            return out

        def _norm_dist(interp):
            dists = []
            for di, dj in pairs:
                if di not in days or dj not in days:
                    continue
                f, g   = interp[di], interp[dj]
                denom  = ((f + g) / 2.0) ** 2
                numer  = (f - g) ** 2
                mask   = denom > 0          # skip points where both are zero
                if mask.sum() == 0:
                    continue
                dists.append(float(np.mean(numer[mask] / denom[mask])))
            return float(np.mean(dists)) if dists else np.nan

        obs_interp = _interp(obs_means)
        obs_score  = _norm_dist(obs_interp)

        boot_scores = []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            try:
                boot_scores.append(_norm_dist(_interp(bm)))
            except Exception:
                pass

        ci_lo, ci_hi = (np.percentile(boot_scores, [2.5, 97.5])
                        if boot_scores else (np.nan, np.nan))

        rows.append({
            "gene":    gene,
            "norm_dist": round(obs_score, 5),
            "CI95_lo":   round(float(ci_lo), 5),
            "CI95_hi":   round(float(ci_hi), 5),
        })

    return (pd.DataFrame(rows)
            .sort_values("norm_dist", ascending=False)
            .reset_index(drop=True))


def _pair_scores(raw_i, frac_i, di, dj):
    """Return (r_raw, r_frac, rmse_raw, rmse_frac) for one pair."""
    r_raw,  _ = pearsonr(raw_i[di],  raw_i[dj])
    r_frac, _ = pearsonr(frac_i[di], frac_i[dj])
    rmse_raw  = float(np.sqrt(np.mean((raw_i[di]  - raw_i[dj])**2)))
    rmse_frac = float(np.sqrt(np.mean((frac_i[di] - frac_i[dj])**2)))
    return float(r_raw), float(r_frac), rmse_raw, rmse_frac


def score_both_methods(all_data, genes, timepoints, ring_width_um,
                       n_boot=N_BOOT, seed=42):
    """
    For each gene compute Δr and Δ(RMSE) with 95% bootstrap CI.

    Δr      = mean_r_frac  − mean_r_raw   (higher = more scale-invariant)
    Δ(RMSE) = mean_RMSE_µm − mean_RMSE_frac (higher = more scale-invariant)

    Both use raw CPM profiles (no y-normalisation).
    Grid size = min(n_rings across pair) for honest degrees of freedom.
    Bootstrap: parametric, sample ring means from Normal(observed, SEM).
    """
    rng   = np.random.default_rng(seed)
    pairs = [(timepoints[i], timepoints[j])
             for i in range(len(timepoints))
             for j in range(i + 1, len(timepoints))]

    summary_rows, pair_rows = [], []

    for gene in genes:
        obs_means, n_rings_day, sem_day = {}, {}, {}

        for day in timepoints:
            d = all_data[day]
            if gene not in d["means"]:
                continue
            m = d["means"][gene][::-1]   # basal-first
            s = d["stds"][gene][::-1]
            n = d["ns"][gene][::-1]
            obs_means[day]   = m
            n_rings_day[day] = len(m)
            sem_day[day]     = s / np.sqrt(np.maximum(n, 1))

        days = [d for d in timepoints if d in obs_means]
        if len(days) < 2:
            continue

        # ── observed scores ────────────────────────────────────────────────
        obs_r_raw, obs_r_frac = [], []
        obs_rmse_raw, obs_rmse_frac = [], []

        for di, dj in pairs:
            if di not in days or dj not in days:
                continue
            n_eff = max(n_rings_day[di], n_rings_day[dj])
            raw_i, frac_i = _build_interp_grids(obs_means, [di, dj],
                                                 ring_width_um, n_eff,
                                                 zero_pad=True)
            r_r, r_f, rmse_r, rmse_f = _pair_scores(raw_i, frac_i, di, dj)
            obs_r_raw.append(r_r);    obs_r_frac.append(r_f)
            obs_rmse_raw.append(rmse_r); obs_rmse_frac.append(rmse_f)
            pair_rows.append({
                "gene": gene, "pair": f"{di}_vs_{dj}", "n_rings": n_eff,
                "r_raw":    round(r_r,    3), "r_frac":    round(r_f,    3),
                "delta_r":  round(r_f - r_r, 3),
                "rmse_raw": round(rmse_r, 4), "rmse_frac": round(rmse_f, 4),
                "delta_rmse": round(rmse_r - rmse_f, 4),
            })

        obs_delta_r    = np.mean(obs_r_frac)   - np.mean(obs_r_raw)
        obs_delta_rmse = np.mean(obs_rmse_raw) - np.mean(obs_rmse_frac)

        # ── bootstrap CIs ─────────────────────────────────────────────────
        boot_dr, boot_drmse = [], []
        for _ in range(n_boot):
            bm = {day: np.maximum(
                rng.normal(obs_means[day], sem_day[day]), 0.0) for day in days}
            br_r, br_f, br_rr, br_rf = [], [], [], []
            for di, dj in pairs:
                if di not in days or dj not in days:
                    continue
                n_eff = max(n_rings_day[di], n_rings_day[dj])
                try:
                    ri, fi = _build_interp_grids(bm, [di, dj], ring_width_um, n_eff,
                                                 zero_pad=True)
                    r_r, r_f, rmse_r, rmse_f = _pair_scores(ri, fi, di, dj)
                    br_r.append(r_r); br_f.append(r_f)
                    br_rr.append(rmse_r); br_rf.append(rmse_f)
                except Exception:
                    pass
            if br_r:
                boot_dr.append(float(np.mean(br_f))  - float(np.mean(br_r)))
                boot_drmse.append(float(np.mean(br_rr)) - float(np.mean(br_rf)))

        def ci(vals):
            return (np.percentile(vals, [2.5, 97.5]) if vals
                    else (np.nan, np.nan))

        dr_lo,    dr_hi    = ci(boot_dr)
        drmse_lo, drmse_hi = ci(boot_drmse)

        summary_rows.append({
            "gene":          gene,
            # Δr
            "mean_r_raw":    round(float(np.mean(obs_r_raw)),   3),
            "mean_r_frac":   round(float(np.mean(obs_r_frac)),  3),
            "delta_r":       round(float(obs_delta_r),          3),
            "dr_CI95_lo":    round(float(dr_lo),                3),
            "dr_CI95_hi":    round(float(dr_hi),                3),
            # Δ(RMSE)
            "mean_rmse_raw": round(float(np.mean(obs_rmse_raw)),  4),
            "mean_rmse_frac":round(float(np.mean(obs_rmse_frac)), 4),
            "delta_rmse":    round(float(obs_delta_rmse),         4),
            "drmse_CI95_lo": round(float(drmse_lo),              4),
            "drmse_CI95_hi": round(float(drmse_hi),              4),
        })

    df = (pd.DataFrame(summary_rows)
          .sort_values("delta_rmse", ascending=False)
          .reset_index(drop=True))
    return df, pd.DataFrame(pair_rows)


# ── PLOTTING HELPERS ──────────────────────────────────────────────────────────

def ribbon_ax(ax, data, tps, colors, labels, ring_width_um,
              normalize_x=False, title=""):
    for tp in tps:
        m = data[tp]["means"].get(tp if tp in data else tp)
        # unpack properly
        m_arr = data[tp]["means"].get(list(data[tp]["means"].keys())[0])
        break

    for tp in tps:
        if tp not in data or not data[tp]["means"]:
            continue
        gene_means = data[tp]["means"]
        # this helper receives pre-selected single-gene arrays
        m = data[tp]["mean"]
        s = data[tp]["sem"]

        m_bf = m[::-1]
        s_bf = s[::-1]
        n_rings = len(m_bf)

        if normalize_x:
            x = np.linspace(0, 1, n_rings)
        else:
            x = (np.arange(1, n_rings + 1) - 0.5) * ring_width_um

        peak = float(m_bf.max()) or 1.0
        y  = m_bf / peak
        ye = s_bf / peak

        col = colors[tp]
        ax.fill_between(x, y - ye, y + ye, color=col, alpha=0.20, linewidth=0)
        ax.plot(x, y, color=col, linewidth=2.0,
                label=f"{labels[tp]} (n={n_rings}r)")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.legend(fontsize=6, frameon=False)
    xlabel = ("Fractional position (0=basal, 1=lumen)"
              if normalize_x else "Distance from basal border (µm)")
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel("Fraction of max CPM  (±SEM)", fontsize=7)


def plot_profiles(all_data_real, synth_data, genes, ring_width_um,
                  width_um, width_rings):
    """3-row ribbon plot per gene: real | synthetic µm | synthetic fractional."""
    n = len(genes)
    fig, axes = plt.subplots(3, n, figsize=(4.0 * n, 4.0 * 3), squeeze=False)

    for gi, gene in enumerate(genes):
        # prepare per-gene data bundles for ribbon_ax
        def bundle(src_dict, tps):
            out = {}
            for tp in tps:
                if gene not in src_dict[tp]["means"]:
                    continue
                m = src_dict[tp]["means"][gene]
                s = src_dict[tp]["stds"][gene]
                n_arr = src_dict[tp]["ns"][gene]
                out[tp] = {"mean": m,
                            "sem":  s / np.sqrt(np.maximum(n_arr, 1))}
            return out

        real_b  = bundle(all_data_real, TIMEPOINTS)
        synth_b = bundle(synth_data, SYNTH_TPS)

        for tp, d in real_b.items():
            axes[0][gi]._ribbon_data = real_b
            break

        # Row 0: real data, µm axis
        _ribbon(axes[0][gi], real_b,  TIMEPOINTS,
                COLORS_REAL, TP_LABEL_REAL, ring_width_um,
                normalize_x=False, title=gene)

        d0_max = float(synth_b["day0"]["mean"].max()) if "day0" in synth_b else None

        # Row 1: synthetic, µm axis
        _ribbon(axes[1][gi], synth_b, SYNTH_TPS,
                COLORS_SYNTH, TP_LABEL_SYNTH, ring_width_um,
                normalize_x=False, title=gene, shared_max=d0_max)

        # Row 2: synthetic, fractional axis
        _ribbon(axes[2][gi], synth_b, SYNTH_TPS,
                COLORS_SYNTH, TP_LABEL_SYNTH, ring_width_um,
                normalize_x=True, title=gene, shared_max=d0_max)

    row_labels = [
        "REAL data  (µm x-axis)",
        "SYNTHETIC truncated  (µm x-axis)",
        "SYNTHETIC truncated  (fractional x-axis)",
    ]
    for ri, rl in enumerate(row_labels):
        axes[ri][0].set_ylabel(f"Mean CPM  (±SEM)\n[{rl}]", fontsize=7)
    for ri in range(3):
        for ci in range(1, n):
            leg = axes[ri][ci].get_legend()
            if leg:
                leg.remove()

    _w_str = "  ".join(f"{tp}={width_um[tp]:.0f}µm ({width_rings[tp]}r)"
                       for tp in TIMEPOINTS if tp in width_um)
    fig.suptitle(
        f"{REP} — scale invariance validation  |  REAL vs SYNTHETIC (truncated)\n{_w_str}",
        fontsize=8, fontweight="bold", y=1.01,
    )
    fig.tight_layout()
    return fig


def _ribbon(ax, bundle, tps, colors, labels, ring_width_um,
            normalize_x=False, title="", shared_max=None):
    """
    Draw ribbon profiles from a {tp: {mean, sem}} bundle using raw CPM values.
    shared_max is ignored (kept for API compatibility).
    """
    for tp in tps:
        if tp not in bundle:
            continue
        m_bf = bundle[tp]["mean"][::-1]
        s_bf = bundle[tp]["sem"][::-1]
        n_r  = len(m_bf)

        x = (np.linspace(0, 1, n_r) if normalize_x
             else (np.arange(1, n_r + 1) - 0.5) * ring_width_um)

        col = colors[tp]
        ax.fill_between(x, m_bf - s_bf, m_bf + s_bf, color=col, alpha=0.20, linewidth=0)
        ax.plot(x, m_bf, color=col, linewidth=2.0,
                label=f"{labels[tp]} (n={n_r}r)")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.legend(fontsize=6, frameon=False)
    ax.set_xlabel(
        "Fractional position (0=basal, 1=lumen)" if normalize_x
        else "Distance from basal border (µm)", fontsize=7)
    ax.set_ylabel("Mean CPM  (±SEM)", fontsize=7)


def plot_score_comparison(real_df, synth_df, genes_ordered, method):
    """
    Side-by-side bar chart for one method ('delta_r' or 'delta_rmse').
    method: 'r' or 'rmse'
    """
    dk  = f"delta_{method}"
    lok = f"d{method}_CI95_lo"
    hik = f"d{method}_CI95_hi"

    rm  = real_df.set_index("gene")
    sm  = synth_df.set_index("gene")
    genes = [g for g in genes_ordered if g in rm.index and g in sm.index]

    dr = np.array([rm.loc[g, dk]  for g in genes])
    ds = np.array([sm.loc[g, dk]  for g in genes])
    lr = np.array([rm.loc[g, lok] for g in genes])
    hr = np.array([rm.loc[g, hik] for g in genes])
    ls = np.array([sm.loc[g, lok] for g in genes])
    hs = np.array([sm.loc[g, hik] for g in genes])

    x, w = np.arange(len(genes)), 0.35
    fig, ax = plt.subplots(figsize=(max(7, 0.9 * len(genes)), 4))

    ax.bar(x - w/2, dr, w, label="Real (day0/3/5)",         color="#4477AA", alpha=0.85)
    ax.bar(x + w/2, ds, w, label="Synthetic (truncated)",   color="#EE6677", alpha=0.85)
    ax.errorbar(x - w/2, dr,
                yerr=[np.clip(dr - lr, 0, None), np.clip(hr - dr, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)
    ax.errorbar(x + w/2, ds,
                yerr=[np.clip(ds - ls, 0, None), np.clip(hs - ds, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)

    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=40, ha="right")

    label = "Δr = r_frac − r_raw" if method == "r" else "Δ(RMSE) = RMSE_µm − RMSE_frac"
    ax.set_ylabel(f"{label}  (±95% CI)")
    ax.set_title(
        f"{REP} — {label}\n"
        f"Synthetic = day0 truncated to day3/day5 widths (NOT scale-invariant by construction)\n"
        "FP! = synthetic CI entirely > 0 (false positive)",
        fontsize=8,
    )
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    n_fp = 0
    for gi, gene in enumerate(genes):
        if ls[gi] > 0:
            n_fp += 1
            ax.text(x[gi] + w/2, hs[gi] + abs(hs[gi]) * 0.05 + 0.002,
                    "FP!", ha="center", va="bottom",
                    fontsize=7, color="red", fontweight="bold")

    fig.tight_layout()
    return fig, n_fp, [genes[i] for i in range(len(genes)) if ls[i] > 0]


# ── MAIN ──────────────────────────────────────────────────────────────────────

print(f"\n{'='*60}")
print(f"Scale invariance validation — {REP}")
print(f"{'='*60}")

UM_PER_PX     = load_um_per_hires_px()
ring_width_um = (KERNEL_SIZE / 2) * UM_PER_PX
print(f"µm/hires-px : {UM_PER_PX:.4f}")
print(f"Ring width  : {ring_width_um:.1f} µm")

# Load real data
print("\nLoading real ring profiles …")
all_data_real = {}
found_genes   = set(GENES_10)
for day in TIMEPOINTS:
    print(f"  {day} …")
    _, means, stds, ns, found = load_ring_profiles(day, GENES_10)
    all_data_real[day] = {"means": means, "stds": stds, "ns": ns}
    found_genes &= set(found)

GENES_USE = [g for g in GENES_10 if g in found_genes]
print(f"\nGenes available in all timepoints: {GENES_USE}")

# Width stats
print("\nLoading width stats …")
width_rings, width_um_d = {}, {}
for day in TIMEPOINTS:
    n_r, m_um = load_mean_width_rings(day, ring_width_um)
    width_rings[day] = n_r
    width_um_d[day]  = m_um
    n_actual = len(all_data_real[day]["means"].get(GENES_USE[0], []))
    print(f"  {day}: {m_um:.0f} µm → {n_r} rings target  (actual: {n_actual})")

# Build synthetic data
print("\nBuilding synthetic truncated data …")
_synth_map = {"day0": "day0", "day3": "day3_trunc", "day5": "day5_trunc"}
SYNTH_TPS  = [_synth_map[tp] for tp in TIMEPOINTS]
synth_data = {"day0": all_data_real["day0"]}
for tp in TIMEPOINTS:
    if tp != "day0":
        key = _synth_map[tp]
        synth_data[key] = make_truncated_data(all_data_real["day0"], width_rings[tp])
for key in SYNTH_TPS:
    n = len(synth_data[key]["means"].get(GENES_USE[0], []))
    print(f"  {key}: {n} rings")

# Score real data
print(f"\nScoring REAL data ({N_BOOT} bootstrap iterations) …")
real_df, real_pairs = score_both_methods(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)
real_euc  = score_euclidean_dist(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)
real_zpd  = score_euclidean_zeropad(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)

print("\nReal data — Euclidean distance scores (fractional; lower = more scale-invariant):")
print(real_euc.to_string(index=False))
print("\nReal data — Zero-padded Euclidean distance (absolute space; higher = more scale-invariant):")
print(real_zpd.to_string(index=False))

print("\nReal data — Δr / Δ(RMSE) scores:")
cols = ["gene", "delta_r", "dr_CI95_lo", "dr_CI95_hi",
        "delta_rmse", "drmse_CI95_lo", "drmse_CI95_hi"]
print(real_df[cols].to_string(index=False))

# Score real data — ΔI and ΔI_alt
real_dI_alt = score_delta_I_alt(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)
print("\nReal data — ΔI_alt (mean-normalised; lower = more scale-invariant):")
print(real_dI_alt.to_string(index=False))

real_dI = score_delta_I(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)
print("\nReal data — ΔI (lower = more scale-invariant):")
print(real_dI.to_string(index=False))

# Score real data — normalized fractional method
real_nfr = score_normdist_fractional(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)
print("\nReal data — Normalized distance fractional x (lower = more scale-invariant):")
print(real_nfr.to_string(index=False))

# Score real data — normalized method
real_nrm = score_euclidean_normalized(
    all_data_real, GENES_USE, TIMEPOINTS, ring_width_um, n_boot=N_BOOT)
print("\nReal data — Normalized Euclidean (pointwise; higher = more scale-invariant):")
print(real_nrm.to_string(index=False))

# Score synthetic data
print(f"\nScoring SYNTHETIC data ({N_BOOT} bootstrap iterations) …")
synth_df, synth_pairs = score_both_methods(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)
synth_euc = score_euclidean_dist(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)
synth_zpd = score_euclidean_zeropad(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)
synth_dI_alt = score_delta_I_alt(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)
print("\nSynthetic data — ΔI_alt:")
print(synth_dI_alt.to_string(index=False))

synth_dI = score_delta_I(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)
print("\nSynthetic data — ΔI:")
print(synth_dI.to_string(index=False))

synth_nfr = score_normdist_fractional(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)
synth_nrm = score_euclidean_normalized(
    synth_data, GENES_USE, SYNTH_TPS, ring_width_um, n_boot=N_BOOT)

print("\nSynthetic data — Euclidean distance (fractional):")
print(synth_euc.to_string(index=False))
print("\nSynthetic data — Zero-padded Euclidean distance (absolute space):")
print(synth_zpd.to_string(index=False))
print("\nSynthetic data — Normalized distance fractional x:")
print(synth_nfr.to_string(index=False))
print("\nSynthetic data — Normalized Euclidean (zero-padded absolute):")
print(synth_nrm.to_string(index=False))

# Save CSVs
OUT_DIR.mkdir(parents=True, exist_ok=True)
real_df.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_real_scores.csv",  index=False)
synth_df.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_scores.csv", index=False)
real_euc.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_real_euclidean.csv",  index=False)
synth_euc.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_euclidean.csv", index=False)
real_zpd.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_real_zeropad.csv",  index=False)
synth_zpd.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_zeropad.csv", index=False)
real_dI.to_csv(OUT_DIR      / f"{REP}_{TP_SUFFIX}_real_deltaI.csv",       index=False)
synth_dI.to_csv(OUT_DIR    / f"{REP}_{TP_SUFFIX}_synth_deltaI.csv",      index=False)
real_dI_alt.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_real_deltaI_alt.csv",   index=False)
synth_dI_alt.to_csv(OUT_DIR/ f"{REP}_{TP_SUFFIX}_synth_deltaI_alt.csv",  index=False)
real_nrm.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_real_normdist.csv",  index=False)
synth_nrm.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_normdist.csv", index=False)
real_nfr.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_real_normdist_frac.csv",  index=False)
synth_nfr.to_csv(OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_normdist_frac.csv", index=False)
print(f"\nCSVs saved to {OUT_DIR}")

# ── Normdist-frac input curves ─────────────────────────────────────────────
# Shows exactly the curves compared in score_normdist_fractional:
# each profile mapped to its own [0,1] fractional axis, raw CPM.
# Row 0 = real data; Row 1 = synthetic truncated data.
print("\nPlotting normdist-frac input curves …")
n_g = len(GENES_USE)
fig_fc, axes_fc = plt.subplots(2, n_g, figsize=(4.0 * n_g, 4.0 * 2), squeeze=False)

for gi, gene in enumerate(GENES_USE):
    for row, (src, tps, colors, labels) in enumerate([
        (all_data_real, TIMEPOINTS,
         COLORS_REAL,  TP_LABEL_REAL),
        (synth_data,   SYNTH_TPS,
         COLORS_SYNTH, TP_LABEL_SYNTH),
    ]):
        ax = axes_fc[row][gi]
        for tp in tps:
            m = src[tp]["means"].get(gene)
            s = src[tp]["stds"].get(gene)
            n_arr = src[tp]["ns"].get(gene)
            if m is None:
                continue
            m_bf  = m[::-1]
            sem   = s[::-1] / np.sqrt(np.maximum(n_arr[::-1], 1))
            n_r   = len(m_bf)
            x     = np.linspace(0, 1, n_r)   # each profile to its own [0,1]
            col   = colors[tp]
            lbl   = labels[tp]
            ax.fill_between(x, m_bf - sem, m_bf + sem, color=col, alpha=0.20, linewidth=0)
            ax.plot(x, m_bf, color=col, linewidth=2.0, label=f"{lbl} ({n_r}r)")
        ax.set_title(gene, fontsize=8, fontweight="bold")
        ax.set_xlabel("Fractional position (0=basal, 1=lumen)", fontsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if gi == 0:
            ax.set_ylabel("Mean CPM  (±SEM)", fontsize=7)
        leg = ax.legend(fontsize=6, frameon=False)
        if gi > 0 and leg:
            leg.remove()

axes_fc[0][0].set_ylabel("REAL — Mean CPM  (±SEM)", fontsize=7)
axes_fc[1][0].set_ylabel("SYNTHETIC — Mean CPM  (±SEM)", fontsize=7)

fig_fc.suptitle(
    f"{REP} — normdist_frac input curves  |  each profile mapped to own [0,1], raw CPM\n"
    "Score = mean (f−g)² / ((f+g)/2)²  at each fractional position  (lower = more scale-invariant)",
    fontsize=8, fontweight="bold", y=1.01,
)
fig_fc.tight_layout()
fc_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_normdist_frac_curves.svg"
fig_fc.savefig(fc_path, bbox_inches="tight", format="svg")
plt.close(fig_fc)
print(f"  Saved: {fc_path}")

# ── Synthetic verification plot ────────────────────────────────────────────
# Shows synthetic profiles on a SHARED absolute µm x-axis, normalized by
# day0's own max so all three curves should overlap and simply end at
# different x positions — confirming the truncation was done correctly.
print("\nPlotting synthetic verification …")
n_g = len(GENES_USE)
fig_v, axes_v = plt.subplots(1, n_g, figsize=(4.0 * n_g, 4.0), squeeze=False)

for gi, gene in enumerate(GENES_USE):
    ax = axes_v[0][gi]

    _verif_tps = [("day0", COLORS_SYNTH["day0"], "day0 (full)")] + [
        (stp, COLORS_SYNTH[stp], f"day0 → {tp} width ({width_rings[tp]}r)")
        for tp, stp in _synth_map.items() if tp != "day0" and stp in synth_data
    ]
    for tp, col, lbl in _verif_tps:
        m = synth_data[tp]["means"].get(gene)
        s = synth_data[tp]["stds"].get(gene)
        n_arr = synth_data[tp]["ns"].get(gene)
        if m is None:
            continue
        m_bf  = m[::-1]           # basal-first
        sem   = (s[::-1] / np.sqrt(np.maximum(n_arr[::-1], 1)))
        n_r   = len(m_bf)
        x     = (np.arange(1, n_r + 1) - 0.5) * ring_width_um

        ax.fill_between(x, m_bf - sem, m_bf + sem, color=col, alpha=0.20, linewidth=0)
        ax.plot(x, m_bf, color=col, linewidth=2.0, label=lbl)

    ax.set_title(gene, fontsize=8, fontweight="bold")
    ax.set_xlabel("Distance from basal border (µm)", fontsize=7)
    ax.set_ylabel("Mean CPM  (±SEM)" if gi == 0 else "", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    leg = ax.legend(fontsize=6, frameon=False)
    if gi > 0 and leg:
        leg.remove()

fig_v.suptitle(
    f"{REP} — synthetic truncation verification (raw CPM)\n"
    "All curves = same day0 data, truncated → they should OVERLAP and simply END at different x positions\n"
    "  ".join(f"{tp}={width_um_d[tp]:.0f}µm" for tp in TIMEPOINTS if tp in width_um_d),
    fontsize=8, fontweight="bold", y=1.02,
)
fig_v.tight_layout()
verif_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_verification.svg"
fig_v.savefig(verif_path, bbox_inches="tight", format="svg")
plt.close(fig_v)
print(f"  Saved: {verif_path}")

# Plot profiles
print("\nPlotting profiles …")
fig_prof = plot_profiles(all_data_real, synth_data, GENES_USE,
                         ring_width_um, width_um_d, width_rings)
prof_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_synth_profiles.svg"
fig_prof.savefig(prof_path, bbox_inches="tight", format="svg")
plt.close(fig_prof)
print(f"  Saved: {prof_path}")

# Score comparison plots — Δr and Δ(RMSE)
print("Plotting score comparisons …")
genes_ordered = real_df["gene"].tolist()   # sorted by delta_rmse descending

for method, label in [("r", "delta_r"), ("rmse", "delta_rmse")]:
    fig_sc, n_fp, fp_genes = plot_score_comparison(
        real_df, synth_df, genes_ordered, method)
    sc_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_{label}.svg"
    fig_sc.savefig(sc_path, bbox_inches="tight", format="svg")
    plt.close(fig_sc)
    print(f"  [{label}]  FP: {n_fp}/{len(genes_ordered)}  "
          f"{'→ ' + str(fp_genes) if fp_genes else '→ none ✓'}  |  {sc_path.name}")

# Euclidean distance comparison plot
# genes ordered by real euc_dist ascending (most scale-invariant first)
genes_euc = real_euc["gene"].tolist()
rm_e  = real_euc.set_index("gene")
sm_e  = synth_euc.set_index("gene")
genes_euc = [g for g in genes_euc if g in sm_e.index]

dr_e = np.array([rm_e.loc[g, "euc_dist"] for g in genes_euc])
ds_e = np.array([sm_e.loc[g, "euc_dist"] for g in genes_euc])
lr_e = np.array([rm_e.loc[g, "CI95_lo"]  for g in genes_euc])
hr_e = np.array([rm_e.loc[g, "CI95_hi"]  for g in genes_euc])
ls_e = np.array([sm_e.loc[g, "CI95_lo"]  for g in genes_euc])
hs_e = np.array([sm_e.loc[g, "CI95_hi"]  for g in genes_euc])

x_e, w_e = np.arange(len(genes_euc)), 0.35
fig_euc, ax_euc = plt.subplots(figsize=(max(7, 0.9 * len(genes_euc)), 4.5))
ax_euc.bar(x_e - w_e/2, dr_e, w_e, label="Real (day0/3/5)",       color="#4477AA", alpha=0.85)
ax_euc.bar(x_e + w_e/2, ds_e, w_e, label="Synthetic (truncated)", color="#EE6677", alpha=0.85)
ax_euc.errorbar(x_e - w_e/2, dr_e,
                yerr=[np.clip(dr_e - lr_e, 0, None), np.clip(hr_e - dr_e, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)
ax_euc.errorbar(x_e + w_e/2, ds_e,
                yerr=[np.clip(ds_e - ls_e, 0, None), np.clip(hs_e - ds_e, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)

ax_euc.set_xticks(x_e)
ax_euc.set_xticklabels(genes_euc, rotation=40, ha="right")
ax_euc.set_ylabel("Mean pairwise Euclidean distance in fractional space\n(raw CPM; lower = more scale-invariant)")
ax_euc.set_title(
    f"{REP} — Euclidean distance score (lower = scale-invariant)\n"
    "Synthetic = day0 truncated to day3/day5 widths (not scale-invariant by construction)\n"
    "FP! = synthetic score ≤ real score (method cannot discriminate)",
    fontsize=8,
)
ax_euc.legend(fontsize=8, frameon=False)
ax_euc.spines["top"].set_visible(False)
ax_euc.spines["right"].set_visible(False)

# Mark false positives: synthetic CI_hi ≤ real euc_dist (synth looks as scale-invariant as real)
fp_euc = []
for gi, gene in enumerate(genes_euc):
    if hs_e[gi] <= dr_e[gi]:   # synth upper CI ≤ real point estimate
        fp_euc.append(gene)
        ax_euc.text(x_e[gi] + w_e/2, hs_e[gi] + 0.002,
                    "FP!", ha="center", va="bottom",
                    fontsize=7, color="red", fontweight="bold")

fig_euc.tight_layout()
euc_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_euclidean.svg"
fig_euc.savefig(euc_path, bbox_inches="tight", format="svg")
plt.close(fig_euc)
print(f"  [euclidean]  FP: {len(fp_euc)}/{len(genes_euc)}  "
      f"{'→ ' + str(fp_euc) if fp_euc else '→ none ✓'}  |  {euc_path.name}")

# Zero-pad comparison plot
# For zero-pad: higher = more scale-invariant.
# FP = synth score ≥ real score (method can't tell apart)
genes_zpd = real_zpd["gene"].tolist()
rm_z  = real_zpd.set_index("gene")
sm_z  = synth_zpd.set_index("gene")
genes_zpd = [g for g in genes_zpd if g in sm_z.index]

dr_z = np.array([rm_z.loc[g, "euc_zpad"] for g in genes_zpd])
ds_z = np.array([sm_z.loc[g, "euc_zpad"] for g in genes_zpd])
lr_z = np.array([rm_z.loc[g, "CI95_lo"]  for g in genes_zpd])
hr_z = np.array([rm_z.loc[g, "CI95_hi"]  for g in genes_zpd])
ls_z = np.array([sm_z.loc[g, "CI95_lo"]  for g in genes_zpd])
hs_z = np.array([sm_z.loc[g, "CI95_hi"]  for g in genes_zpd])

x_z, w_z = np.arange(len(genes_zpd)), 0.35
fig_zpd, ax_zpd = plt.subplots(figsize=(max(7, 0.9 * len(genes_zpd)), 4.5))
ax_zpd.bar(x_z - w_z/2, dr_z, w_z, label="Real (day0/3/5)",       color="#4477AA", alpha=0.85)
ax_zpd.bar(x_z + w_z/2, ds_z, w_z, label="Synthetic (truncated)", color="#EE6677", alpha=0.85)
ax_zpd.errorbar(x_z - w_z/2, dr_z,
                yerr=[np.clip(dr_z - lr_z, 0, None), np.clip(hr_z - dr_z, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)
ax_zpd.errorbar(x_z + w_z/2, ds_z,
                yerr=[np.clip(ds_z - ls_z, 0, None), np.clip(hs_z - ds_z, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)

ax_zpd.set_xticks(x_z)
ax_zpd.set_xticklabels(genes_zpd, rotation=40, ha="right")
ax_zpd.set_ylabel("Mean pairwise Euclidean distance in absolute ring space\n(raw CPM, zero-padded; higher = more scale-invariant)")
ax_zpd.set_title(
    f"{REP} — Zero-padded Euclidean distance (higher = scale-invariant)\n"
    "Profiles aligned at basal end; lumenal overhang of shorter timepoints = 0\n"
    "FP! = synthetic score ≥ real score CI_lo (method cannot discriminate)",
    fontsize=8,
)
ax_zpd.legend(fontsize=8, frameon=False)
ax_zpd.spines["top"].set_visible(False)
ax_zpd.spines["right"].set_visible(False)

fp_zpd = []
for gi, gene in enumerate(genes_zpd):
    # FP if synthetic lower CI ≥ real point estimate (synth at least as "scale-invariant")
    if ls_z[gi] >= dr_z[gi]:
        fp_zpd.append(gene)
        ax_zpd.text(x_z[gi] + w_z/2, hs_z[gi] + 0.002,
                    "FP!", ha="center", va="bottom",
                    fontsize=7, color="red", fontweight="bold")

fig_zpd.tight_layout()
zpd_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_zeropad.svg"
fig_zpd.savefig(zpd_path, bbox_inches="tight", format="svg")
plt.close(fig_zpd)
print(f"  [zeropad]    FP: {len(fp_zpd)}/{len(genes_zpd)}  "
      f"{'→ ' + str(fp_zpd) if fp_zpd else '→ none ✓'}  |  {zpd_path.name}")

# Normalized Euclidean comparison plot
genes_nrm = real_nrm["gene"].tolist()
rm_n  = real_nrm.set_index("gene")
sm_n  = synth_nrm.set_index("gene")
genes_nrm = [g for g in genes_nrm if g in sm_n.index]

dr_n = np.array([rm_n.loc[g, "norm_dist"] for g in genes_nrm])
ds_n = np.array([sm_n.loc[g, "norm_dist"] for g in genes_nrm])
lr_n = np.array([rm_n.loc[g, "CI95_lo"]   for g in genes_nrm])
hr_n = np.array([rm_n.loc[g, "CI95_hi"]   for g in genes_nrm])
ls_n = np.array([sm_n.loc[g, "CI95_lo"]   for g in genes_nrm])
hs_n = np.array([sm_n.loc[g, "CI95_hi"]   for g in genes_nrm])

x_n, w_n = np.arange(len(genes_nrm)), 0.35
fig_nrm, ax_nrm = plt.subplots(figsize=(max(7, 0.9 * len(genes_nrm)), 4.5))
ax_nrm.bar(x_n - w_n/2, dr_n, w_n, label="Real (day0/3/5)",       color="#4477AA", alpha=0.85)
ax_nrm.bar(x_n + w_n/2, ds_n, w_n, label="Synthetic (truncated)", color="#EE6677", alpha=0.85)
ax_nrm.errorbar(x_n - w_n/2, dr_n,
                yerr=[np.clip(dr_n - lr_n, 0, None), np.clip(hr_n - dr_n, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)
ax_nrm.errorbar(x_n + w_n/2, ds_n,
                yerr=[np.clip(ds_n - ls_n, 0, None), np.clip(hs_n - ds_n, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)

ax_nrm.set_xticks(x_n)
ax_nrm.set_xticklabels(genes_nrm, rotation=40, ha="right")
ax_nrm.set_ylabel("Mean (f−g)² / ((f+g)/2)²  (zero-padded; higher = more scale-invariant)")
ax_nrm.set_title(
    f"{REP} — Normalized Euclidean distance (higher = scale-invariant)\n"
    "Pointwise (f−g)²/((f+g)/2)², zero-padded absolute space\n"
    "FP! = synthetic CI_lo ≥ real point estimate",
    fontsize=8,
)
ax_nrm.legend(fontsize=8, frameon=False)
ax_nrm.spines["top"].set_visible(False)
ax_nrm.spines["right"].set_visible(False)

fp_nrm = []
for gi, gene in enumerate(genes_nrm):
    if ls_n[gi] >= dr_n[gi]:
        fp_nrm.append(gene)
        ax_nrm.text(x_n[gi] + w_n/2, hs_n[gi] + abs(hs_n[gi]) * 0.05 + 0.002,
                    "FP!", ha="center", va="bottom",
                    fontsize=7, color="red", fontweight="bold")

fig_nrm.tight_layout()
nrm_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_normdist.svg"
fig_nrm.savefig(nrm_path, bbox_inches="tight", format="svg")
plt.close(fig_nrm)
print(f"  [norm_dist]  FP: {len(fp_nrm)}/{len(genes_nrm)}  "
      f"{'→ ' + str(fp_nrm) if fp_nrm else '→ none ✓'}  |  {nrm_path.name}")

# Normalized fractional comparison plot (lower = SI)
genes_nfr = real_nfr["gene"].tolist()
rm_nf  = real_nfr.set_index("gene")
sm_nf  = synth_nfr.set_index("gene")
genes_nfr = [g for g in genes_nfr if g in sm_nf.index]

dr_nf = np.array([rm_nf.loc[g, "normdist_frac"] for g in genes_nfr])
ds_nf = np.array([sm_nf.loc[g, "normdist_frac"] for g in genes_nfr])
lr_nf = np.array([rm_nf.loc[g, "CI95_lo"] for g in genes_nfr])
hr_nf = np.array([rm_nf.loc[g, "CI95_hi"] for g in genes_nfr])
ls_nf = np.array([sm_nf.loc[g, "CI95_lo"] for g in genes_nfr])
hs_nf = np.array([sm_nf.loc[g, "CI95_hi"] for g in genes_nfr])

x_nf, w_nf = np.arange(len(genes_nfr)), 0.35
fig_nfr, ax_nfr = plt.subplots(figsize=(max(7, 0.9 * len(genes_nfr)), 4.5))
ax_nfr.bar(x_nf - w_nf/2, dr_nf, w_nf, label="Real (day0/3/5)",       color="#4477AA", alpha=0.85)
ax_nfr.bar(x_nf + w_nf/2, ds_nf, w_nf, label="Synthetic (truncated)", color="#EE6677", alpha=0.85)
ax_nfr.errorbar(x_nf - w_nf/2, dr_nf,
                yerr=[np.clip(dr_nf - lr_nf, 0, None), np.clip(hr_nf - dr_nf, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)
ax_nfr.errorbar(x_nf + w_nf/2, ds_nf,
                yerr=[np.clip(ds_nf - ls_nf, 0, None), np.clip(hs_nf - ds_nf, 0, None)],
                fmt="none", color="black", capsize=3, linewidth=1.0)

ax_nfr.set_xticks(x_nf)
ax_nfr.set_xticklabels(genes_nfr, rotation=40, ha="right")
ax_nfr.set_ylabel("Mean (f−g)² / ((f+g)/2)²  (fractional x; lower = more scale-invariant)")
ax_nfr.set_title(
    f"{REP} — Normalized distance, fractional x (lower = scale-invariant)\n"
    "Each profile mapped to own [0,1]; raw CPM; no zero-padding\n"
    "FP! = synthetic CI_hi ≥ real point estimate",
    fontsize=8,
)
ax_nfr.legend(fontsize=8, frameon=False)
ax_nfr.spines["top"].set_visible(False)
ax_nfr.spines["right"].set_visible(False)

fp_nfr = []
for gi, gene in enumerate(genes_nfr):
    if ds_nf[gi] <= dr_nf[gi]:   # synth point est ≤ real (synth looks as SI as real)
        fp_nfr.append(gene)
        ax_nfr.text(x_nf[gi] + w_nf/2, hs_nf[gi] + abs(hs_nf[gi]) * 0.05 + 0.002,
                    "FP!", ha="center", va="bottom",
                    fontsize=7, color="red", fontweight="bold")

fig_nfr.tight_layout()
nfr_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_normdist_frac.svg"
fig_nfr.savefig(nfr_path, bbox_inches="tight", format="svg")
plt.close(fig_nfr)
print(f"  [normdist_frac] FP: {len(fp_nfr)}/{len(genes_nfr)}  "
      f"{'→ ' + str(fp_nfr) if fp_nfr else '→ none ✓'}  |  {nfr_path.name}")

# ΔI comparison plot (lower = SI)
genes_dI = real_dI["gene"].tolist()
rm_dI  = real_dI.set_index("gene")
sm_dI  = synth_dI.set_index("gene")
genes_dI = [g for g in genes_dI if g in sm_dI.index]

dr_dI = np.array([rm_dI.loc[g, "delta_I"] for g in genes_dI])
ds_dI = np.array([sm_dI.loc[g, "delta_I"] for g in genes_dI])
lr_dI = np.array([rm_dI.loc[g, "CI95_lo"] for g in genes_dI])
hr_dI = np.array([rm_dI.loc[g, "CI95_hi"] for g in genes_dI])
ls_dI = np.array([sm_dI.loc[g, "CI95_lo"] for g in genes_dI])
hs_dI = np.array([sm_dI.loc[g, "CI95_hi"] for g in genes_dI])

x_dI, w_dI = np.arange(len(genes_dI)), 0.35
fig_dI, ax_dI = plt.subplots(figsize=(max(7, 0.9 * len(genes_dI)), 4.5))
ax_dI.bar(x_dI - w_dI/2, dr_dI, w_dI, label="Real (day0/3/5)",       color="#4477AA", alpha=0.85)
ax_dI.bar(x_dI + w_dI/2, ds_dI, w_dI, label="Synthetic (truncated)", color="#EE6677", alpha=0.85)
ax_dI.errorbar(x_dI - w_dI/2, dr_dI,
               yerr=[np.clip(dr_dI - lr_dI, 0, None), np.clip(hr_dI - dr_dI, 0, None)],
               fmt="none", color="black", capsize=3, linewidth=1.0)
ax_dI.errorbar(x_dI + w_dI/2, ds_dI,
               yerr=[np.clip(ds_dI - ls_dI, 0, None), np.clip(hs_dI - ds_dI, 0, None)],
               fmt="none", color="black", capsize=3, linewidth=1.0)

ax_dI.set_xticks(x_dI)
ax_dI.set_xticklabels(genes_dI, rotation=40, ha="right")
ax_dI.set_ylabel("ΔI = ½⟨log₂(1 + σ²_between/σ²_meas)⟩  (lower = more scale-invariant)")
ax_dI.set_title(
    f"{REP} — ΔI score (Nikolic et al. 2024 adapted; lower = scale-invariant)\n"
    "Fractional x, raw CPM; σ²_between = between-timepoint variance at each x_s\n"
    "FP! = synthetic CI_hi ≥ real point estimate",
    fontsize=8,
)
ax_dI.legend(fontsize=8, frameon=False)
ax_dI.spines["top"].set_visible(False)
ax_dI.spines["right"].set_visible(False)

fp_dI = []
for gi, gene in enumerate(genes_dI):
    if ds_dI[gi] <= dr_dI[gi]:   # synth point est ≤ real (synth looks as SI as real)
        fp_dI.append(gene)
        ax_dI.text(x_dI[gi] + w_dI/2, hs_dI[gi] + abs(hs_dI[gi]) * 0.05 + 0.002,
                   "FP!", ha="center", va="bottom",
                   fontsize=7, color="red", fontweight="bold")

fig_dI.tight_layout()
dI_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_deltaI.svg"
fig_dI.savefig(dI_path, bbox_inches="tight", format="svg")
plt.close(fig_dI)
print(f"  [delta_I]       FP: {len(fp_dI)}/{len(genes_dI)}  "
      f"{'→ ' + str(fp_dI) if fp_dI else '→ none ✓'}  |  {dI_path.name}")

# ΔI_alt comparison plot
genes_da = real_dI_alt["gene"].tolist()
rm_da = real_dI_alt.set_index("gene")
sm_da = synth_dI_alt.set_index("gene")
genes_da = [g for g in genes_da if g in sm_da.index]

dr_da = np.array([rm_da.loc[g, "dI_alt"] for g in genes_da])
ds_da = np.array([sm_da.loc[g, "dI_alt"] for g in genes_da])
lr_da = np.array([rm_da.loc[g, "CI95_lo"] for g in genes_da])
hr_da = np.array([rm_da.loc[g, "CI95_hi"] for g in genes_da])
ls_da = np.array([sm_da.loc[g, "CI95_lo"] for g in genes_da])
hs_da = np.array([sm_da.loc[g, "CI95_hi"] for g in genes_da])

x_da, w_da = np.arange(len(genes_da)), 0.35
fig_da, ax_da = plt.subplots(figsize=(max(7, 0.9 * len(genes_da)), 4.5))
ax_da.bar(x_da - w_da/2, dr_da, w_da, label="Real (day0/3/5)",       color="#4477AA", alpha=0.85)
ax_da.bar(x_da + w_da/2, ds_da, w_da, label="Synthetic (truncated)", color="#EE6677", alpha=0.85)
ax_da.errorbar(x_da - w_da/2, dr_da,
               yerr=[np.clip(dr_da - lr_da, 0, None), np.clip(hr_da - dr_da, 0, None)],
               fmt="none", color="black", capsize=3, linewidth=1.0)
ax_da.errorbar(x_da + w_da/2, ds_da,
               yerr=[np.clip(ds_da - ls_da, 0, None), np.clip(hs_da - ds_da, 0, None)],
               fmt="none", color="black", capsize=3, linewidth=1.0)
ax_da.set_xticks(x_da)
ax_da.set_xticklabels(genes_da, rotation=40, ha="right")
ax_da.set_ylabel("ΔI_alt = ½⟨log₂(1 + σ²_between/μ²)⟩  (lower = more scale-invariant)")
ax_da.set_title(
    f"{REP} — ΔI_alt (mean-normalised; lower = scale-invariant)\n"
    "Fractional x, raw CPM; normalised by local mean² (not SEM²)\n"
    "FP! = synthetic CI_hi ≥ real point estimate",
    fontsize=8,
)
ax_da.legend(fontsize=8, frameon=False)
ax_da.spines["top"].set_visible(False)
ax_da.spines["right"].set_visible(False)

fp_da = []
for gi, gene in enumerate(genes_da):
    if ds_da[gi] <= dr_da[gi]:   # synth point est ≤ real (synth looks as SI as real)
        fp_da.append(gene)
        ax_da.text(x_da[gi] + w_da/2, hs_da[gi] + abs(hs_da[gi]) * 0.05 + 0.002,
                   "FP!", ha="center", va="bottom",
                   fontsize=7, color="red", fontweight="bold")

fig_da.tight_layout()
da_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_comparison_deltaI_alt.svg"
fig_da.savefig(da_path, bbox_inches="tight", format="svg")
plt.close(fig_da)
print(f"  [dI_alt]        FP: {len(fp_da)}/{len(genes_da)}  "
      f"{'→ ' + str(fp_da) if fp_da else '→ none ✓'}  |  {da_path.name}")

# Summary
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
sm_real  = real_df.set_index("gene")
sm_synth = synth_df.set_index("gene")

for method in ("delta_r", "delta_rmse"):
    lok_key = "dr_CI95_lo" if method == "delta_r" else "drmse_CI95_lo"
    vals_r  = real_df[method].values
    vals_s  = synth_df[method].values
    fp_g    = [g for g in GENES_USE
               if g in sm_synth.index and sm_synth.loc[g, lok_key] > 0]
    print(f"\n{method}:")
    print(f"  Real mean:      {vals_r.mean():.4f}")
    print(f"  Synthetic mean: {vals_s.mean():.4f}")
    print(f"  False positives (synth CI > 0): {len(fp_g)}/{len(GENES_USE)}"
          + (f"  → {fp_g}" if fp_g else "  → none ✓"))

print(f"\neuclidean_dist fractional (lower = scale-invariant):")
print(f"  Real mean:      {dr_e.mean():.5f}")
print(f"  Synthetic mean: {ds_e.mean():.5f}")
print(f"  False positives (synth upper CI ≤ real): {len(fp_euc)}/{len(genes_euc)}"
      + (f"  → {fp_euc}" if fp_euc else "  → none ✓"))

print(f"\neuclidean_dist zero-padded absolute (higher = scale-invariant):")
print(f"  Real mean:      {dr_z.mean():.5f}")
print(f"  Synthetic mean: {ds_z.mean():.5f}")
print(f"  False positives (synth CI_lo ≥ real): {len(fp_zpd)}/{len(genes_zpd)}"
      + (f"  → {fp_zpd}" if fp_zpd else "  → none ✓"))

print(f"\nnormalized Euclidean (f−g)²/((f+g)/2)², zero-padded (higher = scale-invariant):")
print(f"  Real mean:      {dr_n.mean():.5f}")
print(f"  Synthetic mean: {ds_n.mean():.5f}")
print(f"  False positives (synth CI_lo ≥ real): {len(fp_nrm)}/{len(genes_nrm)}"
      + (f"  → {fp_nrm}" if fp_nrm else "  → none ✓"))

print(f"\nnormdist fractional x (f−g)²/((f+g)/2)², each profile to own [0,1] (lower = scale-invariant):")
print(f"  Real mean:      {dr_nf.mean():.5f}")
print(f"  Synthetic mean: {ds_nf.mean():.5f}")
print(f"  False positives (synth point est ≤ real): {len(fp_nfr)}/{len(genes_nfr)}"
      + (f"  → {fp_nfr}" if fp_nfr else "  → none ✓"))

print(f"\nΔI (Nikolic adapted, fractional x, lower = scale-invariant):")
print(f"  Real mean:      {dr_dI.mean():.5f}")
print(f"  Synthetic mean: {ds_dI.mean():.5f}")
print(f"  False positives (synth point est ≤ real): {len(fp_dI)}/{len(genes_dI)}"
      + (f"  → {fp_dI}" if fp_dI else "  → none ✓"))

print(f"\nΔI_alt (mean-normalised, fractional x, lower = scale-invariant):")
print(f"  Real mean:      {dr_da.mean():.5f}")
print(f"  Synthetic mean: {ds_da.mean():.5f}")
print(f"  False positives (synth point est ≤ real): {len(fp_da)}/{len(genes_da)}"
      + (f"  → {fp_da}" if fp_da else "  → none ✓"))

# Summary CSV table
print("\nBuilding summary CSV table …")
_rs = real_df.set_index("gene")
_ss = synth_df.set_index("gene")
_re = real_euc.set_index("gene")
_se = synth_euc.set_index("gene")
_rz = real_zpd.set_index("gene")
_sz = synth_zpd.set_index("gene")
_rn = real_nrm.set_index("gene")
_sn = synth_nrm.set_index("gene")

tbl_rows = []
for g in GENES_USE:
    tbl_rows.append({
        "gene":             g,
        "Δr_real":          round(_rs.loc[g, "delta_r"],    3),
        "Δr_synth":         round(_ss.loc[g, "delta_r"],    3),
        "ΔRMSE_real":       round(_rs.loc[g, "delta_rmse"], 1),
        "ΔRMSE_synth":      round(_ss.loc[g, "delta_rmse"], 1),
        "Euc_frac_real":    round(_re.loc[g, "euc_dist"],   1),
        "Euc_frac_synth":   round(_se.loc[g, "euc_dist"],   1),
        "Euc_zpad_real":    round(_rz.loc[g, "euc_zpad"],   1),
        "Euc_zpad_synth":   round(_sz.loc[g, "euc_zpad"],   1),
        "Norm_dist_real":      round(_rn.loc[g, "norm_dist"],      4),
        "Norm_dist_synth":     round(_sn.loc[g, "norm_dist"],      4),
        "Normdist_frac_real":  round(real_nfr.set_index("gene").loc[g, "normdist_frac"], 4),
        "Normdist_frac_synth": round(synth_nfr.set_index("gene").loc[g, "normdist_frac"], 4),
        "DeltaI_real":         round(real_dI.set_index("gene").loc[g,     "delta_I"], 5),
        "DeltaI_synth":        round(synth_dI.set_index("gene").loc[g,    "delta_I"], 5),
        "DeltaI_alt_real":     round(real_dI_alt.set_index("gene").loc[g,  "dI_alt"], 5),
        "DeltaI_alt_synth":    round(synth_dI_alt.set_index("gene").loc[g, "dI_alt"], 5),
    })
tbl_df = pd.DataFrame(tbl_rows)
tbl_path = OUT_DIR / f"{REP}_{TP_SUFFIX}_score_summary_table.csv"
tbl_df.to_csv(tbl_path, index=False)
print(f"  Saved: {tbl_path.name}")

print(f"\nOutputs in: {OUT_DIR}")
