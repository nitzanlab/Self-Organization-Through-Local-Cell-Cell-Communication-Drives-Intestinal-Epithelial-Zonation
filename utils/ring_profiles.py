"""Ring-profile helpers vendored from the scale-invariance analysis.

Extracted verbatim so the paper repo has no dependency on an external working
repository. Source: monolayer_analysis/crohn's scale invariance/crohn_io.py
"""
import numpy as np
import pandas as pd

def ring_profiles(X_genes_by_bins, barcodes, features_df, ring_df, genes,
                  exclude_last_ring=True, min_bin_counts=1):
    """
    CPM-normalise every bin, then compute mean CPM per ring for each requested gene.
    Returns (ring_ids, means, stds, ns, found_genes, ns_all) ordered lumen->basal.

    min_bin_counts : only bins whose TOTAL counts >= this are averaged into a ring
                     (default 1 = exclude empty bins). This matters near the luminal
                     border: if the drawn border sits slightly off the tissue, the
                     first rings can contain empty background bins that would other-
                     wise dilute the mean. `ns` = bins with counts (used for mean/SEM);
                     `ns_all` = all bins assigned to the ring (so you can see the gap).

    X_genes_by_bins : csr (genes x bins), columns aligned to `barcodes`
    ring_df         : [barcode, ring]  (subset/aligned to `barcodes` beforehand)
    """
    # align ring labels to matrix columns
    ring_by_bc = dict(zip(ring_df["barcode"].astype(str), ring_df["ring"]))
    rings = np.array([ring_by_bc.get(str(b), np.nan) for b in barcodes], dtype=float)
    valid = ~np.isnan(rings)

    # CPM per bin (bins x genes); keep the RAW per-bin totals for the count filter
    X_bg = X_genes_by_bins.T.tocsr().astype(np.float64)
    raw_totals = np.asarray(X_bg.sum(axis=1)).ravel()
    totals = raw_totals.copy(); totals[totals == 0] = 1.0
    row_idx = np.repeat(np.arange(X_bg.shape[0]), np.diff(X_bg.indptr))
    X_bg.data *= (1e6 / totals[row_idx])

    all_genes = features_df["feature_name"].astype(str).tolist()
    g2i = {g: i for i, g in enumerate(all_genes)}
    found = [g for g in genes if g in g2i]
    missing = [g for g in genes if g not in g2i]
    if missing:
        print(f"  genes not in dataset: {missing}")

    sel = [g2i[g] for g in found]
    X_sel = X_bg[:, sel].toarray()          # bins x n_found

    rr = rings[valid]
    X_sel = X_sel[valid]
    keep = raw_totals[valid] >= min_bin_counts     # bins with real counts
    unique_rings = sorted(np.unique(rr).astype(int).tolist())
    if exclude_last_ring and len(unique_rings) > 1:
        unique_rings = unique_rings[:-1]

    means, stds, ns, ns_all = {}, {}, {}, {}
    for gi, gene in enumerate(found):
        expr = X_sel[:, gi]
        m, s, n, na = [], [], [], []
        for r in unique_rings:
            in_r = (rr == r)
            counted = in_r & keep                  # average only over counted bins
            vals = expr[counted]
            m.append(float(vals.mean()) if vals.size else 0.0)
            s.append(float(vals.std(ddof=1)) if vals.size > 1 else 0.0)
            n.append(int(counted.sum()))
            na.append(int(in_r.sum()))
        means[gene] = np.array(m)
        stds[gene] = np.array(s)
        ns[gene] = np.array(n)
        ns_all[gene] = np.array(na)

    return unique_rings, means, stds, ns, found, ns_all

def smooth1d(y, sigma):
    """Gaussian-smooth a 1-D profile along the ring axis for display. sigma in rings;
    0/None = no smoothing. Edge-safe (mode='nearest')."""
    y = np.asarray(y, float)
    if not sigma or sigma <= 0 or y.size < 3:
        return y
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(y, sigma=float(sigma), mode="nearest")

def profiles_to_dataframe(unique_rings, means, stds, ns, found, ring_width_um, ns_all=None):
    """Long-format tidy DataFrame for saving to CSV.

    n_bins       = bins with counts averaged into the ring (denominator of the mean)
    n_bins_all   = all bins assigned to the ring
    frac_empty   = 1 - n_bins/n_bins_all  (high near the luminal border => the drawn
                   border sits off the tissue; those rings' means use only counted bins)
    """
    rows = []
    n_rings = len(unique_rings)
    for gene in found:
        # rings are ordered lumen(1) -> basal(n); give distance from BOTH borders
        for i, r in enumerate(unique_rings):
            n_ct = int(ns[gene][i])
            n_all = int(ns_all[gene][i]) if ns_all is not None else n_ct
            rows.append({
                "gene": gene,
                "ring": r,
                "ring_order_lumen_to_basal": i + 1,
                "dist_from_lumen_um": (i + 0.5) * ring_width_um,
                "dist_from_basal_um": (n_rings - i - 0.5) * ring_width_um,
                "mean_cpm": means[gene][i],
                "std_cpm": stds[gene][i],
                "n_bins": n_ct,
                "n_bins_all": n_all,
                "frac_empty": round(1 - n_ct / n_all, 3) if n_all else 0.0,
            })
    return pd.DataFrame(rows)
