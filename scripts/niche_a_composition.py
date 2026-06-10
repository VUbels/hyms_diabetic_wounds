"""
###############################################################################
SPATIAL XENIUM 5K PIPELINE — niche_a  (revised)
Niche composition + condition
###############################################################################
What changed vs the previous version
###############################################################################
  * condition score (condition signature, sc.tl.score_genes)
    computed once and cached so all three Step-2.5 scripts share it.
  * Exact Mann-Whitney (not the normal-approximation `ranksums`) for the
    sample-level proportion tests, which is the only honest test at n=3 vs n=2.
  * Linear mixed models (statsmodels) with the SAMPLE as a random intercept:
        - per cell type   : prop ~ C(condition) + C(niche) + (1|sample)
        - per niche (cond) : score ~ C(condition) + (1|sample)   [cell-level]
        - per cell type   : score ~ C(condition) + (1|sample)   [cell-level]
    The cell-level condition mixed models are properly powered (many cells per
    sample) and partition out sample variance — this is where defensible
    significance comes from given the limited number of biological replicates.
  * Shared, similarity-ordered niche colour map (compositionally similar niches
    get adjacent hues) and reuse of the R pipeline's cell-type colour map, so
    the SAME niche/cell type is the SAME colour in every figure and in both
    conditions.
  * Significance markers (*, **, ***) drawn on every comparison where a test
    supports them; captions state the test and the n explicitly.

Inputs
###############################################################################
    nichecompass_results/objects/nichecompass_integrated.h5ad   (override with --h5ad)

Key outputs (niche_composition/)
###############################################################################
    condition_by_celltype.pdf / .csv
    condition_by_niche.pdf / .csv
    condition_vs_composition_shift.pdf / .csv
    condition_gp_association.pdf / .csv
    niche_composition_combined.pdf
    differential.pdf  (+ niche_composition_differential.csv, mixedlm_*.csv)
    niche_proportions_per_sample.pdf
###############################################################################
"""

import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

import common_py_functions as cf
from common_py_functions import *  # shared keys, labels, colours, helpers

import argparse as _argparse
_parser = _argparse.ArgumentParser(description=__doc__)
cf.add_common_args(_parser)
_ARGS = _parser.parse_args()
cf.apply_common_args(_ARGS)

###############################################################################
# CONFIG
###############################################################################
INPUT_H5AD = _ARGS.h5ad or cf.DEFAULT_H5AD
OUTPUT_DIR = _ARGS.output_dir or "niche_composition"


# Optionally drop a differently-sourced slide here, e.g. ["Sample_3"].
EXCLUDE_SAMPLES = [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]

NEIGHBOR_K = 15          # local cell-type abundance neighbourhood
MIN_CELLS_NICHE = 10     # ignore a niche in a subset below this many cells

os.makedirs(OUTPUT_DIR, exist_ok=True)

###############################################################################
# 1. LOAD
###############################################################################
print(f"[LOAD] {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

for k in [SPATIAL_KEY]:
    if k not in adata.obsm:
        raise KeyError(f"Missing obsm['{k}']")
for k in [CELLTYPE_KEY, NICHE_KEY, SAMPLE_KEY, CONDITION_KEY]:
    if k not in adata.obs.columns:
        raise KeyError(f"Missing obs['{k}']")

if EXCLUDE_SAMPLES:
    keep = ~adata.obs[SAMPLE_KEY].astype(str).isin(EXCLUDE_SAMPLES)
    print(f"  Excluding samples {EXCLUDE_SAMPLES}: dropping {(~keep).sum()} cells")
    adata = adata[keep].copy()

adata.obs[NICHE_KEY] = adata.obs[NICHE_KEY].astype(str)
adata.obs[CELLTYPE_KEY] = adata.obs[CELLTYPE_KEY].astype(str)
adata.obs[CONDITION_KEY] = adata.obs[CONDITION_KEY].astype(str)
adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype(str)

all_niches = sorted(adata.obs[NICHE_KEY].unique(), key=niche_sort_key)
all_cts = sorted(adata.obs[CELLTYPE_KEY].unique())
all_samples = sorted(adata.obs[SAMPLE_KEY].unique())
sample_cond_map = adata.obs.groupby(SAMPLE_KEY)[CONDITION_KEY].first().to_dict()
n_test = sum(1 for v in sample_cond_map.values() if v == "TEST")
n_cntrl = sum(1 for v in sample_cond_map.values() if v == "CNTRL")

print(f"  {adata.n_obs} cells | {len(all_niches)} niches | {len(all_cts)} cell types")
print(f"  Samples: {len(all_samples)}  (TEST={n_test}, CNTRL={n_cntrl})")
N_NOTE = f"n={n_test} TEST vs n={n_cntrl} CNTRL slides"

###############################################################################
# 2. CONDITION SCORE + STYLES
###############################################################################
print("\n[COND] Scoring condition signature...")
cond_genes_used = compute_condition_score(adata)
score_thr = add_condition_flag(adata)
HAS_COND = CONDITION_SCORE_KEY in adata.obs and not adata.obs[CONDITION_SCORE_KEY].isna().all()
if HAS_COND:
    print(f"  Condition-high-cell threshold (q{int(CONDITION_SCORE_QUANTILE*100)}) = {score_thr:.4f}; "
          f"{int(adata.obs[CONDITION_FLAG_KEY].sum())} condition-high cells "
          f"({100*adata.obs[CONDITION_FLAG_KEY].mean():.1f}%)")

niche_colors, niche_order = build_or_load_niche_style(adata, NICHE_KEY, CELLTYPE_KEY)
ct_colors = build_or_load_celltype_colors(adata, CELLTYPE_KEY)
# niche_order is the global similarity order; restrict to those present here
niche_order = [n for n in niche_order if n in all_niches]

# GP activity matrix (signed) for condition association
gp_df = None
if HAS_COND and LATENT_KEY in adata.obsm:
    latent = np.asarray(adata.obsm[LATENT_KEY])
    gp_names = list(adata.uns.get(ACTIVE_GP_KEY, []))
    if len(gp_names) == latent.shape[1]:
        gp_df = pd.DataFrame(latent, index=adata.obs_names, columns=gp_names)
    else:
        gp_df = pd.DataFrame(latent, index=adata.obs_names,
                             columns=[f"latent_{i}" for i in range(latent.shape[1])])
        print(f"  [GP][WARN] active GP names ({len(gp_names)}) != latent dim "
              f"({latent.shape[1]}); using generic latent labels")

###############################################################################
# 3. LOCAL CELL TYPE ABUNDANCE (per cell, self EXCLUDED)
###############################################################################
print(f"\n[NEIGHBORS] Local cell-type abundance (k={NEIGHBOR_K}, self-excluded)...")
abundance = np.zeros((adata.n_obs, len(all_cts)))
ct_to_idx = {ct: i for i, ct in enumerate(all_cts)}
for slide in all_samples:
    mask = (adata.obs[SAMPLE_KEY] == slide).values
    idx = np.where(mask)[0]
    coords = adata.obsm[SPATIAL_KEY][mask]
    ct_labels = adata.obs[CELLTYPE_KEY].values[mask]
    if len(coords) < 2:
        continue
    tree = cKDTree(coords)
    k = min(NEIGHBOR_K + 1, len(coords))           # +1 for self
    _, nbrs = tree.query(coords, k=k)
    if nbrs.ndim == 1:
        nbrs = nbrs[:, None]
    nbrs = nbrs[:, 1:]                              # drop self
    for i, neigh in enumerate(nbrs):
        for ct in ct_labels[neigh]:
            abundance[idx[i], ct_to_idx[ct]] += 1
    denom = max(nbrs.shape[1], 1)
    abundance[idx, :] /= denom
abundance_df = pd.DataFrame(abundance, columns=all_cts, index=adata.obs_names)

###############################################################################
# 4. PER-CONDITION NICHE x CELLTYPE PROPORTION + CORRELATION
###############################################################################
def compute_niche_ct_panel(adata_sub, abundance_sub, niches_all, cts_all):
    niches = adata_sub.obs[NICHE_KEY].values
    cts = adata_sub.obs[CELLTYPE_KEY].values
    prop = pd.DataFrame(0.0, index=niches_all, columns=cts_all)
    corr = pd.DataFrame(np.nan, index=niches_all, columns=cts_all)
    for niche in niches_all:
        nm = niches == niche
        if nm.sum() < MIN_CELLS_NICHE:
            continue
        for ct in cts_all:
            prop.at[niche, ct] = (cts[nm] == ct).mean()
        ind = nm.astype(float)
        if ind.std() == 0:
            continue
        for ct in cts_all:
            v = abundance_sub[ct].values
            if v.std() > 0:
                corr.at[niche, ct] = pearsonr(ind, v)[0]
    return prop, corr


per_condition = {}
for cond in CONDITION_ORDER:
    m = (adata.obs[CONDITION_KEY] == cond).values
    if m.sum() == 0:
        continue
    sub = adata[m]
    prop, corr = compute_niche_ct_panel(sub, abundance_df.loc[sub.obs_names],
                                        all_niches, all_cts)
    per_condition[cond] = {"prop": prop, "corr": corr}
    prop.to_csv(os.path.join(OUTPUT_DIR, f"niche_composition_prop_{cond}.csv"))
    corr.to_csv(os.path.join(OUTPUT_DIR, f"niche_composition_corr_{cond}.csv"))

###############################################################################
# 5. NICHE ENRICHMENT BY CONDITION (sample-level, EXACT Mann-Whitney)
###############################################################################
print("\n[ENRICH] Niche abundance TEST vs CNTRL (exact Mann-Whitney, sample-level)...")
niche_by_sample = pd.crosstab(adata.obs[NICHE_KEY], adata.obs[SAMPLE_KEY],
                              normalize="columns").reindex(all_niches)
niche_by_sample.to_csv(os.path.join(OUTPUT_DIR, "niche_sample_frequency.csv"))

enr_rows = []
for niche in all_niches:
    props = niche_by_sample.loc[niche]
    test_props = np.array([props[s] for s in props.index if sample_cond_map[s] == "TEST"])
    cntrl_props = np.array([props[s] for s in props.index if sample_cond_map[s] == "CNTRL"])
    m_test, m_cntrl = (test_props.mean() if len(test_props) else 0.0), (cntrl_props.mean() if len(cntrl_props) else 0.0)
    p = np.nan
    if len(test_props) >= 1 and len(cntrl_props) >= 1 and (test_props.std() > 0 or cntrl_props.std() > 0) \
            and len(test_props) + len(cntrl_props) >= 4:
        try:
            p = mannwhitneyu(test_props, cntrl_props, alternative="two-sided", method="exact")[1]
        except Exception:
            p = np.nan
    enr_rows.append({"niche": niche, "mean_prop_TEST": m_test, "mean_prop_CNTRL": m_cntrl,
                     "sd_TEST": test_props.std() if len(test_props) > 1 else np.nan,
                     "sd_CNTRL": cntrl_props.std() if len(cntrl_props) > 1 else np.nan,
                     "log2FC_TEST_vs_CNTRL": np.log2((m_test + 1e-6) / (m_cntrl + 1e-6)),
                     "mw_p": p, "n_TEST": len(test_props), "n_CNTRL": len(cntrl_props)})
enrichment_df = pd.DataFrame(enr_rows)
vp = enrichment_df["mw_p"].dropna()
enrichment_df["mw_padj"] = np.nan
if len(vp) > 0:
    enrichment_df.loc[vp.index, "mw_padj"] = multipletests(vp, method="fdr_bh")[1]
enrichment_df.to_csv(os.path.join(OUTPUT_DIR, "niche_condition_enrichment.csv"), index=False)

###############################################################################
# 6. DIFFERENTIAL NICHE x CELLTYPE COMPOSITION
#    - sample-level exact Mann-Whitney (transparent, honest floor)
#    - per-cell-type linear mixed model (sample random intercept): the
#      principled, pooled test of a global compositional shift per cell type.
###############################################################################
print("\n[DIFF] Niche x cell-type composition: exact MW + mixed model...")

# per (sample, niche) within-niche cell-type proportions
rows = []
for sample in all_samples:
    sdat = adata.obs.loc[adata.obs[SAMPLE_KEY] == sample, [NICHE_KEY, CELLTYPE_KEY]]
    for niche in all_niches:
        nd = sdat[sdat[NICHE_KEY] == niche]
        n_in = len(nd)
        for ct in all_cts:
            rows.append({"sample": sample, "condition": sample_cond_map[sample],
                         "niche": niche, "celltype": ct,
                         "prop": (nd[CELLTYPE_KEY] == ct).sum() / max(n_in, 1),
                         "n_in_niche": n_in})
sn_ct = pd.DataFrame(rows)
sn_ct.to_csv(os.path.join(OUTPUT_DIR, "niche_celltype_proportions_per_sample.csv"),
             index=False)

# (a) per niche x celltype exact MW
diff_rows = []
for niche in all_niches:
    for ct in all_cts:
        d = sn_ct[(sn_ct["niche"] == niche) & (sn_ct["celltype"] == ct)]
        test_props = d[d["condition"] == "TEST"]["prop"].values
        cntrl_props = d[d["condition"] == "CNTRL"]["prop"].values
        m_test, m_cntrl = test_props.mean() if len(test_props) else 0, cntrl_props.mean() if len(cntrl_props) else 0
        p = np.nan
        if len(test_props) + len(cntrl_props) >= 4 and (np.std(test_props) > 0 or np.std(cntrl_props) > 0):
            try:
                p = mannwhitneyu(test_props, cntrl_props, alternative="two-sided", method="exact")[1]
            except Exception:
                p = np.nan
        diff_rows.append({"niche": niche, "celltype": ct,
                          "mean_prop_TEST": m_test, "mean_prop_CNTRL": m_cntrl,
                          "log2FC_TEST_vs_CNTRL": np.log2((m_test + 1e-6) / (m_cntrl + 1e-6)),
                          "mw_p": p})
diff_df = pd.DataFrame(diff_rows)
diff_df["mw_padj"] = np.nan
for niche in all_niches:                                   # FDR within niche
    m = (diff_df["niche"] == niche) & diff_df["mw_p"].notna()
    if m.sum() > 0:
        diff_df.loc[m, "mw_padj"] = multipletests(diff_df.loc[m, "mw_p"],
                                                  method="fdr_bh")[1]

# (b) per cell-type mixed model: prop ~ C(condition) + C(niche) + (1|sample)
ct_mm_rows = []
for ct in all_cts:
    d = sn_ct[sn_ct["celltype"] == ct]
    res = mixedlm_condition(d, "prop", sample_col="sample", extra_fixed=["niche"])
    ct_mm_rows.append({"celltype": ct, **res})
ct_mm = pd.DataFrame(ct_mm_rows)
vp = ct_mm["p"].dropna()
ct_mm["padj"] = np.nan
if len(vp) > 0:
    ct_mm.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
ct_mm.to_csv(os.path.join(OUTPUT_DIR, "mixedlm_celltype_global_shift.csv"), index=False)
diff_df = diff_df.merge(ct_mm[["celltype", "p", "padj"]].rename(
    columns={"p": "celltype_mm_p", "padj": "celltype_mm_padj"}), on="celltype", how="left")
diff_df.to_csv(os.path.join(OUTPUT_DIR, "niche_composition_differential.csv"), index=False)
ct_sig = {r.celltype: r.padj for r in ct_mm.itertuples()}

###############################################################################
# 7. CONDITION: per cell type and per niche (cell-level mixed models)
###############################################################################
cond_ct = cond_niche = None
gp_cond = None
if HAS_COND:
    print("\n[COND] Condition by cell type and niche (cell-level mixed models)...")
    cell_df = adata.obs[[CONDITION_SCORE_KEY, CELLTYPE_KEY, NICHE_KEY,
                         SAMPLE_KEY, CONDITION_KEY, CONDITION_FLAG_KEY]].copy()
    cell_df.columns = ["score", "celltype", "niche", "sample", "condition", "condition-high"]

    # per cell type
    rows = []
    for ct in all_cts:
        d = cell_df[cell_df["celltype"] == ct]
        res = mixedlm_condition(d, "score", sample_col="sample")
        rows.append({"celltype": ct, "mean_score": d["score"].mean(),
                     "mean_score_TEST": d[d.condition == "TEST"]["score"].mean(),
                     "mean_score_CNTRL": d[d.condition == "CNTRL"]["score"].mean(),
                     "frac_condition_high": d["condition-high"].mean(), "n_cells": len(d), **res})
    cond_ct = pd.DataFrame(rows)
    vp = cond_ct["p"].dropna()
    cond_ct["padj"] = np.nan
    if len(vp) > 0:
        cond_ct.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
    cond_ct = cond_ct.sort_values("mean_score", ascending=False)
    cond_ct.to_csv(os.path.join(OUTPUT_DIR, "condition_by_celltype.csv"), index=False)

    # per niche
    rows = []
    for niche in all_niches:
        d = cell_df[cell_df["niche"] == niche]
        res = mixedlm_condition(d, "score", sample_col="sample")
        rows.append({"niche": niche, "mean_score": d["score"].mean(),
                     "mean_score_TEST": d[d.condition == "TEST"]["score"].mean(),
                     "mean_score_CNTRL": d[d.condition == "CNTRL"]["score"].mean(),
                     "frac_condition_high": d["condition-high"].mean(), "n_cells": len(d), **res})
    cond_niche = pd.DataFrame(rows)
    vp = cond_niche["p"].dropna()
    cond_niche["padj"] = np.nan
    if len(vp) > 0:
        cond_niche.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
    cond_niche.to_csv(os.path.join(OUTPUT_DIR, "condition_by_niche.csv"), index=False)

    # GP association (per-cell Spearman; exploratory ranking)
    if gp_df is not None:
        s = adata.obs[CONDITION_SCORE_KEY].values
        rows = []
        for gp in gp_df.columns:
            rho, pv = spearmanr(gp_df[gp].values, s)
            rows.append({"gene_program": gp, "spearman_rho": rho, "p": pv})
        gp_cond = pd.DataFrame(rows).dropna(subset=["spearman_rho"])
        if len(gp_cond):
            gp_cond["padj"] = multipletests(gp_cond["p"].fillna(1.0),
                                           method="fdr_bh")[1]
            gp_cond = gp_cond.reindex(gp_cond["spearman_rho"].abs()
                                    .sort_values(ascending=False).index)
        gp_cond.to_csv(os.path.join(OUTPUT_DIR, "condition_gp_association.csv"),
                      index=False)

###############################################################################
# 8. PLOTS
###############################################################################
def _legend_conditions(fig, loc="upper right"):
    handles = [plt.Line2D([0], [0], marker="o", color="w", markersize=8,
                          markerfacecolor=COND_COLORS[c], label=COND_SHORT[c])
               for c in CONDITION_ORDER]
    fig.legend(handles=handles, loc=loc, frameon=False)


###############################################################################
# 8.1 combined TEST/CNTRL composition dot plot (shared niche colours, sig col)
###############################################################################
def plot_combined_composition(out_path, figsize=(11, 7)):
    conds = [c for c in CONDITION_ORDER if c in per_condition]
    if len(conds) < 1:
        return
    niches_info = [n for n in niche_order
                   if any(per_condition[c]["prop"].loc[n].sum() > 0 for c in conds)]
    if not niches_info:
        return
    prop_ref = per_condition[conds[0]]["prop"].reindex(niches_info)
    ct_order = prop_ref.max(axis=0).sort_values(ascending=False).index.tolist()
    n_n, n_c = len(niches_info), len(ct_order)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, len(conds) + 1, width_ratios=[n_c] * len(conds) + [1.4],
                          height_ratios=[20, 1], hspace=0.35, wspace=0.18)
    axes = []
    for i in range(len(conds) + 1):
        axes.append(fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None))

    for ax, cond in zip(axes[:-1], conds):
        prop = per_condition[cond]["prop"].reindex(niches_info)[ct_order]
        corr = per_condition[cond]["corr"].reindex(niches_info)[ct_order]
        for yi, niche in enumerate(niches_info):
            for xi, ct in enumerate(ct_order):
                p = prop.at[niche, ct]
                if p > 0:
                    r = corr.at[niche, ct]
                    ax.scatter(xi, yi, s=p * 420,
                               c=[(0 if np.isnan(r) else r)], cmap="RdBu_r",
                               vmin=-0.5, vmax=0.5, edgecolors="black",
                               linewidths=0.3, alpha=0.9)
        ax.set_xticks(range(n_c))
        ax.set_xticklabels([f"{ct}{(' '+pval_stars(ct_sig.get(ct, np.nan))).rstrip()}"
                            for ct in ct_order], rotation=45, ha="right")
        ax.set_yticks(range(n_n))
        ax.set_yticklabels([f"N{n}" for n in niches_info])
        color_yticklabels(ax, niches_info, niche_colors)
        ax.set_xlim(-0.5, n_c - 0.5)
        ax.set_ylim(-0.5, n_n - 0.5)
        ax.set_title(COND_SHORT[cond], fontsize=11)
        ax.invert_yaxis()
        ax.grid(alpha=0.2, ls="--", lw=0.3)
    for ax in axes[1:-1]:
        plt.setp(ax.get_yticklabels(), visible=False)

    # differential column (niche enrichment TEST/CNTRL, stars from exact MW)
    axd = axes[-1]
    for yi, niche in enumerate(niches_info):
        row = enrichment_df[enrichment_df["niche"] == niche]
        if len(row) == 0:
            continue
        lfc = row["log2FC_TEST_vs_CNTRL"].iloc[0]
        padj = row["mw_padj"].iloc[0]
        if np.isfinite(lfc):
            sig = (not np.isnan(padj)) and padj < 0.05
            axd.scatter(0, yi, s=min(abs(lfc) * 90, 320) if sig else max(abs(lfc) * 45, 18),
                        c=[lfc], cmap="RdBu_r", vmin=-2, vmax=2,
                        edgecolors="black" if sig else "gray",
                        linewidths=1.2 if sig else 0.3)
            star = pval_stars(padj)
            if star:
                axd.text(0.32, yi, star, va="center", fontsize=9)
    axd.set_xticks([0]); axd.set_xticklabels(["TEST/CNTRL"])
    axd.set_xlim(-0.5, 0.6); axd.set_ylim(-0.5, n_n - 0.5)
    axd.invert_yaxis(); axd.grid(alpha=0.2, ls="--", lw=0.3)
    axd.set_title("Niche\nabundance", fontsize=10)
    plt.setp(axd.get_yticklabels(), visible=False)

    cax = fig.add_subplot(gs[1, :len(conds)])
    sm = mpl.cm.ScalarMappable(cmap="RdBu_r",
                               norm=mpl.colors.Normalize(vmin=-0.5, vmax=0.5))
    fig.colorbar(sm, cax=cax, orientation="horizontal",
                 label="Pearson r (local cell-type abundance vs niche membership)")
    fig.suptitle(f"Niche composition by condition  ·  dot size = within-niche "
                 f"proportion  ·  {N_NOTE}", fontsize=12, y=1.02)
    plt.savefig(out_path)
    plt.close()


print("\n[PLOT] composition dot plot...")
plot_combined_composition(os.path.join(OUTPUT_DIR, "niche_composition_combined.pdf"))

###############################################################################
# 8.2 differential composition heatmap (log2FC, celltype-level sig)
###############################################################################
print("[PLOT] differential composition heatmap...")
lfc_piv = diff_df.pivot(index="niche", columns="celltype",
                        values="log2FC_TEST_vs_CNTRL").reindex(niche_order)[all_cts]
padj_piv = diff_df.pivot(index="niche", columns="celltype",
                         values="mw_padj").reindex(niche_order)[all_cts]
annot = padj_piv.apply(lambda col: col.map(lambda p: pval_stars(p)))
fig, ax = plt.subplots(figsize=(max(5, len(all_cts) * 0.34),
                                max(3, len(niche_order) * 0.30)))
sns.heatmap(lfc_piv, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            annot=annot.values, fmt="", annot_kws={"size": 7},
            linewidths=0.3, ax=ax,
            cbar_kws={"label": "log2(TEST/CNTRL) within-niche proportion", "shrink": 0.6})
ax.set_yticklabels([f"N{n}" for n in niche_order], rotation=0)
color_yticklabels(ax, niche_order, niche_colors)
ax.set_xticklabels([f"{ct}{(' '+pval_stars(ct_sig.get(ct, np.nan))).rstrip()}"
                    for ct in all_cts], rotation=45, ha="right")
ax.set_title(f"Differential niche composition: TEST vs CNTRL\n"
             f"cell stars = exact MW (FDR/niche); x-label stars = mixed-model "
             f"global shift (FDR)\n{N_NOTE} — interpret as exploratory", fontsize=8)
ax.set_xlabel("Cell type"); ax.set_ylabel("Niche")
plt.savefig(os.path.join(OUTPUT_DIR, "differential.pdf"))
plt.close()

###############################################################################
# 8.3 per-sample niche proportions (each dot = sample)
###############################################################################
print("[PLOT] per-sample niche proportions...")
ncol = min(4, len(all_niches)); nrow = int(np.ceil(len(all_niches) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2.5 * nrow), squeeze=False)
long_ns = niche_by_sample.reset_index().melt(id_vars=NICHE_KEY, var_name="sample",
                                             value_name="prop")
long_ns["condition"] = long_ns["sample"].map(sample_cond_map)
rng = np.random.default_rng(0)
for i, niche in enumerate(niche_order):
    ax = axes[i // ncol][i % ncol]
    nd = long_ns[long_ns[NICHE_KEY] == niche]
    for ci, cond in enumerate(CONDITION_ORDER):
        cd = nd[nd["condition"] == cond]
        ax.scatter(rng.uniform(-0.13, 0.13, len(cd)) + ci, cd["prop"],
                   c=COND_COLORS[cond], s=42, edgecolors="black", linewidths=0.5, zorder=3)
        if len(cd):
            ax.plot([ci - 0.2, ci + 0.2], [cd["prop"].mean()] * 2,
                    color=COND_COLORS[cond], lw=2, zorder=2)
    row = enrichment_df[enrichment_df["niche"] == niche]
    p = row["mw_p"].iloc[0] if len(row) else np.nan
    star = pval_stars(p)
    ttl = f"N{niche}" + (f"  {star}" if star else (f"  p={p:.2f}" if np.isfinite(p) else ""))
    ax.set_title(ttl, fontsize=8, color=niche_colors.get(niche, "black"))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["TEST", "CNTRL"])
    ax.set_ylabel("Proportion", fontsize=7)
    ax.margins(y=0.15)
for j in range(len(all_niches), nrow * ncol):
    axes[j // ncol][j % ncol].axis("off")
fig.suptitle(f"Per-sample niche proportions (each dot = one slide)  ·  "
             f"exact MW  ·  {N_NOTE}", fontsize=10, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_proportions_per_sample.pdf"))
plt.close()

###############################################################################
# 9. CONDITION FIGURES
###############################################################################
def split_box_strip(ax, df, group_col, group_order, color_map, sig_map=None,
                    value_col="score"):
    """Per-group TEST/CNTRL split: box + jittered per-sample-aware points + stars."""
    rng = np.random.default_rng(0)
    for gi, g in enumerate(group_order):
        gd = df[df[group_col] == g]
        for ci, cond in enumerate(CONDITION_ORDER):
            cd = gd[gd["condition"] == cond][value_col].values
            if len(cd) == 0:
                continue
            xc = gi + (ci - 0.5) * 0.34
            bp = ax.boxplot([cd], positions=[xc], widths=0.28, showfliers=False,
                            patch_artist=True, manage_ticks=False)
            for box in bp["boxes"]:
                box.set(facecolor=COND_COLORS[cond], alpha=0.55, edgecolor="black",
                        linewidth=0.6)
            for med in bp["medians"]:
                med.set(color="black", linewidth=1.0)
        # significance bracket from sig_map (mixed-model padj)
        if sig_map is not None:
            star = pval_stars(sig_map.get(g, np.nan))
            if star:
                ytop = gd[value_col].quantile(0.97)
                ax.plot([gi - 0.17, gi + 0.17], [ytop, ytop], color="black", lw=0.8)
                ax.text(gi, ytop, star, ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(group_order)))


if HAS_COND:
    # 9.1 condition by cell type
    print("[PLOT] condition by cell type...")
    order_ct = cond_ct["celltype"].tolist()
    cdf = adata.obs[[CONDITION_SCORE_KEY, CELLTYPE_KEY, CONDITION_KEY]].copy()
    cdf.columns = ["score", "celltype", "condition"]
    fig, ax = plt.subplots(figsize=(max(7, len(order_ct) * 0.7), 4.2))
    split_box_strip(ax, cdf, "celltype", order_ct, ct_colors,
                    sig_map={r.celltype: r.padj for r in cond_ct.itertuples()})
    ax.set_xticklabels(order_ct, rotation=45, ha="right")
    # colour the x tick labels by cell-type colour
    for tick, ct in zip(ax.get_xticklabels(), order_ct):
        tick.set_color(ct_colors.get(ct, "black"))
    ax.set_ylabel("condition score")
    ax.set_title(f"Condition by cell type (TEST vs CNTRL)\nstars = mixed-model "
                 f"condition effect, FDR; {N_NOTE}", fontsize=9)
    _legend_conditions(fig)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "condition_by_celltype.pdf"))
    plt.close()

    # 9.2 condition by niche
    print("[PLOT] condition by niche...")
    order_ni = cond_niche.sort_values("mean_score", ascending=False)["niche"].tolist()
    cdf = adata.obs[[CONDITION_SCORE_KEY, NICHE_KEY, CONDITION_KEY]].copy()
    cdf.columns = ["score", "niche", "condition"]
    fig, ax = plt.subplots(figsize=(max(7, len(order_ni) * 0.7), 4.2))
    split_box_strip(ax, cdf, "niche", order_ni, niche_colors,
                    sig_map={r.niche: r.padj for r in cond_niche.itertuples()})
    ax.set_xticklabels([f"N{n}" for n in order_ni])
    for tick, n in zip(ax.get_xticklabels(), order_ni):
        tick.set_color(niche_colors.get(n, "black"))
    ax.set_ylabel("condition score")
    ax.set_title(f"Condition by niche (TEST vs CNTRL)\nstars = mixed-model condition "
                 f"effect, FDR; {N_NOTE}", fontsize=9)
    _legend_conditions(fig)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "condition_by_niche.pdf"))
    plt.close()

    # 9.3 is the niche distribution shift associated with condition?
    print("[PLOT] condition vs compositional shift...")
    merged = cond_niche.merge(enrichment_df[["niche", "log2FC_TEST_vs_CNTRL", "mw_padj"]],
                             on="niche", how="left")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    # (left) niche condition level vs niche abundance log2FC
    ax = axes[0]
    for r in merged.itertuples():
        ax.scatter(r.mean_score, r.log2FC_TEST_vs_CNTRL, s=90,
                   color=niche_colors.get(r.niche, "gray"),
                   edgecolors="black", linewidths=0.5, zorder=3)
        ax.annotate(f"N{r.niche}", (r.mean_score, r.log2FC_TEST_vs_CNTRL),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, ls="--", color="gray", lw=0.6)
    valid = merged.dropna(subset=["mean_score", "log2FC_TEST_vs_CNTRL"])
    if len(valid) >= 3:
        rho, pv = spearmanr(valid["mean_score"], valid["log2FC_TEST_vs_CNTRL"])
        ax.set_title(f"Niche condition vs abundance shift\nSpearman ρ={rho:.2f}, "
                     f"p={pv:.2f}", fontsize=9)
    ax.set_xlabel("Mean niche condition score")
    ax.set_ylabel("log2(TEST/CNTRL) niche abundance")
    # (right) per-niche TEST-CNTRL condition effect vs abundance log2FC
    ax = axes[1]
    for r in merged.itertuples():
        eff = getattr(r, "effect_CNTRL_vs_TEST", np.nan)
        if np.isnan(eff):
            continue
        sig = (not np.isnan(getattr(r, "padj", np.nan))) and r.padj < 0.05
        ax.scatter(eff, r.log2FC_TEST_vs_CNTRL, s=110 if sig else 70,
                   color=niche_colors.get(r.niche, "gray"),
                   edgecolors="black" if sig else "gray",
                   linewidths=1.3 if sig else 0.5, zorder=3)
        ax.annotate(f"N{r.niche}" + ("*" if sig else ""),
                    (eff, r.log2FC_TEST_vs_CNTRL), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, ls="--", color="gray", lw=0.6)
    ax.axvline(0, ls="--", color="gray", lw=0.6)
    ax.set_xlabel("Condition TEST−CNTRL effect (mixed model)")
    ax.set_ylabel("log2(TEST/CNTRL) niche abundance")
    ax.set_title("Is the abundance shift coupled to a\ncondition shift? "
                 "(* niche condition FDR<0.05)", fontsize=9)
    fig.suptitle(f"Niche redistribution vs condition  ·  {N_NOTE}",
                 fontsize=11, y=1.03)
    merged.to_csv(os.path.join(OUTPUT_DIR, "condition_vs_composition_shift.csv"),
                  index=False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "condition_vs_composition_shift.pdf"))
    plt.close()

    # 9.4 condition-associated gene programs
    if gp_cond is not None and len(gp_cond) > 0:
        print("[PLOT] condition-associated gene programs...")
        top = pd.concat([gp_cond.head(15), gp_cond.tail(15)]).drop_duplicates("gene_program")
        top = top.sort_values("spearman_rho")
        fig, ax = plt.subplots(figsize=(7, max(4, len(top) * 0.28)))
        colors = ["#C94A4A" if v > 0 else "#4A7EB8" for v in top["spearman_rho"]]
        ax.barh(range(len(top)), top["spearman_rho"], color=colors,
                edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([gp[:48] for gp in top["gene_program"]], fontsize=6)
        for yi, padj in enumerate(top["padj"]):
            star = pval_stars(padj)
            if star:
                v = top["spearman_rho"].iloc[yi]
                ax.text(v + (0.01 if v >= 0 else -0.01), yi, star,
                        va="center", ha="left" if v >= 0 else "right", fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("Spearman ρ (GP activity vs condition score)")
        ax.set_title("Condition-associated gene programs\n(per-cell ρ; stars FDR; "
                     "exploratory)", fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "condition_gp_association.pdf"))
        plt.close()

###############################################################################
# 10. DONE
###############################################################################
print(f"\n[DONE] niche_a complete. Outputs in {OUTPUT_DIR}/")
print(f"  Statistics: niche abundance & within-niche composition use the SLIDE")
print(f"  as the unit (exact Mann-Whitney). Cell-type compositional shift and all")
print(f"  condition comparisons use linear mixed models with sample as a random")
print(f"  intercept. {N_NOTE}; treat per-pair proportion p-values as exploratory.")
if HAS_COND:
    print(f"  Condition: {len(cond_genes_used)} condition genes used "
          f"(see {SHARED_DIR}/condition_genes_used.txt).")