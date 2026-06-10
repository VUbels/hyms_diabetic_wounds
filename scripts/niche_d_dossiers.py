"""
###############################################################################
SPATIAL XENIUM 5K PIPELINE — niche_d
Condition niche dossier (focused, narrative)
###############################################################################
Purpose
###############################################################################
Instead of plotting every niche, this selects the handful of niches that carry
signal and lays each one out as a single row that reads left-to-right as the
biological chain:

    which niche  ·  is it TEST- or CNTRL-skewed, and how condition-high (TEST vs CNTRL)
        -> what cell types is it built from
            -> which gene programs track its condition

It is run TWICE, writing to two separate folders so each story is self-contained:
  * niche_dossier_condition/  — niches ranked by ABSOLUTE condition level
                                 (the most condition-high niches, regardless of which
                                 condition they sit in)
  * niche_dossier_abundance/   — niches ranked by ABUNDANCE SHIFT
                                 (|log2(TEST/CNTRL)| — the niches that most change
                                 their prevalence between conditions)

Honesty
###############################################################################
With n=3 TEST vs n=2 CNTRL slides almost nothing is individually significant. The
condition Δ (TEST−CNTRL) shows an FDR star only when the cell-level mixed model
clears it; otherwise it reads "ns". Abundance is one observation per slide, so
its exact-MW p is floored at ~0.20 and is shown for transparency, not as proof.
GP correlations are per-cell Spearman and are exploratory. Everything here is an
effect-size / ranking view, labelled as such.

Inputs
###############################################################################
    nichecompass_results/objects/nichecompass_integrated.h5ad   (override with --h5ad)
    shared caches in niche_analysis_shared/ (written by step 2.5a; rebuilt if
    absent so this can also run first)

Outputs (per folder)
###############################################################################
    niche_dossier.pdf          the featured-niche dossier
    niche_overview_ranked.pdf  all niches ranked by the selection metric
    niche_dossier_table.csv    per-niche numbers behind the figure
###############################################################################
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import rankdata, mannwhitneyu
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
OUTPUT_COND = os.environ.get("DOSSIER_COND_DIR", "niche_dossier_condition")
OUTPUT_ABUND = os.environ.get("DOSSIER_ABUND_DIR", "niche_dossier_abundance")


EXCLUDE_SAMPLES = [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]

N_FEATURE = int(os.environ.get("DOSSIER_N_FEATURE", "8"))
TOP_CELLTYPES = 5
TOP_GPS = 5
MIN_CELLS_NICHE = 30          # a niche must have >= this many cells to be featured
MIN_CELLS_GP = 20             # min cells to compute within-niche GP correlation
ABUND_PSEUDO = 1e-4           # pseudocount for abundance log2 fold change

###############################################################################
# 1. LOAD + CONDITION + SHARED STYLES
###############################################################################
print(f"[LOAD] {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
for k in (NICHE_KEY, CELLTYPE_KEY, CONDITION_KEY, SAMPLE_KEY):
    adata.obs[k] = adata.obs[k].astype(str)
if EXCLUDE_SAMPLES:
    keep = ~adata.obs[SAMPLE_KEY].isin(EXCLUDE_SAMPLES)
    print(f"  Excluding {EXCLUDE_SAMPLES}: dropping {(~keep).sum()} cells")
    adata = adata[keep].copy()

print("\n[COND] Scoring condition signature...")
compute_condition_score(adata)
add_condition_flag(adata)
niche_colors, _ = build_or_load_niche_style(adata, NICHE_KEY, CELLTYPE_KEY)
ct_colors = build_or_load_celltype_colors(adata, CELLTYPE_KEY)

all_niches = sorted(adata.obs[NICHE_KEY].unique(), key=niche_sort_key)
samp_cond = adata.obs.groupby(SAMPLE_KEY)[CONDITION_KEY].first().to_dict()
test_samples = [s for s, c in samp_cond.items() if c == "TEST"]
cntrl_samples = [s for s, c in samp_cond.items() if c == "CNTRL"]
n_test, n_cntrl = len(test_samples), len(cntrl_samples)
N_NOTE = f"n={n_test} TEST vs n={n_cntrl} CNTRL slides"

# latent / gene-program activities
gp_names = list(adata.uns.get(ACTIVE_GP_KEY, []))
L = np.asarray(adata.obsm[LATENT_KEY]) if LATENT_KEY in adata.obsm else None
if L is not None and (not gp_names or len(gp_names) != L.shape[1]):
    gp_names = [f"GP_{i}" for i in range(L.shape[1])]
HAS_GP = L is not None and len(gp_names) > 0


###############################################################################
# 2. PER-NICHE METRICS
###############################################################################
print("\n[METRICS] Per-niche condition / abundance / composition / GP...")
score = adata.obs[CONDITION_SCORE_KEY].values
niche_arr = adata.obs[NICHE_KEY].values
cond_arr = adata.obs[CONDITION_KEY].values

niche_by_sample = pd.crosstab(adata.obs[NICHE_KEY], adata.obs[SAMPLE_KEY]) \
    .reindex(all_niches, fill_value=0)
sample_tot = adata.obs[SAMPLE_KEY].value_counts()
prop_table = niche_by_sample.div(sample_tot, axis=1)        # niche x sample
ct_by_niche = pd.crosstab(adata.obs[NICHE_KEY], adata.obs[CELLTYPE_KEY]) \
    .reindex(all_niches, fill_value=0)
ct_prop = ct_by_niche.div(ct_by_niche.sum(axis=1).replace(0, np.nan), axis=0)

cell_df = pd.DataFrame({"score": score, "niche": niche_arr,
                        "sample": adata.obs[SAMPLE_KEY].values, "condition": cond_arr})

by_niche, rows = {}, []
for niche in all_niches:
    m = niche_arr == niche
    n_cells = int(m.sum())
    s_all = score[m]
    s_test = score[m & (cond_arr == "TEST")]
    s_cntrl = score[m & (cond_arr == "CNTRL")]
    mean_score = float(np.nanmean(s_all)) if n_cells else np.nan
    mean_test = float(np.nanmean(s_test)) if len(s_test) else np.nan
    mean_cntrl = float(np.nanmean(s_cntrl)) if len(s_cntrl) else np.nan

    res = mixedlm_condition(cell_df[cell_df["niche"] == niche], "score",
                            sample_col="sample")

    test_props = prop_table.loc[niche, test_samples].values.astype(float)
    cntrl_props = prop_table.loc[niche, cntrl_samples].values.astype(float)
    mean_p_test, mean_p_cntrl = float(test_props.mean()), float(cntrl_props.mean())
    log2fc = float(np.log2((mean_p_test + ABUND_PSEUDO) / (mean_p_cntrl + ABUND_PSEUDO)))
    mw_p = np.nan
    if len(test_props) >= 1 and len(cntrl_props) >= 1 and \
            (test_props.std() + cntrl_props.std()) > 0:
        try:
            mw_p = float(mannwhitneyu(test_props, cntrl_props,
                                      alternative="two-sided")[1])
        except Exception:
            pass

    comp = [(c, float(ct_prop.loc[niche, c])) for c in
            ct_prop.loc[niche].sort_values(ascending=False).index[:TOP_CELLTYPES]
            if ct_prop.loc[niche, c] > 0]

    gp = []
    if HAS_GP and n_cells >= MIN_CELLS_GP:
        rho = spearman_cols(L[m], s_all)
        order = np.argsort(-np.abs(np.nan_to_num(rho)))
        gp = [(gp_names[j], float(rho[j])) for j in order[:TOP_GPS]
              if np.isfinite(rho[j])]

    d = {"niche": niche, "n_cells": n_cells, "mean_score": mean_score,
         "mean_score_TEST": mean_test, "mean_score_CNTRL": mean_cntrl,
         "score_delta_obs": (mean_test - mean_cntrl) if np.isfinite(mean_test)
         and np.isfinite(mean_cntrl) else np.nan,
         "score_effect_mixed": res["effect_CNTRL_vs_TEST"], "score_p": res["p"],
         "mean_prop_TEST": mean_p_test, "mean_prop_CNTRL": mean_p_cntrl,
         "abund_log2fc": log2fc, "abund_mw_p": mw_p,
         "comp": comp, "gp": gp}
    by_niche[niche] = d
    rows.append({k: v for k, v in d.items() if k not in ("comp", "gp")})

df = pd.DataFrame(rows)
for pc, pad in [("score_p", "score_padj"), ("abund_mw_p", "abund_padj")]:
    v = df[pc].notna()
    df[pad] = np.nan
    if v.sum():
        df.loc[v, pad] = multipletests(df.loc[v, pc], method="fdr_bh")[1]
score_padj = dict(zip(df["niche"], df["score_padj"]))

###############################################################################
# 3. FIGURES
###############################################################################
def draw_dossier(featured, out_dir, mode_label):
    K = len(featured)
    if K == 0:
        print(f"  [WARN] no niches to feature for {mode_label}")
        return
    vmax = max([max(by_niche[n]["mean_score_TEST"] if np.isfinite(by_niche[n]["mean_score_TEST"]) else 0,
                    by_niche[n]["mean_score_CNTRL"] if np.isfinite(by_niche[n]["mean_score_CNTRL"]) else 0)
                for n in featured] + [1e-3]) * 1.18
    fig = plt.figure(figsize=(13, 1.3 * K + 1.3))
    gs = fig.add_gridspec(K, 3, width_ratios=[1.15, 1.25, 1.25], wspace=0.45,
                          hspace=0.7, top=0.9, bottom=0.05, left=0.04, right=0.985)
    for r, nid in enumerate(featured):
        d = by_niche[nid]
        # status column
        axs = fig.add_subplot(gs[r, 0]); axs.axis("off")
        axs.set_xlim(0, 1); axs.set_ylim(0, 1)
        axs.add_patch(mpatches.Rectangle((0.0, 0.80), 0.05, 0.17,
                      color=niche_colors.get(nid, "#999"), clip_on=False,
                      transform=axs.transAxes))
        axs.text(0.08, 0.885, f"N{nid}", fontsize=13, fontweight="bold",
                 va="center", transform=axs.transAxes)
        is_test = d["abund_log2fc"] > 0
        pcol = COND_COLORS["TEST"] if is_test else COND_COLORS["CNTRL"]
        axs.text(0.0, 0.60, f"{'TEST' if is_test else 'CNTRL'}-enriched · "
                 f"log2FC {d['abund_log2fc']:+.1f}", fontsize=9, color=pcol,
                 transform=axs.transAxes)
        y, x0, w = 0.33, 0.0, 0.60
        axs.plot([x0, x0 + w], [y, y], color="0.85", lw=1, transform=axs.transAxes)
        mh, ml = d["mean_score_TEST"], d["mean_score_CNTRL"]
        if np.isfinite(mh) and np.isfinite(ml):
            xh, xl = x0 + w * mh / vmax, x0 + w * ml / vmax
            axs.plot([xl, xh], [y, y], color="0.55", lw=2, transform=axs.transAxes)
        for val, c in [(ml, COND_COLORS["CNTRL"]),
                       (mh, COND_COLORS["TEST"])]:
            if np.isfinite(val):
                axs.scatter([x0 + w * val / vmax], [y], s=55, color=c, zorder=3,
                            transform=axs.transAxes)
        axs.text(x0, 0.16, "0", fontsize=8, color="0.5", transform=axs.transAxes)
        axs.text(x0 + w, 0.16, f"{vmax:.2f}", fontsize=8, color="0.5", ha="right",
                 transform=axs.transAxes)
        star = pval_stars(score_padj.get(nid, np.nan))
        dd = d["score_delta_obs"]
        axs.text(0.0, 0.0, f"Δ {dd:+.3f} · {star or 'ns'}" if np.isfinite(dd)
                 else "Δ n/a", fontsize=9, color="0.25", transform=axs.transAxes)
        # composition column
        axc = fig.add_subplot(gs[r, 1])
        comp = d["comp"][:TOP_CELLTYPES][::-1]
        names = [c for c, _ in comp]; vals = [p for _, p in comp]
        if vals:
            axc.barh(range(len(names)), vals, color=[ct_colors.get(c, "#bbb")
                     for c in names], edgecolor="white", linewidth=0.5)
            axc.set_yticks(range(len(names))); axc.set_yticklabels(names, fontsize=9)
            axc.set_xlim(0, max(vals) * 1.22)
            for i, v in enumerate(vals):
                axc.text(v, i, f" {v*100:.0f}%", va="center", fontsize=8, color="0.3")
        axc.tick_params(length=0)
        for sp_ in ("top", "right", "left"):
            axc.spines[sp_].set_visible(False)
        axc.set_xticks([])
        # gene-program column
        axg = fig.add_subplot(gs[r, 2])
        gps = d["gp"][:TOP_GPS][::-1]
        gn = [shorten_gp(g) for g, _ in gps]; gr = [v for _, v in gps]
        if gr:
            axg.barh(range(len(gn)), gr, edgecolor="white", linewidth=0.5,
                     color=[COND_COLORS["TEST"] if v > 0
                            else COND_COLORS["CNTRL"] for v in gr])
            axg.set_yticks(range(len(gn))); axg.set_yticklabels(gn, fontsize=8)
            mm = max([abs(v) for v in gr] + [0.05])
            axg.set_xlim(-mm * 1.3, mm * 1.3)
            for i, v in enumerate(gr):
                axg.text(v + (0.005 if v >= 0 else -0.005), i, f"{v:+.2f}",
                         va="center", ha="left" if v >= 0 else "right",
                         fontsize=7.5, color="0.3")
        axg.axvline(0, color="0.6", lw=0.8)
        axg.tick_params(length=0)
        for sp_ in ("top", "right", "left"):
            axg.spines[sp_].set_visible(False)
        axg.set_xticks([])
        if r == 0:
            axs.text(0.0, 1.10, "niche · where · condition (TEST ● vs CNTRL ●)",
                     fontsize=10, transform=axs.transAxes)
            axc.set_title("composition (within-niche %)", fontsize=10, loc="left")
            axg.set_title("condition-linked gene programs (ρ)", fontsize=10,
                          loc="left")
    fig.suptitle(f"Condition niche dossier — featured by {mode_label}\n"
                 f"{N_NOTE} · Δ = observed TEST−CNTRL mean (star: mixed-model FDR<0.05) "
                 f"· GP ρ per-cell Spearman (exploratory)", fontsize=12)
    fig.savefig(os.path.join(out_dir, "niche_dossier.pdf"), bbox_inches="tight")
    plt.close(fig)


def draw_overview(metric, featured, out_dir):
    feat = set(featured)
    if metric == "mean_score":
        d = df.sort_values("mean_score")
        xv = d["mean_score"].values
        xlabel = "Mean condition score"
        title = "All niches ranked by condition level"
    else:
        d = df.sort_values("abund_log2fc")
        xv = d["abund_log2fc"].values
        xlabel = "Abundance log2(TEST/CNTRL)"
        title = "All niches ranked by abundance shift"
    ids = d["niche"].tolist()
    colors = [COND_COLORS["TEST"] if l > 0
              else COND_COLORS["CNTRL"] for l in d["abund_log2fc"]]
    fig, ax = plt.subplots(figsize=(7, max(5, 0.22 * len(ids))))
    yy = list(range(len(ids)))
    ax.hlines(yy, 0, xv, color=colors, lw=2, alpha=0.85)
    ax.scatter(xv, yy, color=colors, s=28, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"N{n}" for n in ids], fontsize=7)
    for tick, n in zip(ax.get_yticklabels(), ids):
        if n in feat:
            tick.set_fontweight("bold"); tick.set_color("black")
        else:
            tick.set_color("0.6")
    if metric != "mean_score":
        ax.axvline(0, color="0.6", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_title(f"{title} · featured set in bold\n{N_NOTE}", fontsize=10)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    ax.scatter([], [], color=COND_COLORS["TEST"], label="TEST-enriched")
    ax.scatter([], [], color=COND_COLORS["CNTRL"], label="CNTRL-enriched")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(os.path.join(out_dir, "niche_overview_ranked.pdf"),
                bbox_inches="tight")
    plt.close(fig)


def run_mode(mode):
    pool = df[df["n_cells"] >= MIN_CELLS_NICHE]
    if mode == "condition":
        out, metric, label = OUTPUT_COND, "mean_score", "absolute condition level"
        featured = pool.sort_values("mean_score", ascending=False) \
            .head(N_FEATURE)["niche"].tolist()
    else:
        out, metric = OUTPUT_ABUND, "abund_log2fc"
        label = "abundance shift (both directions)"
        nh = N_FEATURE - N_FEATURE // 2
        test_up = pool.sort_values("abund_log2fc", ascending=False).head(nh)["niche"].tolist()
        cntrl_up = pool.sort_values("abund_log2fc", ascending=True).head(N_FEATURE // 2)["niche"].tolist()
        chosen = list(dict.fromkeys(test_up + cntrl_up))
        featured = pool[pool["niche"].isin(chosen)] \
            .sort_values("abund_log2fc", ascending=False)["niche"].tolist()
    os.makedirs(out, exist_ok=True)
    draw_dossier(featured, out, label)
    draw_overview(metric, featured, out)
    df.to_csv(os.path.join(out, "niche_dossier_table.csv"), index=False)
    print(f"  [{mode}] featured: {['N' + str(n) for n in featured]}  ->  {out}/")


print("\n[DOSSIER] Building both selections...")
run_mode("condition")
run_mode("abundance")
print(f"\n[DONE] niche_d complete. {N_NOTE}.")
print(f"  niche_dossier_condition/  (ranked by absolute condition)")
print(f"  niche_dossier_abundance/   (ranked by |abundance shift|)")