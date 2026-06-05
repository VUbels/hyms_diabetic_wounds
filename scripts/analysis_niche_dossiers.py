"""
==============================================================================
SPATIAL TRANSCRIPTOMICS PIPELINE — STEP 2.5d
Condition score niche dossier (focused, narrative)
==============================================================================

Purpose
-------
Instead of plotting every niche, this selects the handful of niches that carry
signal and lays each one out as a single row that reads left-to-right as the
biological chain:

    which niche  ·  is it HS- or LS-skewed, and how senescent (Ctrl vs Test)
        -> what cell types is it built from
            -> which gene programs track its condition score

It is run TWICE, writing to two separate folders so each story is self-contained:
  * niche_dossier_score/  — niches ranked by ABSOLUTE condition score level
                                 (the most senescent niches, regardless of which
                                 condition they sit in)
  * niche_dossier_abundance/   — niches ranked by ABUNDANCE SHIFT
                                 (|log2(Ctrl/Test)| — the niches that most change
                                 their prevalence between conditions)

Honesty
-------
With n=3 D vs n=3 L slides almost nothing is individually significant. The
condition score Δ (Ctrl−Test) shows an FDR star only when the cell-level mixed model
clears it; otherwise it reads "ns". Abundance is one observation per slide, so
its exact-MW p is floored at ~0.10 and is shown for transparency, not as proof.
GP correlations are per-cell Spearman and are exploratory. Everything here is an
effect-size / ranking view, labelled as such.

Inputs
------
    nichecompass_results/objects/nichecompass_integrated.h5ad   (or $NICHE_H5AD)
    shared caches in niche_analysis_shared/ (written by step 2.5a; rebuilt if
    absent so this can also run first)

Outputs (per folder)
--------------------
    niche_dossier.pdf          the featured-niche dossier
    niche_overview_ranked.pdf  all niches ranked by the selection metric
    niche_dossier_table.csv    per-niche numbers behind the figure
==============================================================================
"""

import argparse
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

# =============================================================================
# SHARED UTILITIES  (identical block across the three Step-2.5 scripts)
# =============================================================================
SHARED_DIR        = os.environ.get("NICHE_SHARED_DIR", "niche_analysis_shared")
NICHE_STYLE_CSV   = os.path.join(SHARED_DIR, "niche_style.csv")
CONDITION_SCORE_CSV = os.path.join(SHARED_DIR, "condition_score.csv")
CELLTYPE_STYLE_CANDIDATES = [
    "annotated_data/cell_type_colour_map.csv",        # reuse the R-pipeline map
    os.path.join(SHARED_DIR, "celltype_style.csv"),
]
os.makedirs(SHARED_DIR, exist_ok=True)

CONDITION_SCORE_KEY = "condition_score"
CONDITION_FLAG_KEY   = "is_high_condition"
CONDITION_SCORE_QUANTILE  = float(os.environ.get("CONDITION_SCORE_QUANTILE", "0.75"))

# condition score skin-specific condition score signature (HGNC symbols, de-duplicated).
# A handful of legacy aliases may not match current panel symbols
# Path to the comma-separated gene signature file.
# Each gene on a single line or all on one line, separated by commas.
# The scorer prints matched/unmatched genes so the list can be curated
# against the actual panel before use.
CONDITION_GENES_FILE = os.environ.get("CONDITION_GENES_FILE", None)

def _load_condition_genes(path):
    """Load a comma-separated gene list from a text file."""
    if path is None:
        return []
    try:
        with open(path) as _fh:
            raw = _fh.read()
        genes = [g.strip() for g in raw.replace("\n", ",").split(",") if g.strip()]
        print(f"  [SCORE] Loaded {len(genes)} genes from {path}")
        return genes
    except FileNotFoundError:
        print(f"  [SCORE][WARN] Gene file not found: {path} — no genes loaded.")
        return []

CONDITION_GENES: list[str] = _load_condition_genes(CONDITION_GENES_FILE)

# =============================================================================
# ARGUMENT PARSER
# =============================================================================
_parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_parser.add_argument(
    "--genes", default=None, metavar="FILE",
    help="Path to a comma-separated .txt file of gene symbols to score. "
         "Overrides CONDITION_GENES_FILE env var and the script default.")
_parser.add_argument(
    "--h5ad", default=None, metavar="FILE",
    help="Path to the integrated .h5ad file (default: "
         "nichecompass_results/objects/nichecompass_integrated.h5ad)")
_parser.add_argument(
    "--output_score", default=None, metavar="DIR",
    help="Output folder for the condition-score dossier (default: niche_dossier_score)")
_parser.add_argument(
    "--output_abund", default=None, metavar="DIR",
    help="Output folder for the abundance-shift dossier (default: niche_dossier_abundance)")
_parser.add_argument(
    "--exclude_samples", nargs="*", default=None, metavar="ID",
    help="Sample IDs to exclude, space-separated.")
_args = _parser.parse_args()

# Apply CLI values — take priority over env vars and script defaults
if _args.genes:
    CONDITION_GENES_FILE = _args.genes
    CONDITION_GENES[:] = _load_condition_genes(CONDITION_GENES_FILE)


# Condition labels — must match obs[CONDITION_KEY] values exactly.
# Edit these two lines; everything else derives from them.
CONDITION_CNTRL = "lesion"    # control group
CONDITION_TEST  = "diabetic"  # test group

CONDITION_ORDER = [CONDITION_CNTRL, CONDITION_TEST]
COND_COLORS = {CONDITION_CNTRL: "#C94A4A", CONDITION_TEST: "#4A7EB8"}
COND_SHORT      = {CONDITION_CNTRL: "Ctrl", CONDITION_TEST: "Test"}

# Editable text + sensible defaults for publication-grade vector output
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


def niche_sort_key(x):
    try:
        return (0, int(x))
    except Exception:
        return (1, str(x))


def pval_stars(p, ns=""):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ns
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ns


def _hierarchical_order(mat_df):
    """Leaf order of an average-linkage tree on rows (correlation distance)."""
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    rows = list(mat_df.index)
    if len(rows) <= 2:
        return rows
    X = np.nan_to_num(mat_df.values.astype(float))
    try:
        d = pdist(X, metric="correlation")
        d = np.nan_to_num(d, nan=1.0, posinf=1.0, neginf=1.0)
        order = leaves_list(linkage(d, method="average"))
        return [rows[i] for i in order]
    except Exception:
        return rows


def build_or_load_niche_style(adata, niche_key, celltype_key,
                              cache_csv=NICHE_STYLE_CSV, rebuild=False):
    """(color_dict, ordered_niches). Compositionally similar niches get adjacent
    hues; cached so every figure/condition colours the same niche identically."""
    niches_all = sorted(adata.obs[niche_key].astype(str).unique(), key=niche_sort_key)
    if (not rebuild) and os.path.exists(cache_csv):
        sty = pd.read_csv(cache_csv, dtype={"niche": str})
        color = dict(zip(sty["niche"], sty["color"]))
        if set(niches_all).issubset(set(color)):
            ordered = [n for n in sty.sort_values("order")["niche"].tolist()
                       if n in niches_all]
            ordered += [n for n in niches_all if n not in ordered]
            return color, ordered
    comp = pd.crosstab(adata.obs[niche_key].astype(str),
                       adata.obs[celltype_key], normalize="index").reindex(niches_all)
    ordered = _hierarchical_order(comp)
    hues = sns.color_palette("husl", n_colors=max(len(ordered), 3))
    color = {n: mpl.colors.to_hex(hues[i]) for i, n in enumerate(ordered)}
    pd.DataFrame({"niche": ordered, "order": range(len(ordered)),
                  "color": [color[n] for n in ordered]}).to_csv(cache_csv, index=False)
    return color, ordered


def build_or_load_celltype_colors(adata, celltype_key,
                                  candidates=CELLTYPE_STYLE_CANDIDATES):
    cts = sorted(map(str, adata.obs[celltype_key].unique()))
    color = {}
    for path in candidates:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                cc = {c.lower(): c for c in df.columns}
                kcol = cc.get("cell_type", cc.get("celltype"))
                vcol = cc.get("hex", cc.get("color"))
                if kcol and vcol:
                    color = dict(zip(df[kcol].astype(str), df[vcol].astype(str)))
                    break
            except Exception:
                pass
    missing = [c for c in cts if c not in color]
    if missing:
        pal = sns.color_palette("tab20", n_colors=max(len(missing), 3))
        for i, c in enumerate(missing):
            color[c] = mpl.colors.to_hex(pal[i % len(pal)])
        pd.DataFrame({"cell_type": list(color), "hex": list(color.values())}).to_csv(
            os.path.join(SHARED_DIR, "celltype_style.csv"), index=False)
    return {c: color[c] for c in cts}


def compute_condition_score(adata, gene_set=CONDITION_GENES,
                             score_key=CONDITION_SCORE_KEY,
                             cache_csv=CONDITION_SCORE_CSV, counts_layer="counts",
                             rebuild=False):
    """LogNorm -> sc.tl.score_genes(condition score). Cached per-cell so all scripts
    share identical scores. Returns the matched gene list."""
    var = set(map(str, adata.var_names))
    present = [g for g in gene_set if g in var]
    missing = [g for g in gene_set if g not in var]
    print(f"  [SCORE] {len(present)}/{len(gene_set)} condition score genes present; "
          f"{len(missing)} missing")
    if missing:
        print(f"  [SCORE] missing (curate aliases if needed): {', '.join(missing[:25])}"
              + (" ..." if len(missing) > 25 else ""))

    if (not rebuild) and os.path.exists(cache_csv):
        s = pd.read_csv(cache_csv, index_col=0)
        if score_key in s.columns and adata.obs_names.isin(s.index).all():
            adata.obs[score_key] = s.loc[adata.obs_names, score_key].values
            print(f"  [SCORE] loaded cached score from {cache_csv}")
            return present

    if len(present) < 5:
        print("  [SCORE][WARN] <5 signature genes present; condition score skipped.")
        adata.obs[score_key] = np.nan
        return present

    tmp = adata.copy()
    if counts_layer in tmp.layers:
        tmp.X = tmp.layers[counts_layer].copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.tl.score_genes(tmp, present, score_name=score_key, ctrl_size=50,
                      use_raw=False, random_state=0)
    adata.obs[score_key] = tmp.obs[score_key].values
    del tmp
    pd.DataFrame({score_key: adata.obs[score_key].values},
                 index=adata.obs_names).to_csv(cache_csv)
    with open(os.path.join(SHARED_DIR, "condition_genes_used.txt"), "w") as fh:
        fh.write("# matched (%d)\n%s\n\n# missing (%d)\n%s\n" %
                 (len(present), "\n".join(present), len(missing), "\n".join(missing)))
    print(f"  [SCORE] score computed and cached -> {cache_csv}")
    return present


def add_condition_flag(adata, score_key=CONDITION_SCORE_KEY,
                       flag_key=CONDITION_FLAG_KEY, quantile=CONDITION_SCORE_QUANTILE):
    if score_key not in adata.obs or adata.obs[score_key].isna().all():
        adata.obs[flag_key] = False
        return np.nan
    thr = float(np.nanquantile(adata.obs[score_key].values, quantile))
    adata.obs[flag_key] = (adata.obs[score_key].values >= thr)
    return thr


def mixedlm_condition(df, value_col, sample_col="sample_id",
                      condition_col="condition", extra_fixed=None,
                      hi=CONDITION_CNTRL, lo=CONDITION_TEST,
                      max_rows=30000, seed=0):
    """value ~ C(condition)[+extra] + (1|sample). Reference = Test, so the
    condition coefficient is the Ctrl-vs-Test effect (positive => higher in Ctrl).
    Returns dict(effect_Ctrl_vs_Test, p, n, converged)."""
    out = {"effect_Ctrl_vs_Test": np.nan, "p": np.nan, "n": int(len(df)),
           "converged": False}
    try:
        import statsmodels.formula.api as smf
    except Exception:
        return out
    cols = [value_col, condition_col, sample_col] + list(extra_fixed or [])
    d = df[cols].dropna().copy()
    if d[condition_col].nunique() < 2 or d[sample_col].nunique() < 3:
        out["n"] = int(len(d))
        return out
    if len(d) > max_rows:
        d = d.groupby(condition_col, group_keys=False).apply(
            lambda g: g.sample(n=min(len(g), max_rows // 2), random_state=seed))
    d[condition_col] = pd.Categorical(d[condition_col].astype(str), categories=[lo, hi])
    d["_grp"] = d[sample_col].astype(str)
    # Deterministic row order so any caller passing the same cells obtains an
    # identical fit (matters on near-degenerate problems).
    d = d.sort_values(["_grp", condition_col, value_col]).reset_index(drop=True)
    rhs = f"C({condition_col})"
    for c in (extra_fixed or []):
        rhs += f" + C({c})" if str(d[c].dtype) in ("object", "category") else f" + {c}"
    try:
        md = smf.mixedlm(f"{value_col} ~ {rhs}", d, groups=d["_grp"])
        r = md.fit(reml=False, method="lbfgs", maxiter=200, disp=False)
        # Guard against a singular / unidentifiable random effect: when the
        # between-sample variance collapses to ~0 the fixed-effect estimate sits
        # on a flat likelihood ridge and is not trustworthy, so report NaN
        # rather than a spurious effect (keeps results reproducible and honest).
        # On real replicated data this rarely triggers.
        try:
            re_var = float(np.asarray(r.cov_re).ravel()[0])
        except Exception:
            re_var = np.nan
        resid = float(getattr(r, "scale", np.nan))
        singular = ((not np.isfinite(re_var)) or re_var <= 0
                    or (np.isfinite(resid) and resid > 0 and re_var < 1e-6 * resid))
        cn = [c for c in r.params.index if c.startswith(f"C({condition_col})")]
        out["n"] = int(len(d))
        out["converged"] = bool(getattr(r, "converged", False))
        if cn and not singular and out["converged"] and np.isfinite(r.pvalues[cn[0]]):
            out.update(effect_Ctrl_vs_Test=float(r.params[cn[0]]),
                       p=float(r.pvalues[cn[0]]))
    except Exception:
        pass
    return out


def color_yticklabels(ax, labels, color_map):
    for tick, lab in zip(ax.get_yticklabels(), labels):
        tick.set_color(color_map.get(str(lab), "black"))

# =============================================================================
# CONFIG
# =============================================================================
INPUT_H5AD   = os.environ.get(
    "NICHE_H5AD", "nichecompass_results/objects/nichecompass_integrated.h5ad")
OUTPUT_SCORE = os.environ.get("DOSSIER_SCORE_DIR", "niche_dossier_score")
OUTPUT_ABUND = os.environ.get("DOSSIER_ABUND_DIR", "niche_dossier_abundance")

# CLI overrides (parser runs before this block)
if _args.h5ad:         INPUT_H5AD   = _args.h5ad
if _args.output_score: OUTPUT_SCORE = _args.output_score
if _args.output_abund: OUTPUT_ABUND = _args.output_abund
if _args.exclude_samples is not None: EXCLUDE_SAMPLES = _args.exclude_samples

CELLTYPE_KEY = "predicted_cell_type"
NICHE_KEY = "nichecompass_niche"
SAMPLE_KEY = "sample_id"
CONDITION_KEY = "condition"
LATENT_KEY = "nichecompass_latent"
ACTIVE_GP_KEY = "nichecompass_active_gp_names"

EXCLUDE_SAMPLES = [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]

N_FEATURE = int(os.environ.get("DOSSIER_N_FEATURE", "8"))
TOP_CELLTYPES = 5
TOP_GPS = 5
MIN_CELLS_NICHE = 30          # a niche must have >= this many cells to be featured
MIN_CELLS_GP = 20             # min cells to compute within-niche GP correlation
ABUND_PSEUDO = 1e-4           # pseudocount for abundance log2 fold change

# =============================================================================
# 1. LOAD + CONDITION SCORE + SHARED STYLES
# =============================================================================
print(f"[LOAD] {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
for k in (NICHE_KEY, CELLTYPE_KEY, CONDITION_KEY, SAMPLE_KEY):
    adata.obs[k] = adata.obs[k].astype(str)
if EXCLUDE_SAMPLES:
    keep = ~adata.obs[SAMPLE_KEY].isin(EXCLUDE_SAMPLES)
    print(f"  Excluding {EXCLUDE_SAMPLES}: dropping {(~keep).sum()} cells")
    adata = adata[keep].copy()

print("\n[SCORE] Scoring condition score condition score signature...")
compute_condition_score(adata)
add_condition_flag(adata)
niche_colors, _ = build_or_load_niche_style(adata, NICHE_KEY, CELLTYPE_KEY)
ct_colors = build_or_load_celltype_colors(adata, CELLTYPE_KEY)

all_niches = sorted(adata.obs[NICHE_KEY].unique(), key=niche_sort_key)
samp_cond = adata.obs.groupby(SAMPLE_KEY)[CONDITION_KEY].first().to_dict()
ctrl_samples = [s for s, c in samp_cond.items() if c == CONDITION_CNTRL]
test_samples = [s for s, c in samp_cond.items() if c == CONDITION_TEST]
n_ctrl, n_test = len(ctrl_samples), len(test_samples)
N_NOTE = f"n={n_ctrl} HS vs n={n_test} LS slides"

# latent / gene-program activities
gp_names = list(adata.uns.get(ACTIVE_GP_KEY, []))
L = np.asarray(adata.obsm[LATENT_KEY]) if LATENT_KEY in adata.obsm else None
if L is not None and (not gp_names or len(gp_names) != L.shape[1]):
    gp_names = [f"GP_{i}" for i in range(L.shape[1])]
HAS_GP = L is not None and len(gp_names) > 0


def shorten_gp(name, maxlen=18):
    s = str(name)
    base = s.split("_")[0].replace("COMPLEX:", "")
    if "ligand_receptor_target_gene" in s:
        tag = "tg"
    elif "ligand_receptor" in s:
        tag = "LR"
    elif "combined" in s:
        tag = "comb"
    elif s.startswith("Add-on"):
        tag = "addon"
    else:
        tag = ""
    if len(base) > maxlen - 5:
        base = base[:maxlen - 6] + "…"
    return f"{base} {tag}".strip()


def spearman_cols(M, y):
    """Spearman rho between each column of M (n x p) and vector y (n)."""
    yr = rankdata(y).astype(float)
    yr -= yr.mean()
    sy = yr.std()
    out = np.full(M.shape[1], np.nan)
    if sy < 1e-12:
        return out
    yr /= sy
    for j in range(M.shape[1]):
        xr = rankdata(M[:, j]).astype(float)
        sx = xr.std()
        if sx < 1e-12:
            continue
        xr = (xr - xr.mean()) / sx
        out[j] = float(np.mean(xr * yr))
    return out


# =============================================================================
# 2. PER-NICHE METRICS
# =============================================================================
print("\n[METRICS] Per-niche condition score / abundance / composition / GP...")
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
    s_hs = score[m & (cond_arr == CONDITION_CNTRL)]
    s_ls = score[m & (cond_arr == CONDITION_TEST)]
    mean_score = float(np.nanmean(s_all)) if n_cells else np.nan
    mean_ctrl = float(np.nanmean(s_hs)) if len(s_hs) else np.nan
    mean_test = float(np.nanmean(s_ls)) if len(s_ls) else np.nan

    res = mixedlm_condition(cell_df[cell_df["niche"] == niche], "score",
                            sample_col="sample")

    hs_props = prop_table.loc[niche, ctrl_samples].values.astype(float)
    ls_props = prop_table.loc[niche, test_samples].values.astype(float)
    mean_pH, mean_pL = float(hs_props.mean()), float(ls_props.mean())
    log2fc = float(np.log2((mean_pH + ABUND_PSEUDO) / (mean_pL + ABUND_PSEUDO)))
    mw_p = np.nan
    if len(hs_props) >= 1 and len(ls_props) >= 1 and \
            (hs_props.std() + ls_props.std()) > 0:
        try:
            mw_p = float(mannwhitneyu(hs_props, ls_props,
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
         "mean_score_Ctrl": mean_ctrl, "mean_score_Test": mean_test,
         "score_delta_obs": (mean_ctrl - mean_test) if np.isfinite(mean_ctrl)
         and np.isfinite(mean_test) else np.nan,
         "score_effect_mixed": res["effect_Ctrl_vs_Test"], "score_p": res["p"],
         "mean_prop_Ctrl": mean_pH, "mean_prop_Test": mean_pL,
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

# =============================================================================
# 3. FIGURES
# =============================================================================
def draw_dossier(featured, out_dir, mode_label):
    K = len(featured)
    if K == 0:
        print(f"  [WARN] no niches to feature for {mode_label}")
        return
    vmax = max([max(by_niche[n]["mean_score_Ctrl"] if np.isfinite(by_niche[n]["mean_score_Ctrl"]) else 0,
                    by_niche[n]["mean_score_Test"] if np.isfinite(by_niche[n]["mean_score_Test"]) else 0)
                for n in featured] + [1e-3]) * 1.18
    fig = plt.figure(figsize=(13, 1.3 * K + 1.3))
    gs = fig.add_gridspec(K, 3, width_ratios=[1.15, 1.25, 1.25], wspace=0.45,
                          hspace=0.7, top=0.9, bottom=0.05, left=0.04, right=0.985)
    for r, nid in enumerate(featured):
        d = by_niche[nid]
        # --- status column ---
        axs = fig.add_subplot(gs[r, 0]); axs.axis("off")
        axs.set_xlim(0, 1); axs.set_ylim(0, 1)
        axs.add_patch(mpatches.Rectangle((0.0, 0.80), 0.05, 0.17,
                      color=niche_colors.get(nid, "#999"), clip_on=False,
                      transform=axs.transAxes))
        axs.text(0.08, 0.885, f"N{nid}", fontsize=13, fontweight="bold",
                 va="center", transform=axs.transAxes)
        is_hs = d["abund_log2fc"] > 0
        pcol = COND_COLORS[CONDITION_CNTRL] if is_hs else COND_COLORS[CONDITION_TEST]
        axs.text(0.0, 0.60, f"{COND_SHORT[CONDITION_CNTRL] if is_hs else COND_SHORT[CONDITION_TEST]}-enriched · "
                 f"log2FC {d['abund_log2fc']:+.1f}", fontsize=9, color=pcol,
                 transform=axs.transAxes)
        y, x0, w = 0.33, 0.0, 0.60
        axs.plot([x0, x0 + w], [y, y], color="0.85", lw=1, transform=axs.transAxes)
        mh, ml = d["mean_score_Ctrl"], d["mean_score_Test"]
        if np.isfinite(mh) and np.isfinite(ml):
            xh, xl = x0 + w * mh / vmax, x0 + w * ml / vmax
            axs.plot([xl, xh], [y, y], color="0.55", lw=2, transform=axs.transAxes)
        for val, c in [(ml, COND_COLORS[CONDITION_TEST]),
                       (mh, COND_COLORS[CONDITION_CNTRL])]:
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
        # --- composition column ---
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
        # --- gene-program column ---
        axg = fig.add_subplot(gs[r, 2])
        gps = d["gp"][:TOP_GPS][::-1]
        gn = [shorten_gp(g) for g, _ in gps]; gr = [v for _, v in gps]
        if gr:
            axg.barh(range(len(gn)), gr, edgecolor="white", linewidth=0.5,
                     color=[COND_COLORS[CONDITION_CNTRL] if v > 0
                            else COND_COLORS[CONDITION_TEST] for v in gr])
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
            axs.text(0.0, 1.10, "niche · where · condition score (Ctrl ● vs Test ●)",
                     fontsize=10, transform=axs.transAxes)
            axc.set_title("composition (within-niche %)", fontsize=10, loc="left")
            axg.set_title("condition-score-linked gene programs (ρ)", fontsize=10,
                          loc="left")
    fig.suptitle(f"Condition score niche dossier — featured by {mode_label}\n"
                 f"{N_NOTE} · Δ = observed Ctrl−Test mean (star: mixed-model FDR<0.05) "
                 f"· GP ρ per-cell Spearman (exploratory)", fontsize=12)
    fig.savefig(os.path.join(out_dir, "niche_dossier.pdf"), bbox_inches="tight")
    plt.close(fig)


def draw_overview(metric, featured, out_dir):
    feat = set(featured)
    if metric == "mean_score":
        d = df.sort_values("mean_score")
        xv = d["mean_score"].values
        xlabel = "Mean condition score"
        title = "All niches ranked by condition score level"
    else:
        d = df.sort_values("abund_log2fc")
        xv = d["abund_log2fc"].values
        xlabel = "Abundance log2(Ctrl/Test)"
        title = "All niches ranked by abundance shift"
    ids = d["niche"].tolist()
    colors = [COND_COLORS[CONDITION_CNTRL] if l > 0
              else COND_COLORS[CONDITION_TEST] for l in d["abund_log2fc"]]
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
    ax.scatter([], [], color=COND_COLORS[CONDITION_CNTRL], label="Ctrl-enriched")
    ax.scatter([], [], color=COND_COLORS[CONDITION_TEST], label="Test-enriched")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(os.path.join(out_dir, "niche_overview_ranked.pdf"),
                bbox_inches="tight")
    plt.close(fig)


def run_mode(mode):
    pool = df[df["n_cells"] >= MIN_CELLS_NICHE]
    if mode == "condition_score":
        out, metric, label = OUTPUT_SCORE, "mean_score", "absolute condition score level"
        featured = pool.sort_values("mean_score", ascending=False) \
            .head(N_FEATURE)["niche"].tolist()
    else:
        out, metric = OUTPUT_ABUND, "abund_log2fc"
        label = "abundance shift (both directions)"
        nh = N_FEATURE - N_FEATURE // 2
        hs_up = pool.sort_values("abund_log2fc", ascending=False).head(nh)["niche"].tolist()
        ls_up = pool.sort_values("abund_log2fc", ascending=True).head(N_FEATURE // 2)["niche"].tolist()
        chosen = list(dict.fromkeys(hs_up + ls_up))
        featured = pool[pool["niche"].isin(chosen)] \
            .sort_values("abund_log2fc", ascending=False)["niche"].tolist()
    os.makedirs(out, exist_ok=True)
    draw_dossier(featured, out, label)
    draw_overview(metric, featured, out)
    df.to_csv(os.path.join(out, "niche_dossier_table.csv"), index=False)
    print(f"  [{mode}] featured: {['N' + str(n) for n in featured]}  ->  {out}/")


print("\n[DOSSIER] Building both selections...")
run_mode("condition_score")
run_mode("abundance")
print(f"\n[DONE] Step 2.5d complete. {N_NOTE}.")
print(f"  niche_dossier_score/  (ranked by absolute condition score)")
print(f"  niche_dossier_abundance/   (ranked by |abundance shift|)")