"""
==============================================================================
SPATIAL TRANSCRIPTOMICS PIPELINE — STEP 2.5a
Niche composition + condition score
==============================================================================

Overview
--------
A fully general niche-composition and condition-score analysis script.
"Condition" is defined by the user via a gene signature (CONDITION_GENES) and
a binary grouping in obs[CONDITION_KEY] (e.g. "control" vs "test"). Nothing
in this script is hard-wired to a particular biology; swap CONDITION_GENES and
the condition labels to reuse for any scored contrast.

Statistics summary
------------------
  * Condition score (sc.tl.score_genes on CONDITION_GENES) computed once and
    cached in SHARED_DIR so all Step-2.5 scripts share identical values.
  * Exact Mann-Whitney (not the normal-approximation ranksums) for the
    sample-level proportion tests — the only honest test when n_samples is
    small (e.g. 3 vs 2).
  * Linear mixed models (statsmodels) with SAMPLE as a random intercept:
        - per cell type  : prop  ~ C(condition) + C(niche) + (1|sample)
        - per niche      : score ~ C(condition) + (1|sample)  [cell-level]
        - per cell type  : score ~ C(condition) + (1|sample)  [cell-level]
    Cell-level condition mixed models are properly powered (many cells per
    sample) and partition out between-sample variance — this is where
    defensible significance comes from given limited biological replicates.
  * Shared, similarity-ordered niche colour map (compositionally similar niches
    get adjacent hues) and reuse of the R pipeline's cell-type colour map, so
    the same niche / cell type is the same colour in every figure.
  * Significance markers (*, **, ***) on every comparison where a test
    supports them; captions state the test and n explicitly.

Inputs
------
    nichecompass_results/objects/nichecompass_integrated.h5ad   (or $NICHE_H5AD)

Key outputs (niche_composition/)
---------------------------------
    condition_score_by_celltype.pdf / .csv
    condition_score_by_niche.pdf / .csv
    condition_score_vs_composition_shift.pdf / .csv
    condition_score_gp_association.pdf / .csv
    niche_composition_combined.pdf
    differential.pdf  (+ niche_composition_differential.csv, mixedlm_*.csv)
    niche_proportions_per_sample.pdf
==============================================================================
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

# =============================================================================
# SHARED UTILITIES  (identical block across the three Step-2.5 scripts)
# =============================================================================
SHARED_DIR      = os.environ.get("NICHE_SHARED_DIR", "niche_analysis_shared")
NICHE_STYLE_CSV = os.path.join(SHARED_DIR, "niche_style.csv")
CONDITION_SCORE_CSV = os.path.join(SHARED_DIR, "condition_score.csv")
CELLTYPE_STYLE_CANDIDATES = [
    "annotated_data/cell_type_colour_map.csv",        # reuse the R-pipeline map
    os.path.join(SHARED_DIR, "celltype_style.csv"),
]
os.makedirs(SHARED_DIR, exist_ok=True)

# ── Score / flag keys ────────────────────────────────────────────────────────
CONDITION_SCORE_KEY    = "condition_score"   # obs column for the scored signature
CONDITION_FLAG_KEY     = "is_control"        # obs column for the top-quantile flag
CONDITION_SCORE_QUANTILE = float(os.environ.get("CONDITION_SCORE_QUANTILE", "0.75"))

# ── Gene signature ────────────────────────────────────────────────────────────
# Path to a comma-separated .txt file of gene symbols.
# Genes may be on one line or multiple lines; the loader handles both.
# Override at runtime: --genes /path/to/file.txt
CONDITION_GENES_FILE = os.environ.get(
    "CONDITION_GENES_FILE",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "condition_genes.txt"),
)

def _load_condition_genes(path):
    """Load a comma-separated gene list from a text file."""
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
    "--output_dir", default=None, metavar="DIR",
    help="Output directory (default: niche_composition)")
_parser.add_argument(
    "--exclude_samples", nargs="*", default=None, metavar="ID",
    help="Sample IDs to exclude, space-separated.")
_args = _parser.parse_args()

# Apply CLI values — take priority over env vars and script defaults
if _args.genes:
    CONDITION_GENES_FILE = _args.genes
    CONDITION_GENES[:] = _load_condition_genes(CONDITION_GENES_FILE)


# ── Condition labels ──────────────────────────────────────────────────────────
# Must match the values in obs[CONDITION_KEY] exactly.
CONDITION_CNTRL    = "lesion"
CONDITION_TEST    = "diabetic"
CONDITION_ORDER = [CONDITION_CNTRL, CONDITION_TEST]
COND_COLORS     = {CONDITION_CNTRL: "#C94A4A", CONDITION_TEST: "#4A7EB8"}
COND_SHORT      = {CONDITION_CNTRL: "Ctrl",    CONDITION_TEST: "Test"}

# ── Publication-grade vector output ──────────────────────────────────────────
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


# ─────────────────────────────────────────────────────────────────────────────
# Utility functions
# ─────────────────────────────────────────────────────────────────────────────

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
    """Return (color_dict, ordered_niches). Compositionally similar niches get
    adjacent hues; result is cached so every figure colours the same niche
    identically across conditions and scripts."""
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
                             cache_csv=CONDITION_SCORE_CSV,
                             counts_layer="counts", rebuild=False):
    """LogNorm → sc.tl.score_genes(gene_set). Cached per-cell so all scripts
    share identical scores. Returns the matched gene list."""
    var = set(map(str, adata.var_names))
    present = [g for g in gene_set if g in var]
    missing = [g for g in gene_set if g not in var]
    print(f"  [SCORE] {len(present)}/{len(gene_set)} signature genes present; "
          f"{len(missing)} missing")
    if missing:
        print(f"  [SCORE] missing (curate aliases if needed): "
              f"{', '.join(missing[:25])}" + (" ..." if len(missing) > 25 else ""))

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
                       flag_key=CONDITION_FLAG_KEY,
                       quantile=CONDITION_SCORE_QUANTILE):
    """Flag the top-quantile cells as 'high condition'. Returns the threshold."""
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
    """Fit: value ~ C(condition)[+extra] + (1|sample).
    Reference level = LO (test), so the condition coefficient is the
    Ctrl-vs-Test effect (positive => higher in control).
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
    d[condition_col] = pd.Categorical(d[condition_col].astype(str),
                                      categories=[lo, hi])
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
        try:
            re_var = float(np.asarray(r.cov_re).ravel()[0])
        except Exception:
            re_var = np.nan
        resid = float(getattr(r, "scale", np.nan))
        singular = ((not np.isfinite(re_var)) or re_var <= 0
                    or (np.isfinite(resid) and resid > 0
                        and re_var < 1e-6 * resid))
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
INPUT_H5AD  = os.environ.get(
    "NICHE_H5AD", "nichecompass_results/objects/nichecompass_integrated.h5ad")
OUTPUT_DIR  = "niche_composition"

# CLI overrides (parser runs before this block)
if _args.h5ad:           INPUT_H5AD  = _args.h5ad
if _args.output_dir:     OUTPUT_DIR  = _args.output_dir
if _args.exclude_samples is not None: EXCLUDE_SAMPLES = _args.exclude_samples

SPATIAL_KEY   = "X_spatial"
CELLTYPE_KEY  = "predicted_cell_type"
NICHE_KEY     = "nichecompass_niche"
SAMPLE_KEY    = "sample_id"
CONDITION_KEY = "condition"
LATENT_KEY    = "nichecompass_latent"
ACTIVE_GP_KEY = "nichecompass_active_gp_names"

# Optionally exclude specific samples, e.g. EXCLUDE_SAMPLES=["Sample_X"]
EXCLUDE_SAMPLES = [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]

NEIGHBOR_K      = 15    # k for local cell-type abundance neighbourhood
MIN_CELLS_NICHE = 10    # ignore a niche in a subset below this many cells

os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# 1. LOAD
# =============================================================================
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

adata.obs[NICHE_KEY]     = adata.obs[NICHE_KEY].astype(str)
adata.obs[CELLTYPE_KEY]  = adata.obs[CELLTYPE_KEY].astype(str)
adata.obs[CONDITION_KEY] = adata.obs[CONDITION_KEY].astype(str)
adata.obs[SAMPLE_KEY]    = adata.obs[SAMPLE_KEY].astype(str)

all_niches  = sorted(adata.obs[NICHE_KEY].unique(), key=niche_sort_key)
all_cts     = sorted(adata.obs[CELLTYPE_KEY].unique())
all_samples = sorted(adata.obs[SAMPLE_KEY].unique())
sample_cond_map = adata.obs.groupby(SAMPLE_KEY)[CONDITION_KEY].first().to_dict()
n_ctrl = sum(1 for v in sample_cond_map.values() if v == CONDITION_CNTRL)
n_test = sum(1 for v in sample_cond_map.values() if v == CONDITION_TEST)

print(f"  {adata.n_obs} cells | {len(all_niches)} niches | {len(all_cts)} cell types")
print(f"  Samples: {len(all_samples)}  ({COND_SHORT[CONDITION_CNTRL]}={n_ctrl}, "
      f"{COND_SHORT[CONDITION_TEST]}={n_test})")
N_NOTE = (f"n={n_ctrl} {COND_SHORT[CONDITION_CNTRL]} vs "
          f"n={n_test} {COND_SHORT[CONDITION_TEST]} slides")

# =============================================================================
# 2. CONDITION SCORE + STYLES
# =============================================================================
print("\n[SCORE] Computing condition gene signature score...")
score_genes_used = compute_condition_score(adata)
score_thr = add_condition_flag(adata)
HAS_SCORE = (CONDITION_SCORE_KEY in adata.obs
             and not adata.obs[CONDITION_SCORE_KEY].isna().all())
if HAS_SCORE:
    print(f"  High-condition threshold "
          f"(q{int(CONDITION_SCORE_QUANTILE * 100)}) = {score_thr:.4f}; "
          f"{int(adata.obs[CONDITION_FLAG_KEY].sum())} high-condition cells "
          f"({100 * adata.obs[CONDITION_FLAG_KEY].mean():.1f}%)")

niche_colors, niche_order = build_or_load_niche_style(adata, NICHE_KEY, CELLTYPE_KEY)
ct_colors = build_or_load_celltype_colors(adata, CELLTYPE_KEY)
# niche_order is the global similarity order; restrict to those present here
niche_order = [n for n in niche_order if n in all_niches]

# GP activity matrix (signed) for condition-score association
gp_df = None
if HAS_SCORE and LATENT_KEY in adata.obsm:
    latent   = np.asarray(adata.obsm[LATENT_KEY])
    gp_names = list(adata.uns.get(ACTIVE_GP_KEY, []))
    if len(gp_names) == latent.shape[1]:
        gp_df = pd.DataFrame(latent, index=adata.obs_names, columns=gp_names)
    else:
        gp_df = pd.DataFrame(latent, index=adata.obs_names,
                             columns=[f"latent_{i}" for i in range(latent.shape[1])])
        print(f"  [GP][WARN] active GP names ({len(gp_names)}) != latent dim "
              f"({latent.shape[1]}); using generic latent labels")

# =============================================================================
# 3. LOCAL CELL TYPE ABUNDANCE (per cell, self EXCLUDED)
# =============================================================================
print(f"\n[NEIGHBORS] Local cell-type abundance (k={NEIGHBOR_K}, self-excluded)...")
abundance = np.zeros((adata.n_obs, len(all_cts)))
ct_to_idx = {ct: i for i, ct in enumerate(all_cts)}
for slide in all_samples:
    mask  = (adata.obs[SAMPLE_KEY] == slide).values
    idx   = np.where(mask)[0]
    coords    = adata.obsm[SPATIAL_KEY][mask]
    ct_labels = adata.obs[CELLTYPE_KEY].values[mask]
    if len(coords) < 2:
        continue
    tree = cKDTree(coords)
    k    = min(NEIGHBOR_K + 1, len(coords))    # +1 for self
    _, nbrs = tree.query(coords, k=k)
    if nbrs.ndim == 1:
        nbrs = nbrs[:, None]
    nbrs = nbrs[:, 1:]                          # drop self
    for i, neigh in enumerate(nbrs):
        for ct in ct_labels[neigh]:
            abundance[idx[i], ct_to_idx[ct]] += 1
    denom = max(nbrs.shape[1], 1)
    abundance[idx, :] /= denom
abundance_df = pd.DataFrame(abundance, columns=all_cts, index=adata.obs_names)

# =============================================================================
# 4. PER-CONDITION NICHE × CELLTYPE PROPORTION + CORRELATION
# =============================================================================
def compute_niche_ct_panel(adata_sub, abundance_sub, niches_all, cts_all):
    niches = adata_sub.obs[NICHE_KEY].values
    cts    = adata_sub.obs[CELLTYPE_KEY].values
    prop   = pd.DataFrame(0.0,      index=niches_all, columns=cts_all)
    corr   = pd.DataFrame(np.nan,   index=niches_all, columns=cts_all)
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

# =============================================================================
# 5. NICHE ENRICHMENT BY CONDITION (sample-level, EXACT Mann-Whitney)
# =============================================================================
print(f"\n[ENRICH] Niche abundance {COND_SHORT[CONDITION_CNTRL]} vs "
      f"{COND_SHORT[CONDITION_TEST]} (exact Mann-Whitney, sample-level)...")
niche_by_sample = pd.crosstab(adata.obs[NICHE_KEY], adata.obs[SAMPLE_KEY],
                               normalize="columns").reindex(all_niches)
niche_by_sample.to_csv(os.path.join(OUTPUT_DIR, "niche_sample_frequency.csv"))

enr_rows = []
for niche in all_niches:
    props = niche_by_sample.loc[niche]
    ctrl = np.array([props[s] for s in props.index
                     if sample_cond_map[s] == CONDITION_CNTRL])
    test = np.array([props[s] for s in props.index
                     if sample_cond_map[s] == CONDITION_TEST])
    mctrl = ctrl.mean() if len(ctrl) else 0.0
    mtest = test.mean() if len(test) else 0.0
    p = np.nan
    if (len(ctrl) >= 1 and len(test) >= 1
            and (ctrl.std() > 0 or test.std() > 0)
            and len(ctrl) + len(test) >= 4):
        try:
            p = mannwhitneyu(ctrl, test, alternative="two-sided", method="exact")[1]
        except Exception:
            p = np.nan
    enr_rows.append({
        "niche": niche,
        f"mean_prop_{COND_SHORT[CONDITION_CNTRL]}": mctrl,
        f"mean_prop_{COND_SHORT[CONDITION_TEST]}": mtest,
        f"sd_{COND_SHORT[CONDITION_CNTRL]}": ctrl.std() if len(ctrl) > 1 else np.nan,
        f"sd_{COND_SHORT[CONDITION_TEST]}": test.std() if len(test) > 1 else np.nan,
        f"log2FC_{COND_SHORT[CONDITION_CNTRL]}_vs_{COND_SHORT[CONDITION_TEST]}":
            np.log2((mctrl + 1e-6) / (mtest + 1e-6)),
        "mw_p": p,
        f"n_{COND_SHORT[CONDITION_CNTRL]}": len(ctrl),
        f"n_{COND_SHORT[CONDITION_TEST]}": len(test),
    })
enrichment_df = pd.DataFrame(enr_rows)
# Alias used internally for log2FC lookups
enrichment_df["log2FC_Ctrl_vs_Test"] = enrichment_df[
    f"log2FC_{COND_SHORT[CONDITION_CNTRL]}_vs_{COND_SHORT[CONDITION_TEST]}"]
vp = enrichment_df["mw_p"].dropna()
enrichment_df["mw_padj"] = np.nan
if len(vp) > 0:
    enrichment_df.loc[vp.index, "mw_padj"] = multipletests(vp, method="fdr_bh")[1]
enrichment_df.to_csv(
    os.path.join(OUTPUT_DIR, "niche_condition_enrichment.csv"), index=False)

# =============================================================================
# 6. DIFFERENTIAL NICHE × CELLTYPE COMPOSITION
#    - sample-level exact Mann-Whitney (transparent, honest floor)
#    - per-cell-type linear mixed model (sample random intercept): the
#      principled, pooled test of a global compositional shift per cell type.
# =============================================================================
print("\n[DIFF] Niche × cell-type composition: exact MW + mixed model...")

# Per (sample, niche) within-niche cell-type proportions
rows = []
for sample in all_samples:
    sdat = adata.obs.loc[adata.obs[SAMPLE_KEY] == sample,
                         [NICHE_KEY, CELLTYPE_KEY]]
    for niche in all_niches:
        nd   = sdat[sdat[NICHE_KEY] == niche]
        n_in = len(nd)
        for ct in all_cts:
            rows.append({
                "sample":    sample,
                "condition": sample_cond_map[sample],
                "niche":     niche,
                "celltype":  ct,
                "prop":      (nd[CELLTYPE_KEY] == ct).sum() / max(n_in, 1),
                "n_in_niche": n_in,
            })
sn_ct = pd.DataFrame(rows)
sn_ct.to_csv(
    os.path.join(OUTPUT_DIR, "niche_celltype_proportions_per_sample.csv"),
    index=False)

# (a) Per niche × celltype exact Mann-Whitney
diff_rows = []
for niche in all_niches:
    for ct in all_cts:
        d  = sn_ct[(sn_ct["niche"] == niche) & (sn_ct["celltype"] == ct)]
        ctrl = d[d["condition"] == CONDITION_CNTRL]["prop"].values
        test = d[d["condition"] == CONDITION_TEST]["prop"].values
        mctrl = ctrl.mean() if len(ctrl) else 0
        mtest = test.mean() if len(test) else 0
        p = np.nan
        if len(ctrl) + len(test) >= 4 and (np.std(ctrl) > 0 or np.std(test) > 0):
            try:
                p = mannwhitneyu(ctrl, test, alternative="two-sided", method="exact")[1]
            except Exception:
                p = np.nan
        diff_rows.append({
            "niche":               niche,
            "celltype":            ct,
            "mean_prop_Ctrl":      mctrl,
            "mean_prop_Test":      mtest,
            "log2FC_Ctrl_vs_Test": np.log2((mctrl + 1e-6) / (mtest + 1e-6)),
            "mw_p": p,
        })
diff_df = pd.DataFrame(diff_rows)
diff_df["mw_padj"] = np.nan
for niche in all_niches:                           # FDR corrected within niche
    m = (diff_df["niche"] == niche) & diff_df["mw_p"].notna()
    if m.sum() > 0:
        diff_df.loc[m, "mw_padj"] = multipletests(
            diff_df.loc[m, "mw_p"], method="fdr_bh")[1]

# (b) Per cell-type mixed model: prop ~ C(condition) + C(niche) + (1|sample)
ct_mm_rows = []
for ct in all_cts:
    d   = sn_ct[sn_ct["celltype"] == ct]
    res = mixedlm_condition(d, "prop", sample_col="sample", extra_fixed=["niche"])
    ct_mm_rows.append({"celltype": ct, **res})
ct_mm = pd.DataFrame(ct_mm_rows)
vp = ct_mm["p"].dropna()
ct_mm["padj"] = np.nan
if len(vp) > 0:
    ct_mm.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
ct_mm.to_csv(
    os.path.join(OUTPUT_DIR, "mixedlm_celltype_global_shift.csv"), index=False)
diff_df = diff_df.merge(
    ct_mm[["celltype", "p", "padj"]].rename(
        columns={"p": "celltype_mm_p", "padj": "celltype_mm_padj"}),
    on="celltype", how="left")
diff_df.to_csv(
    os.path.join(OUTPUT_DIR, "niche_composition_differential.csv"), index=False)
ct_sig = {r.celltype: r.padj for r in ct_mm.itertuples()}

# =============================================================================
# 7. CONDITION SCORE: per cell type and per niche (cell-level mixed models)
# =============================================================================
score_ct = score_niche = None
gp_score = None
if HAS_SCORE:
    print("\n[SCORE] Condition score by cell type and niche "
          "(cell-level mixed models)...")
    cell_df = adata.obs[[CONDITION_SCORE_KEY, CELLTYPE_KEY, NICHE_KEY,
                         SAMPLE_KEY, CONDITION_KEY, CONDITION_FLAG_KEY]].copy()
    cell_df.columns = ["score", "celltype", "niche",
                       "sample", "condition", "is_control"]

    # ── Per cell type ─────────────────────────────────────────────────────────
    rows = []
    for ct in all_cts:
        d   = cell_df[cell_df["celltype"] == ct]
        res = mixedlm_condition(d, "score", sample_col="sample")
        rows.append({
            "celltype":        ct,
            "mean_score":      d["score"].mean(),
            "mean_score_Ctrl": d[d.condition == CONDITION_CNTRL]["score"].mean(),
            "mean_score_Test": d[d.condition == CONDITION_TEST]["score"].mean(),
            "frac_control":    d["is_control"].mean(),
            "n_cells":         len(d),
            **res,
        })
    score_ct = pd.DataFrame(rows)
    vp = score_ct["p"].dropna()
    score_ct["padj"] = np.nan
    if len(vp) > 0:
        score_ct.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
    score_ct = score_ct.sort_values("mean_score", ascending=False)
    score_ct.to_csv(
        os.path.join(OUTPUT_DIR, "condition_score_by_celltype.csv"), index=False)

    # ── Per niche ─────────────────────────────────────────────────────────────
    rows = []
    for niche in all_niches:
        d   = cell_df[cell_df["niche"] == niche]
        res = mixedlm_condition(d, "score", sample_col="sample")
        rows.append({
            "niche":           niche,
            "mean_score":      d["score"].mean(),
            "mean_score_Ctrl": d[d.condition == CONDITION_CNTRL]["score"].mean(),
            "mean_score_Test": d[d.condition == CONDITION_TEST]["score"].mean(),
            "frac_control":    d["is_control"].mean(),
            "n_cells":         len(d),
            **res,
        })
    score_niche = pd.DataFrame(rows)
    vp = score_niche["p"].dropna()
    score_niche["padj"] = np.nan
    if len(vp) > 0:
        score_niche.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
    score_niche.to_csv(
        os.path.join(OUTPUT_DIR, "condition_score_by_niche.csv"), index=False)

    # ── GP association (per-cell Spearman; exploratory ranking) ───────────────
    if gp_df is not None:
        s = adata.obs[CONDITION_SCORE_KEY].values
        rows = []
        for gp in gp_df.columns:
            rho, pv = spearmanr(gp_df[gp].values, s)
            rows.append({"gene_program": gp, "spearman_rho": rho, "p": pv})
        gp_score = pd.DataFrame(rows).dropna(subset=["spearman_rho"])
        if len(gp_score):
            gp_score["padj"] = multipletests(
                gp_score["p"].fillna(1.0), method="fdr_bh")[1]
            gp_score = gp_score.reindex(
                gp_score["spearman_rho"].abs()
                .sort_values(ascending=False).index)
        gp_score.to_csv(
            os.path.join(OUTPUT_DIR, "condition_score_gp_association.csv"),
            index=False)

# =============================================================================
# 8. PLOTS
# =============================================================================
def _legend_conditions(fig, loc="upper right"):
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markersize=8,
                   markerfacecolor=COND_COLORS[c], label=COND_SHORT[c])
        for c in CONDITION_ORDER
    ]
    fig.legend(handles=handles, loc=loc, frameon=False)


# ── 8.1 Combined Ctrl/Test composition dot plot ───────────────────────────────
def plot_combined_composition(out_path, figsize=(11, 7)):
    conds = [c for c in CONDITION_ORDER if c in per_condition]
    if not conds:
        return
    niches_info = [n for n in niche_order
                   if any(per_condition[c]["prop"].loc[n].sum() > 0
                          for c in conds)]
    if not niches_info:
        return
    prop_ref = per_condition[conds[0]]["prop"].reindex(niches_info)
    ct_order = prop_ref.max(axis=0).sort_values(ascending=False).index.tolist()
    n_n, n_c = len(niches_info), len(ct_order)

    fig = plt.figure(figsize=figsize)
    gs  = fig.add_gridspec(
        2, len(conds) + 1,
        width_ratios=[n_c] * len(conds) + [1.4],
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
        ax.set_xticklabels(
            [f"{ct}{(' ' + pval_stars(ct_sig.get(ct, np.nan))).rstrip()}"
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

    # Differential column: niche enrichment Ctrl/Test, stars from exact MW
    axd = axes[-1]
    for yi, niche in enumerate(niches_info):
        row = enrichment_df[enrichment_df["niche"] == niche]
        if len(row) == 0:
            continue
        lfc  = row["log2FC_Ctrl_vs_Test"].iloc[0]
        padj = row["mw_padj"].iloc[0]
        if np.isfinite(lfc):
            sig = (not np.isnan(padj)) and padj < 0.05
            axd.scatter(
                0, yi,
                s=min(abs(lfc) * 90, 320) if sig else max(abs(lfc) * 45, 18),
                c=[lfc], cmap="RdBu_r", vmin=-2, vmax=2,
                edgecolors="black" if sig else "gray",
                linewidths=1.2 if sig else 0.3)
            star = pval_stars(padj)
            if star:
                axd.text(0.32, yi, star, va="center", fontsize=9)
    axd.set_xticks([0])
    axd.set_xticklabels([f"{COND_SHORT[CONDITION_CNTRL]}/{COND_SHORT[CONDITION_TEST]}"])
    axd.set_xlim(-0.5, 0.6)
    axd.set_ylim(-0.5, n_n - 0.5)
    axd.invert_yaxis()
    axd.grid(alpha=0.2, ls="--", lw=0.3)
    axd.set_title("Niche\nabundance", fontsize=10)
    plt.setp(axd.get_yticklabels(), visible=False)

    cax = fig.add_subplot(gs[1, :len(conds)])
    sm  = mpl.cm.ScalarMappable(cmap="RdBu_r",
                                norm=mpl.colors.Normalize(vmin=-0.5, vmax=0.5))
    fig.colorbar(sm, cax=cax, orientation="horizontal",
                 label="Pearson r (local cell-type abundance vs niche membership)")
    fig.suptitle(
        f"Niche composition by condition  ·  dot size = within-niche "
        f"proportion  ·  {N_NOTE}",
        fontsize=12, y=1.02)
    plt.savefig(out_path)
    plt.close()


print("\n[PLOT] composition dot plot...")
plot_combined_composition(
    os.path.join(OUTPUT_DIR, "niche_composition_combined.pdf"))

# ── 8.2 Differential composition heatmap (log2FC, celltype-level sig) ─────────
print("[PLOT] differential composition heatmap...")
lfc_piv  = diff_df.pivot(index="niche", columns="celltype",
                         values="log2FC_Ctrl_vs_Test").reindex(niche_order)[all_cts]
padj_piv = diff_df.pivot(index="niche", columns="celltype",
                         values="mw_padj").reindex(niche_order)[all_cts]
annot = padj_piv.apply(lambda col: col.map(lambda p: pval_stars(p)))
fig, ax = plt.subplots(figsize=(max(5, len(all_cts) * 0.34),
                                max(3, len(niche_order) * 0.30)))
sns.heatmap(lfc_piv, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            annot=annot.values, fmt="", annot_kws={"size": 7},
            linewidths=0.3, ax=ax,
            cbar_kws={"label": (f"log2({COND_SHORT[CONDITION_CNTRL]}/"
                                f"{COND_SHORT[CONDITION_TEST]}) "
                                f"within-niche proportion"),
                      "shrink": 0.6})
ax.set_yticklabels([f"N{n}" for n in niche_order], rotation=0)
color_yticklabels(ax, niche_order, niche_colors)
ax.set_xticklabels(
    [f"{ct}{(' ' + pval_stars(ct_sig.get(ct, np.nan))).rstrip()}"
     for ct in all_cts], rotation=45, ha="right")
ax.set_title(
    f"Differential niche composition: "
    f"{COND_SHORT[CONDITION_CNTRL]} vs {COND_SHORT[CONDITION_TEST]}\n"
    f"cell stars = exact MW (FDR/niche); x-label stars = mixed-model "
    f"global shift (FDR)\n{N_NOTE} — interpret as exploratory",
    fontsize=8)
ax.set_xlabel("Cell type")
ax.set_ylabel("Niche")
plt.savefig(os.path.join(OUTPUT_DIR, "differential.pdf"))
plt.close()

# ── 8.3 Per-sample niche proportions (each dot = one sample) ─────────────────
print("[PLOT] per-sample niche proportions...")
ncol = min(4, len(all_niches))
nrow = int(np.ceil(len(all_niches) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2.5 * nrow), squeeze=False)
long_ns = (niche_by_sample.reset_index()
           .melt(id_vars=NICHE_KEY, var_name="sample", value_name="prop"))
long_ns["condition"] = long_ns["sample"].map(sample_cond_map)
rng = np.random.default_rng(0)
for i, niche in enumerate(niche_order):
    ax = axes[i // ncol][i % ncol]
    nd = long_ns[long_ns[NICHE_KEY] == niche]
    for ci, cond in enumerate(CONDITION_ORDER):
        cd = nd[nd["condition"] == cond]
        ax.scatter(rng.uniform(-0.13, 0.13, len(cd)) + ci, cd["prop"],
                   c=COND_COLORS[cond], s=42,
                   edgecolors="black", linewidths=0.5, zorder=3)
        if len(cd):
            ax.plot([ci - 0.2, ci + 0.2], [cd["prop"].mean()] * 2,
                    color=COND_COLORS[cond], lw=2, zorder=2)
    row = enrichment_df[enrichment_df["niche"] == niche]
    p    = row["mw_p"].iloc[0] if len(row) else np.nan
    star = pval_stars(p)
    ttl  = (f"N{niche}"
            + (f"  {star}" if star
               else (f"  p={p:.2f}" if np.isfinite(p) else "")))
    ax.set_title(ttl, fontsize=8, color=niche_colors.get(niche, "black"))
    ax.set_xticks([0, 1])
    ax.set_xticklabels([COND_SHORT[CONDITION_CNTRL], COND_SHORT[CONDITION_TEST]])
    ax.set_ylabel("Proportion", fontsize=7)
    ax.margins(y=0.15)
for j in range(len(all_niches), nrow * ncol):
    axes[j // ncol][j % ncol].axis("off")
fig.suptitle(
    f"Per-sample niche proportions (each dot = one slide)  ·  "
    f"exact MW  ·  {N_NOTE}",
    fontsize=10, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_proportions_per_sample.pdf"))
plt.close()

# =============================================================================
# 9. CONDITION SCORE FIGURES
# =============================================================================
def split_box_strip(ax, df, group_col, group_order, color_map,
                    sig_map=None, value_col="score"):
    """Per-group Ctrl/Test split: boxplot + jitter + significance bracket."""
    rng = np.random.default_rng(0)
    for gi, g in enumerate(group_order):
        gd = df[df[group_col] == g]
        for ci, cond in enumerate(CONDITION_ORDER):
            cd = gd[gd["condition"] == cond][value_col].values
            if len(cd) == 0:
                continue
            xc = gi + (ci - 0.5) * 0.34
            bp = ax.boxplot([cd], positions=[xc], widths=0.28,
                            showfliers=False, patch_artist=True,
                            manage_ticks=False)
            for box in bp["boxes"]:
                box.set(facecolor=COND_COLORS[cond], alpha=0.55,
                        edgecolor="black", linewidth=0.6)
            for med in bp["medians"]:
                med.set(color="black", linewidth=1.0)
        # Significance bracket from sig_map (mixed-model padj)
        if sig_map is not None:
            star = pval_stars(sig_map.get(g, np.nan))
            if star:
                ytop = gd[value_col].quantile(0.97)
                ax.plot([gi - 0.17, gi + 0.17], [ytop, ytop],
                        color="black", lw=0.8)
                ax.text(gi, ytop, star, ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(group_order)))


if HAS_SCORE:
    # ── 9.1 Condition score by cell type ──────────────────────────────────────
    print("[PLOT] condition score by cell type...")
    order_ct = score_ct["celltype"].tolist()
    cdf = adata.obs[[CONDITION_SCORE_KEY, CELLTYPE_KEY, CONDITION_KEY]].copy()
    cdf.columns = ["score", "celltype", "condition"]
    fig, ax = plt.subplots(figsize=(max(7, len(order_ct) * 0.7), 4.2))
    split_box_strip(ax, cdf, "celltype", order_ct, ct_colors,
                    sig_map={r.celltype: r.padj
                             for r in score_ct.itertuples()})
    ax.set_xticklabels(order_ct, rotation=45, ha="right")
    for tick, ct in zip(ax.get_xticklabels(), order_ct):
        tick.set_color(ct_colors.get(ct, "black"))
    ax.set_ylabel("Condition score")
    ax.set_title(
        f"Condition score by cell type "
        f"({COND_SHORT[CONDITION_CNTRL]} vs {COND_SHORT[CONDITION_TEST]})\n"
        f"stars = mixed-model condition effect, FDR; {N_NOTE}",
        fontsize=9)
    _legend_conditions(fig)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "condition_score_by_celltype.pdf"))
    plt.close()

    # ── 9.2 Condition score by niche ───────────────────────────────────────────
    print("[PLOT] condition score by niche...")
    order_ni = (score_niche.sort_values("mean_score", ascending=False)
                ["niche"].tolist())
    cdf = adata.obs[[CONDITION_SCORE_KEY, NICHE_KEY, CONDITION_KEY]].copy()
    cdf.columns = ["score", "niche", "condition"]
    fig, ax = plt.subplots(figsize=(max(7, len(order_ni) * 0.7), 4.2))
    split_box_strip(ax, cdf, "niche", order_ni, niche_colors,
                    sig_map={r.niche: r.padj
                             for r in score_niche.itertuples()})
    ax.set_xticklabels([f"N{n}" for n in order_ni])
    for tick, n in zip(ax.get_xticklabels(), order_ni):
        tick.set_color(niche_colors.get(n, "black"))
    ax.set_ylabel("Condition score")
    ax.set_title(
        f"Condition score by niche "
        f"({COND_SHORT[CONDITION_CNTRL]} vs {COND_SHORT[CONDITION_TEST]})\n"
        f"stars = mixed-model condition effect, FDR; {N_NOTE}",
        fontsize=9)
    _legend_conditions(fig)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "condition_score_by_niche.pdf"))
    plt.close()

    # ── 9.3 Is the niche abundance shift coupled to the condition score? ───────
    print("[PLOT] condition score vs compositional shift...")
    merged = score_niche.merge(
        enrichment_df[["niche", "log2FC_Ctrl_vs_Test", "mw_padj"]],
        on="niche", how="left")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    # Left: mean niche score vs niche abundance log2FC
    ax = axes[0]
    for r in merged.itertuples():
        ax.scatter(r.mean_score, r.log2FC_Ctrl_vs_Test, s=90,
                   color=niche_colors.get(r.niche, "gray"),
                   edgecolors="black", linewidths=0.5, zorder=3)
        ax.annotate(f"N{r.niche}", (r.mean_score, r.log2FC_Ctrl_vs_Test),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, ls="--", color="gray", lw=0.6)
    valid = merged.dropna(subset=["mean_score", "log2FC_Ctrl_vs_Test"])
    if len(valid) >= 3:
        rho, pv = spearmanr(valid["mean_score"], valid["log2FC_Ctrl_vs_Test"])
        ax.set_title(
            f"Niche condition score vs abundance shift\n"
            f"Spearman ρ={rho:.2f}, p={pv:.2f}",
            fontsize=9)
    ax.set_xlabel("Mean niche condition score")
    ax.set_ylabel(f"log2({COND_SHORT[CONDITION_CNTRL]}/{COND_SHORT[CONDITION_TEST]}) "
                  f"niche abundance")

    # Right: per-niche Ctrl-Test condition score effect vs abundance log2FC
    ax = axes[1]
    for r in merged.itertuples():
        eff = getattr(r, "effect_Ctrl_vs_Test", np.nan)
        if np.isnan(eff):
            continue
        sig = (not np.isnan(getattr(r, "padj", np.nan))) and r.padj < 0.05
        ax.scatter(eff, r.log2FC_Ctrl_vs_Test,
                   s=110 if sig else 70,
                   color=niche_colors.get(r.niche, "gray"),
                   edgecolors="black" if sig else "gray",
                   linewidths=1.3 if sig else 0.5, zorder=3)
        ax.annotate(f"N{r.niche}" + ("*" if sig else ""),
                    (eff, r.log2FC_Ctrl_vs_Test), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, ls="--", color="gray", lw=0.6)
    ax.axvline(0, ls="--", color="gray", lw=0.6)
    ax.set_xlabel(f"Condition score {COND_SHORT[CONDITION_CNTRL]}−"
                  f"{COND_SHORT[CONDITION_TEST]} effect (mixed model)")
    ax.set_ylabel(f"log2({COND_SHORT[CONDITION_CNTRL]}/{COND_SHORT[CONDITION_TEST]}) "
                  f"niche abundance")
    ax.set_title(
        "Is the abundance shift coupled to a\ncondition score shift? "
        "(* niche score FDR<0.05)",
        fontsize=9)
    fig.suptitle(
        f"Niche redistribution vs condition score  ·  {N_NOTE}",
        fontsize=11, y=1.03)
    merged.to_csv(
        os.path.join(OUTPUT_DIR, "condition_score_vs_composition_shift.csv"),
        index=False)
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUTPUT_DIR, "condition_score_vs_composition_shift.pdf"))
    plt.close()

    # ── 9.4 Condition-score-associated gene programs ───────────────────────────
    if gp_score is not None and len(gp_score) > 0:
        print("[PLOT] condition-score-associated gene programs...")
        top = (pd.concat([gp_score.head(15), gp_score.tail(15)])
               .drop_duplicates("gene_program")
               .sort_values("spearman_rho"))
        fig, ax = plt.subplots(figsize=(7, max(4, len(top) * 0.28)))
        colors = ["#C94A4A" if v > 0 else "#4A7EB8"
                  for v in top["spearman_rho"]]
        ax.barh(range(len(top)), top["spearman_rho"],
                color=colors, edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([gp[:48] for gp in top["gene_program"]], fontsize=6)
        for yi, padj in enumerate(top["padj"]):
            star = pval_stars(padj)
            if star:
                v = top["spearman_rho"].iloc[yi]
                ax.text(v + (0.01 if v >= 0 else -0.01), yi, star,
                        va="center",
                        ha="left" if v >= 0 else "right", fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("Spearman ρ (GP activity vs condition score)")
        ax.set_title(
            "Condition-score-associated gene programs\n"
            "(per-cell ρ; stars FDR; exploratory)",
            fontsize=9)
        plt.tight_layout()
        plt.savefig(
            os.path.join(OUTPUT_DIR, "condition_score_gp_association.pdf"))
        plt.close()

# =============================================================================
# 10. DONE
# =============================================================================
print(f"\n[DONE] Step 2.5a complete. Outputs in {OUTPUT_DIR}/")
print(f"  Statistics: niche abundance & within-niche composition use the SLIDE")
print(f"  as the unit (exact Mann-Whitney). Cell-type compositional shift and all")
print(f"  condition-score comparisons use linear mixed models with sample as a")
print(f"  random intercept. {N_NOTE}; treat per-pair proportion p-values as")
print(f"  exploratory.")
if HAS_SCORE:
    print(f"  Condition score: {len(score_genes_used)} signature genes used "
          f"(see {SHARED_DIR}/condition_genes_used.txt).")