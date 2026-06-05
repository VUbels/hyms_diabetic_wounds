"""
==============================================================================
SPATIAL TRANSCRIPTOMICS PIPELINE — STEP 2.5b  (revised)
Niche spatial architecture: proximity, contact, neighbourhood enrichment
==============================================================================

What changed vs the previous version
------------------------------------
  * Shared niche colour map (same niche = same colour as in 2.5a / 2.5c) and a
    shared Kamada-Kawai layout, so the Ctrl and Test networks are directly
    comparable panel-to-panel and node-to-node.
  * POOLED edge thresholds (computed once over both conditions) instead of
    per-condition quantiles — an edge drawn "strong" means the same physical
    proximity in both panels.
  * Neighbourhood enrichment (squidpy) added as the field-standard niche-
    adjacency metric: per-slide z-scores aggregated to a per-condition mean,
    with Ctrl / Test / Δ heatmaps. Proximity (closest-approach distance) is kept as
    a supplement and is explicitly flagged as density-sensitive; contact
    frequency (spatial intermingling) is promoted alongside it.
  * Mantel interpretation corrected: a SIGNIFICANT Mantel (high r, low p) means
    the architecture is CONSERVED across conditions, not reorganised.
  * A slide-respecting global permutation test of architectural reorganisation
    (PERMANOVA-style pseudo-F on per-slide edge vectors, exact label
    permutation) — the honest global test given the slide is the unit.
  * Optional condition-score-coloured network: nodes shaded by mean condition score score on
    the shared layout, to reveal high-condition-score niche hubs in their spatial context.
  * Significance markers / captions state the test and the n explicitly.

Sample size note: with n=3 HS vs n=2 LS slides, the minimum achievable exact
permutation / Mann-Whitney p for a single edge is 1/10 = 0.10, so no single
edge survives multiple-testing correction. Per-edge results are therefore
reported as an effect-size atlas; FDR-controlled inference is reserved for the
pre-specified FOCUS_EDGES set and for the single global architecture test.

Inputs
------
    nichecompass_results/objects/nichecompass_integrated.h5ad   (or $NICHE_H5AD)

Outputs (niche_network/)
------------------------
    niche_network_{HS,LS,combined}.pdf          per-condition + overlaid (shared layout)
    niche_architecture_change.pdf               prevalence shift vs neighbour rewiring
    niche_egocentric_change.pdf                 LS-vs-HS neighbours for most-changed niches
    niche_network_condition_score.pdf
    niche_network_differential_{mantel,heatmap,top_edges}.pdf
    niche_nhood_enrichment_heatmap.pdf
    niche_proximity_heatmap.pdf / niche_contact_frequency_heatmap.pdf
    niche_network_focus_edges.pdf               (only if FOCUS_EDGES set)
    global_architecture_test.txt / mantel_test_result.txt
    niche_architecture_change.csv               (per-niche prevalence + rewiring table)
    niche_proximity_{per_slide,per_condition,differential}.csv
    niche_nhood_enrichment_{HS,LS,delta}.csv
==============================================================================
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from math import comb
from itertools import combinations
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except Exception:
    HAS_ADJUSTTEXT = False

try:
    import squidpy as sq
    HAS_SQUIDPY = True
except Exception:
    HAS_SQUIDPY = False

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
    "--output_dir", default=None, metavar="DIR",
    help="Output directory (default: niche_network)")
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
INPUT_H5AD  = os.environ.get(
    "NICHE_H5AD", "nichecompass_results/objects/nichecompass_integrated.h5ad")
OUTPUT_DIR  = "niche_network"

# CLI overrides (parser runs before this block)
if _args.h5ad:           INPUT_H5AD  = _args.h5ad
if _args.output_dir:     OUTPUT_DIR  = _args.output_dir
if _args.exclude_samples is not None: EXCLUDE_SAMPLES = _args.exclude_samples

SPATIAL_KEY   = "X_spatial"
NICHE_KEY     = "nichecompass_niche"
CELLTYPE_KEY  = "predicted_cell_type"
SAMPLE_KEY    = "sample_id"
CONDITION_KEY = "condition"
LATENT_KEY    = "nichecompass_latent"
ACTIVE_GP_KEY = "nichecompass_active_gp_names"

CONDITIONS = list(CONDITION_ORDER)

EXCLUDE_SAMPLES = [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]

PROXIMITY_PERCENTILE = 25
CONTACT_K = 8                     # k for contact frequency / neighbourhood graph

# -- Niche filtering --
MIN_CELLS_PER_NICHE_SLIDE = 30
MIN_SLIDES_WITH_NICHE     = 2
MIN_TOTAL_FRACTION        = 0.005

# -- Edge thresholding (POOLED across conditions -> comparable panels) --
STRONG_EDGE_Q = 0.20
WEAK_EDGE_Q   = 0.40
STRONG_EDGE_COLOR = "#1F3D8A"; STRONG_EDGE_WIDTH = 3.5; STRONG_EDGE_ALPHA = 0.95
WEAK_EDGE_COLOR   = "#A8C2F0"; WEAK_EDGE_WIDTH   = 1.2; WEAK_EDGE_ALPHA   = 0.55

# -- Differential views --
N_TOP_DIFF_EDGES = 15
N_PERMUTATIONS   = 9999          # cap; exact permutation used when feasible

# -- Hypothesis-driven follow-up (figure D; activates when non-empty) --
FOCUS_EDGES = []
FOCUS_FDR_THRESHOLD = 0.1

# -- Layout / nodes --
KK_SCALE       = 4.0
NODE_SIZE_MIN  = 600
NODE_SIZE_MAX  = 2200
LABEL_FONTSIZE = 9

# Optional manual niche annotation: {"3": {"name": "...", "color": "#...",
# "group": "..."}}. Colour overrides the shared map; name/group are cosmetic.
NICHE_ANNOTATION = {}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# 1. LOAD
# =============================================================================
print(f"[LOAD] {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
adata.obs[NICHE_KEY] = adata.obs[NICHE_KEY].astype(str)
adata.obs[CELLTYPE_KEY] = adata.obs[CELLTYPE_KEY].astype(str)
adata.obs[CONDITION_KEY] = adata.obs[CONDITION_KEY].astype(str)
adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype(str)

if EXCLUDE_SAMPLES:
    keep = ~adata.obs[SAMPLE_KEY].isin(EXCLUDE_SAMPLES)
    print(f"  Excluding {EXCLUDE_SAMPLES}: dropping {(~keep).sum()} cells")
    adata = adata[keep].copy()

all_niches_raw = sorted(adata.obs[NICHE_KEY].unique(), key=niche_sort_key)
print(f"  {adata.n_obs} cells, {len(all_niches_raw)} niches (pre-filter)")

# Condition score (for the condition-score-coloured network) + shared styles.
# Styles are built on the FULL niche set so colours match the other scripts.
print("\n[SCORE] Scoring condition score condition score signature...")
compute_condition_score(adata)
score_thr = add_condition_flag(adata)
HAS_SCORE = (CONDITION_SCORE_KEY in adata.obs
           and not adata.obs[CONDITION_SCORE_KEY].isna().all())

niche_colors, _ = build_or_load_niche_style(adata, NICHE_KEY, CELLTYPE_KEY)

# =============================================================================
# 2. NICHE FILTERING
# =============================================================================
print("\n[FILTER] Niche cell-count thresholds:")
print(f"  >={MIN_CELLS_PER_NICHE_SLIDE} cells/slide in "
      f">={MIN_SLIDES_WITH_NICHE}/{adata.obs[SAMPLE_KEY].nunique()} slides; "
      f">={MIN_TOTAL_FRACTION*100:.1f}% of total cells")

niche_counts_per_slide = pd.crosstab(adata.obs[NICHE_KEY], adata.obs[SAMPLE_KEY]) \
    .reindex(all_niches_raw, fill_value=0)
total_per_niche = niche_counts_per_slide.sum(axis=1)
total_cells = adata.n_obs

filter_log, retained_niches = [], []
for niche in all_niches_raw:
    cps = niche_counts_per_slide.loc[niche]
    slides_above = int((cps >= MIN_CELLS_PER_NICHE_SLIDE).sum())
    frac = float(total_per_niche.loc[niche] / total_cells)
    kept = (slides_above >= MIN_SLIDES_WITH_NICHE) and (frac >= MIN_TOTAL_FRACTION)
    filter_log.append({"niche": niche, "total_cells": int(total_per_niche.loc[niche]),
                       "total_fraction": round(frac, 4), "slides_above": slides_above,
                       "max_cells_any_slide": int(cps.max()), "retained": bool(kept)})
    if kept:
        retained_niches.append(niche)
pd.DataFrame(filter_log).to_csv(os.path.join(OUTPUT_DIR, "niche_filtering_log.csv"),
                                index=False)
all_niches = retained_niches
print(f"  Retained {len(all_niches)}/{len(all_niches_raw)}: {all_niches}")
if len(all_niches) < 2:
    raise RuntimeError("Fewer than 2 niches survived filtering — relax thresholds.")
adata = adata[adata.obs[NICHE_KEY].isin(all_niches)].copy()

annotated = bool(NICHE_ANNOTATION)


def get_niche_display(niche_id):
    """(name, color, group). Colour from the shared similarity map unless
    NICHE_ANNOTATION overrides it; name/group are cosmetic overrides."""
    n = str(niche_id)
    name, group = f"N{n}", "other"
    color = niche_colors.get(n, "#9E9E9E")
    if annotated and n in NICHE_ANNOTATION:
        ann = NICHE_ANNOTATION[n]
        name = ann.get("name", name)
        color = ann.get("color", color)
        group = ann.get("group", group)
    return name, color, group


display_labels = [get_niche_display(n)[0] for n in all_niches]

# Node sizes scaled by niche abundance (shared across all network panels)
_frac = (adata.obs[NICHE_KEY].value_counts() / adata.n_obs)
_fmax = _frac.reindex(all_niches).max()
niche_sizes = {n: NODE_SIZE_MIN + (NODE_SIZE_MAX - NODE_SIZE_MIN)
               * float(np.sqrt(max(_frac.get(n, 0.0), 0.0) / (_fmax + 1e-12)))
               for n in all_niches}

slide_cond = adata.obs.groupby(SAMPLE_KEY)[CONDITION_KEY].first().to_dict()
slides_list = sorted(slide_cond.keys())
n_slides = len(slides_list)
n_ctrl = sum(1 for c in slide_cond.values() if c == CONDITION_CNTRL)
n_test = n_slides - n_ctrl
N_NOTE = f"n={n_ctrl} HS vs n={n_test} LS slides"

# =============================================================================
# 3. PER-SLIDE PROXIMITY + CONTACT FREQUENCY (vectorised)
# =============================================================================
def niche_proximity_per_slide(adata, spatial_key, niche_key, sample_key,
                              min_cells, percentile, contact_k=CONTACT_K):
    records = []
    for slide in sorted(adata.obs[sample_key].unique()):
        mask = (adata.obs[sample_key] == slide).values
        coords = np.asarray(adata.obsm[spatial_key])[mask]
        niches = adata.obs[niche_key].values[mask].astype(str)
        if len(coords) < 3:
            continue
        cond = adata.obs[CONDITION_KEY].values[mask][0]
        present = sorted([n for n in np.unique(niches)
                          if (niches == n).sum() >= min_cells], key=niche_sort_key)
        if len(present) < 2:
            continue
        code = {n: i for i, n in enumerate(present)}
        cell_code = np.array([code.get(n, -1) for n in niches])

        tree = cKDTree(coords)
        k = min(contact_k + 1, len(coords))
        _, nbrs = tree.query(coords, k=k)
        if nbrs.ndim == 1:
            nbrs = nbrs[:, None]
        nbrs = nbrs[:, 1:]
        neigh_code = cell_code[nbrs]                       # (n_cells, k)
        has = np.zeros((len(coords), len(present)), dtype=bool)
        for j in range(len(present)):
            has[:, j] = np.any(neigh_code == j, axis=1)

        trees = {n: cKDTree(coords[niches == n]) for n in present}
        ncoords = {n: coords[niches == n] for n in present}
        for i, nA in enumerate(present):
            aE = cell_code == i
            for nB in present[i + 1:]:
                j = code[nB]
                bE = cell_code == j
                cab = has[aE, j].mean() if aE.sum() else np.nan
                cba = has[bE, i].mean() if bE.sum() else np.nan
                contact = float(np.nanmean([cab, cba]))
                dAB, _ = trees[nB].query(ncoords[nA], k=1)
                dBA, _ = trees[nA].query(ncoords[nB], k=1)
                prox = float(np.percentile(np.concatenate([dAB, dBA]), percentile))
                records.append({"slide": slide, "condition": cond,
                                "niche_A": nA, "niche_B": nB,
                                "proximity_p25": prox, "contact_freq": contact,
                                "n_A": int(ncoords[nA].shape[0]),
                                "n_B": int(ncoords[nB].shape[0])})
    return pd.DataFrame(records)


print(f"\n[PROXIMITY] p{PROXIMITY_PERCENTILE} NN distance + contact freq per slide...")
prox_per_slide = niche_proximity_per_slide(
    adata, SPATIAL_KEY, NICHE_KEY, SAMPLE_KEY,
    min_cells=MIN_CELLS_PER_NICHE_SLIDE, percentile=PROXIMITY_PERCENTILE)
prox_per_slide.to_csv(os.path.join(OUTPUT_DIR, "niche_proximity_per_slide.csv"),
                      index=False)
print(f"  {len(prox_per_slide)} niche-pair-slide records")

prox_by_cond = (prox_per_slide.groupby(["condition", "niche_A", "niche_B"])
                .agg(mean_proximity=("proximity_p25", "mean"),
                     sd_proximity=("proximity_p25", "std"),
                     mean_contact_freq=("contact_freq", "mean"),
                     sd_contact_freq=("contact_freq", "std"),
                     n_slides=("proximity_p25", "count")).reset_index())
prox_by_cond.to_csv(os.path.join(OUTPUT_DIR, "niche_proximity_per_condition.csv"),
                    index=False)

# =============================================================================
# 4. PER-EDGE DIFFERENTIAL (exact MW + slide-label permutation)
# =============================================================================
print("\n[DIFF] Δ proximity / Δ contact, exact MW + permutation per edge...")
n_arr = comb(n_slides, n_ctrl)
use_exact = n_arr <= N_PERMUTATIONS
if use_exact:
    all_hs_assignments = list(combinations(range(n_slides), n_ctrl))
    print(f"  exact permutation: {len(all_hs_assignments)} arrangements")
else:
    rng = np.random.default_rng(42)
    all_hs_assignments = [tuple(sorted(rng.choice(n_slides, n_ctrl, replace=False)))
                          for _ in range(N_PERMUTATIONS)]
    print(f"  sampled permutation: {N_PERMUTATIONS} arrangements")


def _perm_p(slide_to_val, observed_abs):
    n_ge = tot = 0
    for combo in all_hs_assignments:
        hs_s = [slides_list[i] for i in combo]
        hv = [slide_to_val[s] for s in hs_s if s in slide_to_val]
        lv = [slide_to_val[s] for s in slides_list
              if s not in hs_s and s in slide_to_val]
        if not hv or not lv:
            continue
        if abs(np.mean(hv) - np.mean(lv)) >= observed_abs - 1e-9:
            n_ge += 1
        tot += 1
    return n_ge / tot if tot else np.nan


def compute_edge_stats(pair_df):
    hv = pair_df[pair_df.condition == CONDITION_CNTRL]["proximity_p25"].values
    lv = pair_df[pair_df.condition == CONDITION_TEST]["proximity_p25"].values
    hc = pair_df[pair_df.condition == CONDITION_CNTRL]["contact_freq"].values
    lc = pair_df[pair_df.condition == CONDITION_TEST]["contact_freq"].values
    out = {"Ctrl_mean": hv.mean() if len(hv) else np.nan,
           "Test_mean": lv.mean() if len(lv) else np.nan,
           "n_Ctrl": len(hv), "n_Test": len(lv),
           "delta_proximity": np.nan,
           "Ctrl_contact_freq": hc.mean() if len(hc) else np.nan,
           "Test_contact_freq": lc.mean() if len(lc) else np.nan,
           "delta_contact_freq": np.nan,
           "mw_p": np.nan, "perm_p": np.nan, "contact_perm_p": np.nan}
    if len(hv) == 0 or len(lv) == 0:
        return out
    out["delta_proximity"] = hv.mean() - lv.mean()
    if len(hc) and len(lc):
        out["delta_contact_freq"] = hc.mean() - lc.mean()
    if len(hv) >= 2 and len(lv) >= 2:
        try:
            out["mw_p"] = mannwhitneyu(hv, lv, alternative="two-sided",
                                       method="exact")[1]
        except Exception:
            pass
    out["perm_p"] = _perm_p(dict(zip(pair_df.slide, pair_df.proximity_p25)),
                            abs(out["delta_proximity"]))
    if not np.isnan(out["delta_contact_freq"]):
        out["contact_perm_p"] = _perm_p(dict(zip(pair_df.slide, pair_df.contact_freq)),
                                        abs(out["delta_contact_freq"]))
    return out


diff_rows = []
for (nA, nB), sub in prox_per_slide.groupby(["niche_A", "niche_B"]):
    st = compute_edge_stats(sub)
    st.update({"niche_A": nA, "niche_B": nB})
    diff_rows.append(st)
diff_df = pd.DataFrame(diff_rows)
for pc in ["mw_p", "perm_p", "contact_perm_p"]:
    v = diff_df[pc].notna()
    diff_df[pc.replace("_p", "_padj")] = np.nan
    if v.sum() > 0:
        diff_df.loc[v, pc.replace("_p", "_padj")] = multipletests(
            diff_df.loc[v, pc], method="fdr_bh")[1]
diff_df = diff_df.sort_values("delta_proximity", key=lambda x: x.abs(),
                              ascending=False).reset_index(drop=True)
diff_df.to_csv(os.path.join(OUTPUT_DIR, "niche_proximity_differential.csv"),
               index=False)

# =============================================================================
# 5. GLOBAL ARCHITECTURE TESTS
#    (a) Mantel concordance (high r, low p => CONSERVED)
#    (b) slide-respecting permutation PERMANOVA-style pseudo-F (reorganisation)
# =============================================================================
def build_matrix(prox_subset, niches, col):
    mat = pd.DataFrame(np.nan, index=niches, columns=niches)
    for _, r in prox_subset.iterrows():
        if r["niche_A"] in niches and r["niche_B"] in niches:
            mat.at[r["niche_A"], r["niche_B"]] = r[col]
            mat.at[r["niche_B"], r["niche_A"]] = r[col]
    return mat


hs_mat = build_matrix(prox_by_cond[prox_by_cond.condition == CONDITION_CNTRL],
                      all_niches, "mean_proximity")
ls_mat = build_matrix(prox_by_cond[prox_by_cond.condition == CONDITION_TEST],
                      all_niches, "mean_proximity")
np.fill_diagonal(hs_mat.values, 0.0)
np.fill_diagonal(ls_mat.values, 0.0)


def mantel_test(mat_A, mat_B, n_perm=9999, seed=42):
    n = mat_A.shape[0]
    iu = np.triu_indices(n, k=1)
    a, b = mat_A.values[iu], mat_B.values[iu]
    valid = ~(np.isnan(a) | np.isnan(b))
    a, b = a[valid], b[valid]
    if len(a) < 3:
        return np.nan, np.nan, len(a)
    obs = pearsonr(a, b)[0]
    rng = np.random.default_rng(seed)
    iu_full = np.triu_indices(n, k=1)
    vmask = valid
    n_ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        ap = mat_A.values[np.ix_(perm, perm)][iu_full][vmask]
        try:
            if abs(pearsonr(ap, b)[0]) >= abs(obs) - 1e-9:
                n_ge += 1
        except Exception:
            continue
    return obs, (n_ge + 1) / (n_perm + 1), len(a)


print("\n[MANTEL] Ctrl vs Test proximity-matrix concordance...")
mantel_r, mantel_p, n_pairs = mantel_test(hs_mat, ls_mat, n_perm=N_PERMUTATIONS)
print(f"  Mantel r = {mantel_r:.3f}, p = {mantel_p:.4f} (n_pairs={n_pairs})")
with open(os.path.join(OUTPUT_DIR, "mantel_test_result.txt"), "w") as f:
    f.write("Mantel test: Ctrl vs Test niche proximity matrices\n")
    f.write(f"  n niche pairs : {n_pairs}\n  Pearson r     : {mantel_r:.4f}\n")
    f.write(f"  p-value       : {mantel_p:.4f}\n  permutations  : {N_PERMUTATIONS}\n\n")
    f.write("Interpretation (corrected):\n")
    f.write("  The Mantel test asks whether the two proximity matrices are\n")
    f.write("  CORRELATED. High r with a SIGNIFICANT (low) p therefore means the\n")
    f.write("  niche architecture is CONSERVED between conditions, NOT reorganised.\n")
    f.write("  Reorganisation is indicated by (i) departure of points from the\n")
    f.write("  y=x line in the scatter, (ii) large per-edge Δ values, and (iii)\n")
    f.write("  the global permutation pseudo-F test below.\n")


def global_architecture_permutation(per_slide_df, metric_col):
    piv = per_slide_df.pivot_table(index="slide", columns=["niche_A", "niche_B"],
                                   values=metric_col)
    slides = list(piv.index)
    if len(slides) < 4:
        return None
    full = piv.dropna(axis=1, how="any")
    if full.shape[1] < 3:
        keep = piv.columns[piv.notna().sum(axis=0) >= (len(slides) - 1)]
        full = piv[keep].apply(lambda c: c.fillna(c.mean()), axis=0)
        full = full.dropna(axis=1, how="any")
    labels = np.array([slide_cond[s] for s in slides])
    if full.shape[1] < 3 or len(set(labels)) < 2:
        return None
    X = full.values

    def pseudo_f(lab):
        g = X.mean(axis=0)
        sst = ((X - g) ** 2).sum()
        ssw = sum(((X[lab == k] - X[lab == k].mean(axis=0)) ** 2).sum()
                  for k in np.unique(lab))
        dfb, dfw = len(np.unique(lab)) - 1, len(lab) - len(np.unique(lab))
        return ((sst - ssw) / dfb) / (ssw / dfw) if (ssw > 0 and dfw > 0) else np.nan

    def cdist(lab):
        gs = np.unique(lab)
        return (np.linalg.norm(X[lab == gs[0]].mean(axis=0)
                               - X[lab == gs[1]].mean(axis=0))
                if len(gs) == 2 else np.nan)

    of, od = pseudo_f(labels), cdist(labels)
    nhsi = int((labels == CONDITION_CNTRL).sum())
    perms = list(combinations(range(len(slides)), nhsi))
    cf = cd = tot = 0
    for combo in perms:
        lab = np.array([CONDITION_TEST] * len(slides), dtype=object)
        for i in combo:
            lab[i] = CONDITION_CNTRL
        f, d = pseudo_f(lab), cdist(lab)
        if not np.isnan(f) and f >= of - 1e-9:
            cf += 1
        if not np.isnan(d) and d >= od - 1e-9:
            cd += 1
        tot += 1
    return {"metric": metric_col, "n_edges": int(full.shape[1]), "n_slides": len(slides),
            "pseudo_F": of, "perm_p_F": cf / tot if tot else np.nan,
            "centroid_dist": od, "perm_p_dist": cd / tot if tot else np.nan,
            "n_perm": tot}


print("[GLOBAL] Slide-respecting permutation test of reorganisation...")
glob_results = []
for col in ["proximity_p25", "contact_freq"]:
    r = global_architecture_permutation(prox_per_slide, col)
    if r:
        glob_results.append(r)
        print(f"  {col}: pseudo-F={r['pseudo_F']:.2f}, perm p={r['perm_p_F']:.3f} "
              f"({r['n_edges']} common edges, {r['n_perm']} perms)")
with open(os.path.join(OUTPUT_DIR, "global_architecture_test.txt"), "w") as f:
    f.write(f"Global niche-architecture reorganisation test ({N_NOTE})\n")
    f.write("PERMANOVA-style pseudo-F on per-slide edge vectors; exact label\n")
    f.write("permutation (slide = unit). Larger F / smaller p => Ctrl and Test slides\n")
    f.write("differ more in architecture than expected by chance.\n\n")
    for r in glob_results:
        f.write(f"[{r['metric']}] edges={r['n_edges']}, slides={r['n_slides']}, "
                f"pseudo-F={r['pseudo_F']:.3f}, perm p(F)={r['perm_p_F']:.3f}, "
                f"centroid-dist p={r['perm_p_dist']:.3f}, n_perm={r['n_perm']}\n")
    f.write(f"\nNOTE: with {N_NOTE} the smallest achievable permutation p is "
            f"{1.0/max(comb(n_slides, n_ctrl),1):.2f}.\n")
if glob_results:
    pd.DataFrame(glob_results).to_csv(
        os.path.join(OUTPUT_DIR, "global_architecture_test.csv"), index=False)

# Mantel scatter (corrected caption)
fig, ax = plt.subplots(figsize=(6.5, 6.5))
iu = np.triu_indices(len(all_niches), k=1)
vh, vl = hs_mat.values[iu], ls_mat.values[iu]
ok = ~(np.isnan(vh) | np.isnan(vl))
ax.scatter(vl[ok], vh[ok], alpha=0.7, s=45, edgecolors="black",
           linewidths=0.3, color="#6A8CBF")
if ok.sum():
    lo = float(min(vl[ok].min(), vh[ok].min()))
    hi = float(max(vl[ok].max(), vh[ok].max()))
    ax.plot([lo, hi], [lo, hi], "--", color="gray", alpha=0.6, label="Conserved (y=x)")
ax.set_xlabel("LS proximity (µm)"); ax.set_ylabel("HS proximity (µm)")
gp = glob_results[0]["perm_p_F"] if glob_results else np.nan
ax.set_title(f"Global niche architecture: Ctrl vs Test\nMantel r={mantel_r:.3f}, "
             f"p={mantel_p:.3f} (high r + low p ⇒ conserved)\n"
             f"reorganisation perm pseudo-F p={gp:.3f}", fontsize=9)
ax.legend(loc="upper left"); ax.set_aspect("equal")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_network_differential_mantel.pdf"))
plt.close()

# =============================================================================
# 6. NEIGHBOURHOOD ENRICHMENT (squidpy) per slide -> per condition
# =============================================================================
nhood = {}
if HAS_SQUIDPY:
    print("\n[NHOOD] squidpy neighbourhood enrichment per slide...")

    def nhood_per_condition(adata, niche_key, spatial_key, niches, contact_k):
        per = {c: [] for c in CONDITIONS}
        for slide in sorted(adata.obs[SAMPLE_KEY].unique()):
            sub = adata[adata.obs[SAMPLE_KEY] == slide].copy()
            cond = sub.obs[CONDITION_KEY].values[0]
            sub.obs[niche_key] = sub.obs[niche_key].astype("category")
            sub.obsm["spatial"] = np.asarray(sub.obsm[spatial_key])
            try:
                sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=contact_k)
                sq.gr.nhood_enrichment(sub, cluster_key=niche_key, seed=0,
                                       show_progress_bar=False)
                z = sub.uns[f"{niche_key}_nhood_enrichment"]["zscore"]
                cats = list(sub.obs[niche_key].cat.categories)
                per[cond].append(pd.DataFrame(z, index=cats, columns=cats)
                                 .reindex(index=niches, columns=niches))
            except Exception as e:
                print(f"  [NHOOD][WARN] slide {slide}: {e}")
        out = {}
        for c in CONDITIONS:
            if per[c]:
                out[c] = (pd.concat(per[c]).groupby(level=0).mean()
                          .reindex(index=niches, columns=niches))
        return out

    try:
        nhood = nhood_per_condition(adata, NICHE_KEY, SPATIAL_KEY, all_niches, CONTACT_K)
    except Exception as e:
        print(f"  [NHOOD][WARN] skipped: {e}")
        nhood = {}

    if nhood:
        for c, m in nhood.items():
            m.to_csv(os.path.join(OUTPUT_DIR, f"niche_nhood_enrichment_{COND_SHORT[c]}.csv"))
        if all(c in nhood for c in CONDITIONS):
            dz = nhood[CONDITION_CNTRL] - nhood[CONDITION_TEST]
            dz.to_csv(os.path.join(OUTPUT_DIR, "niche_nhood_enrichment_delta.csv"))
            present = [c for c in CONDITIONS if c in nhood]
            vmax = np.nanmax([np.nanmax(np.abs(nhood[c].values)) for c in present] + [1])
            fig, axes = plt.subplots(1, len(present) + 1,
                                     figsize=(5 * (len(present) + 1), 4.6))
            lbl = [f"N{n}" for n in all_niches]
            for ax, c in zip(axes[:-1], present):
                sns.heatmap(pd.DataFrame(nhood[c].values, index=lbl, columns=lbl),
                            cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax, ax=ax,
                            linewidths=0.3, cbar_kws={"label": "enrichment z"})
                ax.set_title(COND_SHORT[c]); plt.setp(ax.get_xticklabels(),
                                                      rotation=45, ha="right")
            dzmax = np.nanmax(np.abs(dz.values)) if np.isfinite(
                np.nanmax(np.abs(dz.values))) else 1
            sns.heatmap(pd.DataFrame(dz.values, index=lbl, columns=lbl),
                        cmap="PuOr_r", center=0, vmin=-dzmax, vmax=dzmax, ax=axes[-1],
                        linewidths=0.3, cbar_kws={"label": "Δ z (Ctrl−Test)"})
            axes[-1].set_title("Ctrl - Test")
            plt.setp(axes[-1].get_xticklabels(), rotation=45, ha="right")
            fig.suptitle(f"Niche neighbourhood enrichment (squidpy)  ·  "
                         f"per-slide z averaged per condition  ·  {N_NOTE}",
                         fontsize=11, y=1.03)
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, "niche_nhood_enrichment_heatmap.pdf"))
            plt.close()
else:
    print("\n[NHOOD] squidpy not available — skipping neighbourhood enrichment.")

# =============================================================================
# 7. EFFECT-SIZE HEATMAP (Δ proximity)
# =============================================================================
delta_mat = hs_mat - ls_mat
delta_display = pd.DataFrame(delta_mat.values, index=display_labels,
                            columns=display_labels)
abs_max = np.nanmax(np.abs(delta_display.values))
abs_max = abs_max if np.isfinite(abs_max) and abs_max > 0 else 1.0
fig, ax = plt.subplots(figsize=(max(7, len(all_niches) * 0.45),
                                max(6, len(all_niches) * 0.45)))
sns.heatmap(delta_display, cmap="RdBu", center=0, vmin=-abs_max, vmax=abs_max,
            cbar_kws={"label": "Δ proximity Ctrl−Test (µm)  ·  red: closer in HS"},
            linewidths=0.3, ax=ax)
color_yticklabels(ax, all_niches, niche_colors)
ax.set_title(f"Δ proximity between conditions (effect-size atlas)\n"
             f"no statistical filter; {N_NOTE}", fontsize=9)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
plt.setp(ax.get_yticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_network_differential_heatmap.pdf"))
plt.close()

# =============================================================================
# 8. LAYOUT (Kamada-Kawai, shared across conditions) + POOLED thresholds
# =============================================================================
prox_pooled = (prox_per_slide.groupby(["niche_A", "niche_B"])["proximity_p25"]
               .mean().reset_index().rename(columns={"proximity_p25": "mean_proximity"}))
pooled_vals = prox_pooled["mean_proximity"].values
POOLED_STRONG_THR = float(np.quantile(pooled_vals, STRONG_EDGE_Q)) if len(pooled_vals) else 0
POOLED_WEAK_THR = float(np.quantile(pooled_vals, WEAK_EDGE_Q)) if len(pooled_vals) else 0


def compute_kk_layout(prox_pooled_df, niches):
    if len(prox_pooled_df) == 0:
        return {n: (0, i) for i, n in enumerate(niches)}, list(niches)
    G = nx.Graph()
    G.add_nodes_from(niches)
    for _, r in prox_pooled_df.iterrows():
        if r["niche_A"] in niches and r["niche_B"] in niches \
                and r["mean_proximity"] <= POOLED_WEAK_THR:
            G.add_edge(r["niche_A"], r["niche_B"], weight=float(r["mean_proximity"]))
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    main = comps[0] if comps else set()
    disconnected = [n for n in niches if n not in main]
    pos = {}
    if len(main) >= 2:
        sg = G.subgraph(main).copy()
        try:
            pos.update(nx.kamada_kawai_layout(sg, weight="weight", scale=KK_SCALE))
        except Exception:
            pos.update(nx.spring_layout(sg, scale=KK_SCALE, seed=42, k=1.5))
    elif len(main) == 1:
        pos[list(main)[0]] = (0.0, 0.0)
    if disconnected:
        x_max = max((p[0] for p in pos.values()), default=0.0)
        y_mid = float(np.mean([p[1] for p in pos.values()])) if pos else 0.0
        ys = np.linspace(y_mid - 0.3 * KK_SCALE, y_mid + 0.3 * KK_SCALE,
                         max(len(disconnected), 2))
        for n, y in zip(disconnected, ys):
            pos[n] = (x_max + 0.35 * KK_SCALE, y)
    return pos, disconnected


shared_pos, disconnected = compute_kk_layout(prox_pooled, all_niches)
print(f"\n[LAYOUT] KK main component: {len(all_niches) - len(disconnected)} niches; "
      f"pooled thresholds strong<={POOLED_STRONG_THR:.1f} weak<={POOLED_WEAK_THR:.1f} µm")


def _draw_nodes(ax, pos, niches, disconnected, sizes=None, facecolor=None,
                edgecolor=None):
    for n in niches:
        if n not in pos:
            continue
        c = (facecolor[n] if facecolor else get_niche_display(n)[1])
        s = (sizes.get(n) if sizes else NODE_SIZE_MAX)
        disc = n in disconnected
        ax.scatter(pos[n][0], pos[n][1], s=s, c=[c],
                   edgecolors=(edgecolor or ("gray" if disc else "black")),
                   linewidths=1.0 if disc else 1.3, alpha=0.5 if disc else 1.0,
                   zorder=3)


def _place_labels(ax, pos, niches, disconnected):
    texts = []
    for n in niches:
        if n not in pos:
            continue
        name = get_niche_display(n)[0]
        texts.append(ax.text(pos[n][0], pos[n][1] + 0.03 * KK_SCALE, name,
                             ha="center", va="bottom", fontsize=LABEL_FONTSIZE,
                             fontweight="bold", zorder=5,
                             color="gray" if n in disconnected else "black",
                             bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                                       edgecolor="none", alpha=0.7)))
    if HAS_ADJUSTTEXT and len(texts) > 1:
        try:
            adjust_text(texts, arrowprops=dict(arrowstyle="-", color="lightgray",
                        lw=0.5, alpha=0.6), ax=ax)
        except Exception:
            pass


def _frame_axes(ax, pos):
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    if xs:
        ax.set_xlim(min(xs) - ((max(xs) - min(xs)) * 0.15 + 0.3),
                    max(xs) + ((max(xs) - min(xs)) * 0.15 + 0.3))
        ax.set_ylim(min(ys) - ((max(ys) - min(ys)) * 0.20 + 0.3),
                    max(ys) + ((max(ys) - min(ys)) * 0.20 + 0.3))
    ax.axis("off"); ax.set_aspect("equal")


def plot_network_thresholded(prox_cond_df, niches, title, ax, pos, disconnected,
                             strong_thr, weak_thr, show_legend=False):
    if len(prox_cond_df) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes); ax.axis("off"); return
    for _, r in prox_cond_df.iterrows():
        u, v = r["niche_A"], r["niche_B"]
        if u not in pos or v not in pos or u in disconnected or v in disconnected:
            continue
        p = r["mean_proximity"]
        if p <= strong_thr:
            col, w, al, z = STRONG_EDGE_COLOR, STRONG_EDGE_WIDTH, STRONG_EDGE_ALPHA, 2
        elif p <= weak_thr:
            col, w, al, z = WEAK_EDGE_COLOR, WEAK_EDGE_WIDTH, WEAK_EDGE_ALPHA, 1
        else:
            continue
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color=col,
                linewidth=w, alpha=al, zorder=z, solid_capstyle="round")
    _draw_nodes(ax, pos, niches, disconnected, sizes=niche_sizes)
    _place_labels(ax, pos, niches, disconnected)
    ax.set_title(title, fontsize=12); _frame_axes(ax, pos)
    if show_legend:
        handles = [plt.Line2D([0], [0], color=STRONG_EDGE_COLOR,
                              linewidth=STRONG_EDGE_WIDTH, label="Strongly adjacent"),
                   plt.Line2D([0], [0], color=WEAK_EDGE_COLOR,
                              linewidth=WEAK_EDGE_WIDTH, label="Weakly adjacent")]
        if disconnected:
            handles.append(plt.Line2D([0], [0], marker="o", color="w",
                           markerfacecolor="lightgray", markeredgecolor="gray",
                           markersize=10, label=f"Isolated ({len(disconnected)})"))
        ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=True,
                  framealpha=0.9)


print("\n[PLOT] Per-condition networks (pooled thresholds -> comparable)...")
for cond in CONDITIONS:
    sub = prox_by_cond[prox_by_cond.condition == cond]
    if len(sub) == 0:
        continue
    fig, ax = plt.subplots(figsize=(9, 8))
    plot_network_thresholded(sub, all_niches, COND_SHORT[cond], ax, shared_pos,
                             disconnected, POOLED_STRONG_THR, POOLED_WEAK_THR,
                             show_legend=True)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"niche_network_{COND_SHORT[cond]}.pdf"))
    plt.close()

fig, axes = plt.subplots(1, len(CONDITIONS), figsize=(8.5 * len(CONDITIONS), 7.5))
if len(CONDITIONS) == 1:
    axes = [axes]
for ax, cond in zip(axes, CONDITIONS):
    plot_network_thresholded(prox_by_cond[prox_by_cond.condition == cond],
                             all_niches, COND_SHORT[cond], ax, shared_pos,
                             disconnected, POOLED_STRONG_THR, POOLED_WEAK_THR,
                             show_legend=(cond == CONDITIONS[0]))
fig.suptitle(f"Niche spatial proximity by condition (pooled thresholds, shared "
             f"layout)  ·  node size ∝ niche abundance  ·  {N_NOTE}",
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_network_combined.pdf"))
plt.close()

# =============================================================================
# 9. TOP DIFFERENTIAL EDGES
# =============================================================================
print(f"[PLOT] Top {N_TOP_DIFF_EDGES} differential edges by |Δ proximity|...")
diff_top = diff_df.dropna(subset=["delta_proximity"]).head(N_TOP_DIFF_EDGES).copy()
fig, ax = plt.subplots(figsize=(10, 9))
_draw_nodes(ax, shared_pos, all_niches, disconnected, sizes=niche_sizes)
if len(diff_top):
    mx = diff_top["delta_proximity"].abs().max() + 1
    cmap = plt.cm.RdBu_r
    for _, r in diff_top.iterrows():
        u, v = r["niche_A"], r["niche_B"]
        if u not in shared_pos or v not in shared_pos:
            continue
        color = cmap((-r["delta_proximity"] / mx + 1) / 2)
        w = 1.5 + 4 * abs(r["delta_proximity"]) / mx
        p = r["perm_p"] if not np.isnan(r["perm_p"]) else r["mw_p"]
        al = 0.4 if np.isnan(p) else min(0.95, 0.3 + 0.3 * (-np.log10(max(p, 1e-3))))
        ax.plot([shared_pos[u][0], shared_pos[v][0]],
                [shared_pos[u][1], shared_pos[v][1]], color=color, linewidth=w,
                alpha=al, zorder=2, solid_capstyle="round")
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(-mx, mx))
    plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.04,
                 label="Δ proximity Ctrl−Test (µm)  ·  alpha ∝ −log10(perm p)")
_place_labels(ax, shared_pos, all_niches, disconnected); _frame_axes(ax, shared_pos)
ax.set_title(f"Top {len(diff_top)} differential edges (effect-size ranked)\n"
             f"{N_NOTE} — per-edge inference exploratory", fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_network_differential_top_edges.pdf"))
plt.close()

# =============================================================================
# 10. FOCUS EDGES (hypothesis-driven, FDR over the small set)
# =============================================================================
if FOCUS_EDGES:
    print(f"[FOCUS] {len(FOCUS_EDGES)} pre-specified edges...")
    rows = []
    for edge in FOCUS_EDGES:
        nA, nB = sorted([str(x) for x in edge], key=niche_sort_key)
        row = diff_df[(diff_df.niche_A == nA) & (diff_df.niche_B == nB)]
        if len(row):
            rows.append(row.iloc[0])
        else:
            print(f"  [WARN] edge ({nA},{nB}) absent")
    focus_df = pd.DataFrame(rows)
    if len(focus_df):
        for pc in ["mw_p", "perm_p"]:
            v = focus_df[pc].notna()
            focus_df[pc.replace("_p", "_padj_focus")] = np.nan
            if v.sum():
                focus_df.loc[v, pc.replace("_p", "_padj_focus")] = multipletests(
                    focus_df.loc[v, pc], method="fdr_bh")[1]
        focus_df.to_csv(os.path.join(OUTPUT_DIR, "niche_proximity_focus_edges.csv"),
                        index=False)
        fig, ax = plt.subplots(figsize=(10, 9))
        _draw_nodes(ax, shared_pos, all_niches, disconnected, sizes=niche_sizes)
        mx = focus_df["delta_proximity"].abs().max() + 1
        for _, r in focus_df.iterrows():
            u, v = r["niche_A"], r["niche_B"]
            if u not in shared_pos or v not in shared_pos:
                continue
            color = plt.cm.RdBu_r((-r["delta_proximity"] / mx + 1) / 2)
            sig = (not np.isnan(r.get("perm_padj_focus", np.nan))) and \
                  r["perm_padj_focus"] < FOCUS_FDR_THRESHOLD
            ax.plot([shared_pos[u][0], shared_pos[v][0]],
                    [shared_pos[u][1], shared_pos[v][1]], color=color,
                    linewidth=2 + 4 * abs(r["delta_proximity"]) / mx,
                    alpha=0.95 if sig else 0.45, zorder=2, solid_capstyle="round")
        _place_labels(ax, shared_pos, all_niches, disconnected)
        _frame_axes(ax, shared_pos)
        nsig = int((focus_df.get("perm_padj_focus", pd.Series())
                    .lt(FOCUS_FDR_THRESHOLD)).sum())
        ax.set_title(f"Hypothesis-driven edges (FDR over set; {nsig} sig at "
                     f"FDR<{FOCUS_FDR_THRESHOLD})\nsolid=sig, faded=ns", fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "niche_network_focus_edges.pdf"))
        plt.close()

# =============================================================================
# 11. CONDITION-SCORE-COLOURED NETWORK (nodes shaded by mean condition score score)
# =============================================================================
if HAS_SCORE:
    print("[PLOT] Condition-score-coloured network...")
    sen_by_niche = adata.obs.groupby(NICHE_KEY)[CONDITION_SCORE_KEY].mean()
    vals = np.array([sen_by_niche.get(n, np.nan) for n in all_niches], dtype=float)
    finite = vals[np.isfinite(vals)]
    if len(finite):
        norm = mpl.colors.Normalize(vmin=float(np.nanmin(finite)),
                                    vmax=float(np.nanmax(finite)))
        cmap = plt.cm.magma
        facecol = {n: cmap(norm(sen_by_niche.get(n, np.nan)))
                   if np.isfinite(sen_by_niche.get(n, np.nan)) else "#DDDDDD"
                   for n in all_niches}
        fig, ax = plt.subplots(figsize=(10, 9))
        for _, r in prox_pooled.iterrows():
            u, v = r["niche_A"], r["niche_B"]
            if u in shared_pos and v in shared_pos and u not in disconnected \
                    and v not in disconnected and r["mean_proximity"] <= POOLED_WEAK_THR:
                ax.plot([shared_pos[u][0], shared_pos[v][0]],
                        [shared_pos[u][1], shared_pos[v][1]], color="lightgray",
                        linewidth=1.0, alpha=0.6, zorder=1)
        _draw_nodes(ax, shared_pos, all_niches, disconnected, sizes=niche_sizes,
                    facecolor=facecol)
        _place_labels(ax, shared_pos, all_niches, disconnected)
        _frame_axes(ax, shared_pos)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.04,
                     label="Mean condition score")
        ax.set_title("Senescent niche hubs on the shared spatial layout\n"
                     "(node colour = mean condition score; edges = pooled proximity)",
                     fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "niche_network_condition_score.pdf"))
        plt.close()

# =============================================================================
# 12. SUPPLEMENTS: proximity / contact heatmaps
# =============================================================================
# proximity heatmap (density-sensitive — flagged)
fig, axes = plt.subplots(1, len(CONDITIONS) + 1,
                         figsize=(5.5 * (len(CONDITIONS) + 1), 5))
for ax, cond in zip(axes[:-1], CONDITIONS):
    mat = build_matrix(prox_by_cond[prox_by_cond.condition == cond],
                       all_niches, "mean_proximity")
    sns.heatmap(pd.DataFrame(mat.values, index=display_labels, columns=display_labels),
                cmap="viridis_r", ax=ax, cbar_kws={"label": "Mean proximity (µm)"})
    ax.set_title(COND_SHORT[cond]); plt.setp(ax.get_xticklabels(), rotation=45,
                                             ha="right")
sns.heatmap(delta_display, cmap="RdBu", center=0, vmin=-abs_max, vmax=abs_max,
            ax=axes[-1], cbar_kws={"label": "Ctrl−Test proximity (µm)"})
axes[-1].set_title("Ctrl - Test"); plt.setp(axes[-1].get_xticklabels(), rotation=45,
                                        ha="right")
fig.suptitle("Niche proximity (closest-approach; density-sensitive — read with "
             "contact freq / nhood z)", fontsize=10, y=1.03)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_proximity_heatmap.pdf"))
plt.close()


def build_contact_matrix(sub, niches, col="mean_contact_freq"):
    mat = build_matrix(sub, niches, col)
    np.fill_diagonal(mat.values, 1.0)
    return mat


ctrl_c = build_contact_matrix(prox_by_cond[prox_by_cond.condition == CONDITION_CNTRL],
                            all_niches)
test_c = build_contact_matrix(prox_by_cond[prox_by_cond.condition == CONDITION_TEST],
                            all_niches)
dc = ctrl_c - test_c
dc_disp = pd.DataFrame(dc.values, index=display_labels, columns=display_labels)
abs_c = np.nanmax(np.abs(dc_disp.values))
abs_c = abs_c if np.isfinite(abs_c) and abs_c > 0 else 1.0
fig, axes = plt.subplots(1, len(CONDITIONS) + 1,
                         figsize=(5.5 * (len(CONDITIONS) + 1), 5))
for ax, cond in zip(axes[:-1], CONDITIONS):
    mat = build_contact_matrix(prox_by_cond[prox_by_cond.condition == cond], all_niches)
    sns.heatmap(pd.DataFrame(mat.values, index=display_labels, columns=display_labels),
                cmap="YlOrRd", vmin=0, vmax=1, ax=ax,
                cbar_kws={"label": "Contact frequency"})
    ax.set_title(COND_SHORT[cond]); plt.setp(ax.get_xticklabels(), rotation=45,
                                             ha="right")
sns.heatmap(dc_disp, cmap="RdBu_r", center=0, vmin=-abs_c, vmax=abs_c, ax=axes[-1],
            cbar_kws={"label": "Δ contact freq (Ctrl−Test)"})
axes[-1].set_title("Ctrl - Test"); plt.setp(axes[-1].get_xticklabels(), rotation=45,
                                        ha="right")
fig.suptitle(f"Niche contact frequency (spatial intermingling)  ·  {N_NOTE}",
             fontsize=10, y=1.03)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_contact_frequency_heatmap.pdf"))
plt.close()

diff_df.dropna(subset=["delta_contact_freq"]).sort_values(
    "delta_contact_freq", key=lambda x: x.abs(), ascending=False
).head(N_TOP_DIFF_EDGES).to_csv(
    os.path.join(OUTPUT_DIR, "niche_contact_freq_top_differential.csv"), index=False)

# =============================================================================
# 12b. ARCHITECTURE CHANGE: prevalence shift vs neighbour rewiring + egocentric
#      (replaces the hairball; separates "appears/vanishes" from "rewires")
# =============================================================================
import matplotlib.patches as mpatches

samp_tot_net = adata.obs[SAMPLE_KEY].value_counts()
prop_tab_net = (pd.crosstab(adata.obs[NICHE_KEY], adata.obs[SAMPLE_KEY])
                .reindex(all_niches, fill_value=0).div(samp_tot_net, axis=1))
ctrl_slides_n = [s for s in slides_list if slide_cond[s] == CONDITION_CNTRL]
test_slides_n = [s for s in slides_list if slide_cond[s] == CONDITION_TEST]
ncell_cond = (pd.crosstab(adata.obs[NICHE_KEY], adata.obs[CONDITION_KEY])
              .reindex(all_niches, fill_value=0))
ARCH_MIN_CELLS = 30

arch_rows = []
for n in all_niches:
    mh = float(prop_tab_net.loc[n, ctrl_slides_n].mean())
    ml = float(prop_tab_net.loc[n, test_slides_n].mean())
    nh = int(ncell_cond.loc[n].get(CONDITION_CNTRL, 0))
    nl = int(ncell_cond.loc[n].get(CONDITION_TEST, 0))
    both = nh >= ARCH_MIN_CELLS and nl >= ARCH_MIN_CELLS
    h = ctrl_c.loc[n].drop(n).values.astype(float)
    l = test_c.loc[n].drop(n).values.astype(float)
    ok = ~(np.isnan(h) | np.isnan(l))
    rw = np.nan
    if both and ok.sum() >= 2 and np.nansum(h[ok]) > 0 and np.nansum(l[ok]) > 0:
        rw = 1.0 - float(np.dot(h[ok], l[ok]) /
                         (np.linalg.norm(h[ok]) * np.linalg.norm(l[ok]) + 1e-12))
    arch_rows.append({"niche": n, "mean_prop_Ctrl": mh, "mean_prop_Test": ml,
                      "abund_log2fc": float(np.log2((mh + 1e-4) / (ml + 1e-4))),
                      "n_Ctrl": nh, "n_Test": nl, "present_both": both,
                      "rewiring_cosine": rw})
arch = pd.DataFrame(arch_rows)
arch.to_csv(os.path.join(OUTPUT_DIR, "niche_architecture_change.csv"), index=False)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, max(5, 0.26 * len(all_niches))))
a = arch.sort_values("abund_log2fc")
ya = range(len(a))
ca = [COND_COLORS[CONDITION_CNTRL] if v > 0 else COND_COLORS[CONDITION_TEST]
      for v in a["abund_log2fc"]]
axA.hlines(ya, 0, a["abund_log2fc"], color=ca, lw=2, alpha=0.85)
axA.scatter(a["abund_log2fc"], ya, color=ca, s=26, zorder=3)
axA.set_yticks(list(ya)); axA.set_yticklabels([f"N{n}" for n in a["niche"]], fontsize=7)
for t, n in zip(axA.get_yticklabels(), a["niche"]):
    t.set_color(niche_colors.get(n, "black"))
axA.axvline(0, color="0.6", lw=0.8); axA.set_xlabel("abundance log2(Ctrl/Test)")
axA.set_title("Prevalence shift\n(how often the niche occurs: Ctrl vs Test)", fontsize=10)
for s_ in ("top", "right"):
    axA.spines[s_].set_visible(False)
b = arch[arch["rewiring_cosine"].notna()].sort_values("rewiring_cosine")
yb = range(len(b))
cb = [COND_COLORS[CONDITION_CNTRL] if v > 0 else COND_COLORS[CONDITION_TEST]
      for v in b["abund_log2fc"]]
axB.hlines(yb, 0, b["rewiring_cosine"], color=cb, lw=2, alpha=0.85)
axB.scatter(b["rewiring_cosine"], yb, color=cb, s=26, zorder=3)
axB.set_yticks(list(yb)); axB.set_yticklabels([f"N{n}" for n in b["niche"]], fontsize=7)
for t, n in zip(axB.get_yticklabels(), b["niche"]):
    t.set_color(niche_colors.get(n, "black"))
axB.set_xlabel("neighbour rewiring (cosine distance of contact profile)")
n_excl = int((~arch["rewiring_cosine"].notna()).sum())
axB.set_title(f"Neighbourhood rewiring\n(present in both; {n_excl} too rare in one "
              f"→ excluded)", fontsize=10)
for s_ in ("top", "right"):
    axB.spines[s_].set_visible(False)
fig.suptitle(f"Niche architecture change between conditions  ·  {N_NOTE}\n"
             "left: does the niche become more/less common · right: do its spatial "
             "neighbours change · bar colour = HS- (red) or LS- (blue) enriched",
             fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "niche_architecture_change.pdf"))
plt.close()


def top_partners(cmat, niche, k=6):
    s = cmat.loc[niche].drop(niche)
    s = s[s.notna()].sort_values(ascending=False)
    return [(p, float(s[p])) for p in s.index[:k] if s[p] > 0]


ab_rank = arch.assign(_a=arch["abund_log2fc"].abs()) \
    .sort_values("_a", ascending=False)["niche"].tolist()
rw_rank = arch[arch["rewiring_cosine"].notna()] \
    .sort_values("rewiring_cosine", ascending=False)["niche"].tolist()
ego = list(dict.fromkeys(ab_rank[:4] + rw_rank[:3]))[:6]
arch_i = arch.set_index("niche")
ego = sorted(ego, key=lambda n: -abs(float(arch_i.loc[n, "abund_log2fc"])))

if ego:
    K = len(ego)
    fig = plt.figure(figsize=(13, 1.7 * K + 1.4))
    gs = fig.add_gridspec(K, 3, width_ratios=[1.0, 1.3, 1.3], wspace=0.4,
                          hspace=0.85, top=0.91, bottom=0.05, left=0.04, right=0.985)
    for r, n in enumerate(ego):
        row = arch_i.loc[n]
        axs = fig.add_subplot(gs[r, 0]); axs.axis("off")
        axs.set_xlim(0, 1); axs.set_ylim(0, 1)
        axs.add_patch(mpatches.Rectangle((0.0, 0.82), 0.06, 0.16,
                      color=niche_colors.get(n, "#999"), clip_on=False,
                      transform=axs.transAxes))
        axs.text(0.1, 0.9, f"N{n}", fontsize=13, fontweight="bold", va="center",
                 transform=axs.transAxes)
        is_hs = row["abund_log2fc"] > 0
        axs.text(0.0, 0.66, f"{COND_SHORT[CONDITION_CNTRL] if is_hs else COND_SHORT[CONDITION_TEST]}-enriched · "
                 f"log2FC {row['abund_log2fc']:+.1f}", fontsize=9,
                 color=COND_COLORS[CONDITION_CNTRL] if is_hs
                 else COND_COLORS[CONDITION_TEST], transform=axs.transAxes)
        axs.text(0.0, 0.47, f"cells: {COND_SHORT[CONDITION_CNTRL]} {int(row['n_Ctrl'])} · {COND_SHORT[CONDITION_TEST]} {int(row['n_Test'])}",
                 fontsize=9, color="0.3", transform=axs.transAxes)
        if np.isfinite(row["rewiring_cosine"]):
            axs.text(0.0, 0.28, f"rewiring {row['rewiring_cosine']:.2f}",
                     fontsize=9, color="0.3", transform=axs.transAxes)
        else:
            axs.text(0.0, 0.28, "rewiring n/a · rare in one cond.",
                     fontsize=8.5, color="0.45", transform=axs.transAxes)
        lp, hp = top_partners(test_c, n), top_partners(ctrl_c, n)
        mx = max([v for _, v in lp] + [v for _, v in hp] + [0.01])
        for col, parts0, lab in [(1, lp, "LS neighbours"), (2, hp, "HS neighbours")]:
            ax = fig.add_subplot(gs[r, col])
            parts = parts0[::-1]
            if parts:
                names = [f"N{p}" for p, _ in parts]; vals = [v for _, v in parts]
                ax.barh(range(len(names)), vals, edgecolor="white", linewidth=0.5,
                        color=[niche_colors.get(p, "#bbb") for p, _ in parts])
                ax.set_yticks(range(len(names)))
                ax.set_yticklabels(names, fontsize=8)
                for t, (p, _) in zip(ax.get_yticklabels(), parts):
                    t.set_color(niche_colors.get(p, "black"))
                for i, v in enumerate(vals):
                    ax.text(v, i, f" {v:.2f}", va="center", fontsize=7.5, color="0.3")
            else:
                ax.text(0.5, 0.5, "— none —", ha="center", va="center",
                        transform=ax.transAxes, color="0.6", fontsize=9)
            ax.set_xlim(0, mx * 1.2); ax.tick_params(length=0); ax.set_xticks([])
            for s_ in ("top", "right", "left"):
                ax.spines[s_].set_visible(False)
            if r == 0:
                ax.set_title(lab, fontsize=10, loc="left")
        if r == 0:
            axs.text(0.0, 1.08, "niche · prevalence · cells", fontsize=10,
                     transform=axs.transAxes)
    fig.suptitle(f"Egocentric niche neighbourhoods: Test vs Ctrl  ·  {N_NOTE}\n"
                 "bars = contact frequency to each neighbouring niche; compare the "
                 "two columns to see what each niche sits next to", fontsize=11)
    plt.savefig(os.path.join(OUTPUT_DIR, "niche_egocentric_change.pdf"),
                bbox_inches="tight")
    plt.close()

# =============================================================================
# 13. DONE
# =============================================================================
print("\n[DONE] Step 2.5b complete. Outputs in niche_network/")
print(f"  Architecture change: niche_architecture_change.pdf separates prevalence")
print(f"  shifts (niche becomes more/less common) from neighbour rewiring (niche")
print(f"  keeps its prevalence but changes who it sits next to). Egocentric panels")
print(f"  (niche_egocentric_change.pdf) show LS-vs-HS neighbours for the most-changed.")
print(f"  Mantel r={mantel_r:.3f}, p={mantel_p:.3f} -> high r + low p means the")
print(f"  architecture is CONSERVED across conditions (not reorganised).")
if glob_results:
    g = glob_results[0]
    print(f"  Global reorganisation (slide-level perm, {g['metric']}): "
          f"pseudo-F={g['pseudo_F']:.2f}, p={g['perm_p_F']:.3f}.")
print(f"  Edge thresholds are POOLED so Ctrl/Test panels are comparable. Colours and")
print(f"  layout are shared with 2.5a/2.5c. {N_NOTE}; per-edge p-values exploratory")
print(f"  (use FOCUS_EDGES for FDR-controlled hypothesis tests).")