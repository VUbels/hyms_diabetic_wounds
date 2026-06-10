"""
common_py_functions.py
###############################################################################
Single shared module for the whole spatial-analysis pipeline.

Every execution script imports this module and nothing else for shared logic.
It is the one source of truth for what used to be copy-pasted across the
analysis scripts and previously split across more than one shared module:

  * obs / obsm / uns key names
  * condition labels (CNTRL vs TEST), colours, ordering
  * the condition gene set and its cached, provenance-checked per-cell scoring
  * shared niche / cell-type colour maps
  * the sample-aware linear mixed model
  * small plotting / ordering / FDR / effect-size helpers
  * AnnData load + validation + sample/condition bookkeeping
  * gene-program (NicheCompass latent) activity extraction
  * relabel-without-retraining
  * ligand-receptor resources, spatial co-localisation graph + tests
  * pseudobulk aggregation, output/reporting helpers, shared caveats

DESIGN CONTRACT
Nothing here is wired to a particular biology. Two condition groups, a gene
set, a score quantile and all key names are resolved at IMPORT TIME from
environment variables (with sensible defaults) and may be overridden at RUNTIME
via add_common_args() / apply_common_args().

The standardised convention is that obs[CONDITION_KEY] holds the two literal
values "CNTRL" and "TEST", written once up front by set_condition.py. Every
downstream comparison and every log2FC is therefore TEST relative to CNTRL with
no per-script --cntrl/--test flags to keep consistent.

Functions resolve their label/key/gene defaults from the module globals at call
time, so env vars and CLI overrides both work regardless of import style. If you
override the condition labels at runtime, reference derived objects by attribute
(cf.COND_COLORS, cf.CONDITION_ORDER, ...) rather than from-importing them, so you
see the rebuilt values. Pure functions (pval_stars, mixedlm_condition, ...) are
safe to from-import.

Per-project setup is one of:
    export NICHE_COND_CNTRL=CNTRL NICHE_COND_TEST=TEST
    export CONDITION_GENES_FILE=/path/to/genes.txt
or pass --condition-cntrl/--condition-test/--genes on the command line.
###############################################################################
"""

from __future__ import annotations

import os
import glob
import json
import hashlib
import warnings

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt  # noqa: F401  (re-exported convenience)
import seaborn as sns

__all__ = [
    # paths / dirs
    "SHARED_DIR", "NICHE_STYLE_CSV", "CONDITION_SCORE_CSV",
    "CONDITION_SCORE_META", "CELLTYPE_STYLE_CANDIDATES", "DEFAULT_H5AD",
    # obs / obsm / uns keys
    "SPATIAL_KEY", "CELLTYPE_KEY", "NICHE_KEY", "SAMPLE_KEY", "CONDITION_KEY",
    "LATENT_KEY", "ACTIVE_GP_KEY",
    # score keys / params
    "CONDITION_SCORE_KEY", "CONDITION_FLAG_KEY", "CONDITION_SCORE_QUANTILE",
    "CONDITION_GENES_FILE", "CONDITION_GENES",
    # condition labels (derived)
    "CONDITION_CNTRL", "CONDITION_TEST", "CONDITION_ORDER",
    "COND_COLORS", "COND_SHORT",
    # style
    "PUBLICATION_RCPARAMS", "set_plot_style",
    # small helpers
    "niche_sort_key", "pval_stars", "color_yticklabels", "shorten_gp",
    "hierarchical_order", "bh_fdr", "spearman_cols",
    # gene list
    "load_gene_list", "harmonize_genes_to_panel",
    # config / argparse
    "add_common_args", "apply_common_args", "resolve_exclude",
    # styles
    "build_or_load_niche_style", "build_or_load_celltype_colors",
    # scoring
    "compute_condition_score", "add_condition_flag",
    # stats
    "mixedlm_condition", "cliffs_delta", "mannwhitney_test",
    # io / bookkeeping
    "load_adata", "sample_condition_map", "condition_counts",
    "get_gp_activity", "relabel_from_source",
    # generic path resolution + counts handling
    "resolve_input_path", "resolve_counts_layer", "make_lognorm_view",
    # ligand-receptor + spatial co-localisation
    "liana_resource", "split_complex", "ligand_receptor_universe",
    "build_spatial_graph", "colocalization_table", "stable_colocalized_pairs",
    "condition_split_grouping",
    # pseudobulk + output / reporting
    "pseudobulk_matrix", "ensure_out", "banner", "write_report", "CAVEATS",
]


###############################################################################
# PATHS / SHARED CACHE LOCATION
###############################################################################
SHARED_DIR          = os.environ.get("NICHE_SHARED_DIR", "niche_analysis_shared")
NICHE_STYLE_CSV     = os.path.join(SHARED_DIR, "niche_style.csv")
CONDITION_SCORE_CSV = os.path.join(SHARED_DIR, "condition_score.csv")
CONDITION_SCORE_META = os.path.join(SHARED_DIR, "condition_score.meta.json")
CELLTYPE_STYLE_CANDIDATES = [
    "annotated_data/cell_type_colour_map.csv",          # reuse the R-pipeline map
    os.path.join(SHARED_DIR, "celltype_style.csv"),
]
DEFAULT_H5AD = os.environ.get(
    "NICHE_H5AD", "nichecompass_results/objects/nichecompass_integrated.h5ad")

os.makedirs(SHARED_DIR, exist_ok=True)


###############################################################################
# OBS / OBSM / UNS KEYS  (override via env if your converter names them differently)
###############################################################################
SPATIAL_KEY   = os.environ.get("NICHE_SPATIAL_KEY",   "X_spatial")
CELLTYPE_KEY  = os.environ.get("NICHE_CELLTYPE_KEY",  "predicted_cell_type")
NICHE_KEY     = os.environ.get("NICHE_NICHE_KEY",     "nichecompass_niche")
SAMPLE_KEY    = os.environ.get("NICHE_SAMPLE_KEY",    "sample_id")
CONDITION_KEY = os.environ.get("NICHE_CONDITION_KEY", "condition")
LATENT_KEY    = os.environ.get("NICHE_LATENT_KEY",    "nichecompass_latent")
ACTIVE_GP_KEY = os.environ.get("NICHE_ACTIVE_GP_KEY", "nichecompass_active_gp_names")


###############################################################################
# SCORE KEYS / PARAMS
###############################################################################
CONDITION_SCORE_KEY = os.environ.get("NICHE_SCORE_KEY", "condition_score")
# Canonical flag-column name. (The old scripts disagreed: composition used
# "is_control", the others "is_high_condition". This is the single value now;
# it marks the top-quantile cells of the score = "high condition".)
CONDITION_FLAG_KEY  = os.environ.get("NICHE_FLAG_KEY", "is_high_condition")
CONDITION_SCORE_QUANTILE = float(os.environ.get("CONDITION_SCORE_QUANTILE", "0.75"))


###############################################################################
# CONDITION LABELS  (+ derived colours / ordering)
###############################################################################
CONDITION_CNTRL = os.environ.get("NICHE_COND_CNTRL", "CNTRL")   # reference / baseline
CONDITION_TEST  = os.environ.get("NICHE_COND_TEST",  "TEST")    # contrast group

# Colours are configurable; defaults keep TEST warm/red and CNTRL cool/blue,
# matching the figures the niche scripts already produce.
_CNTRL_COLOR = os.environ.get("NICHE_COND_CNTRL_COLOR", "#4A7EB8")
_TEST_COLOR  = os.environ.get("NICHE_COND_TEST_COLOR",  "#C94A4A")


def _rebuild_condition_derived():
    """(Re)build CONDITION_ORDER / COND_COLORS / COND_SHORT from the labels."""
    global CONDITION_ORDER, COND_COLORS, COND_SHORT
    CONDITION_ORDER = [CONDITION_CNTRL, CONDITION_TEST]
    COND_COLORS = {CONDITION_CNTRL: _CNTRL_COLOR, CONDITION_TEST: _TEST_COLOR}
    # Short labels: env override, else "Ctrl"/"Test".
    COND_SHORT = {
        CONDITION_CNTRL: os.environ.get("NICHE_COND_CNTRL_SHORT", "Ctrl"),
        CONDITION_TEST:  os.environ.get("NICHE_COND_TEST_SHORT",  "Test"),
    }


CONDITION_ORDER: list[str] = []
COND_COLORS: dict[str, str] = {}
COND_SHORT: dict[str, str] = {}
_rebuild_condition_derived()


###############################################################################
# GENE SIGNATURE
###############################################################################
# Default file lives next to this module unless overridden by env / CLI.
CONDITION_GENES_FILE = os.environ.get(
    "CONDITION_GENES_FILE",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "condition_genes.txt"),
)


def load_gene_list(path: str | None, verbose: bool = True) -> list[str]:
    """Load gene symbols from a text file (comma- and/or newline-separated).

    Returns [] (with a warning) if the path is missing or None, so that a score
    is simply skipped downstream rather than crashing the run. Pass verbose=False
    to suppress logging (used for the silent import-time default load).
    """
    if not path:
        if verbose:
            print("  [SCORE][WARN] No gene file configured - condition score disabled.")
        return []
    try:
        with open(path) as fh:
            raw = fh.read()
    except FileNotFoundError:
        if verbose:
            print(f"  [SCORE][WARN] Gene file not found: {path} - score disabled.")
        return []
    genes = [g.strip() for g in raw.replace("\n", ",").split(",") if g.strip()]
    # de-duplicate, preserve order
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    if verbose:
        print(f"  [SCORE] Loaded {len(out)} genes from {path}")
    return out


CONDITION_GENES: list[str] = load_gene_list(CONDITION_GENES_FILE, verbose=False)


def harmonize_genes_to_panel(genes, var_names, alias_map=None):
    """Map a gene list onto the symbols actually present in ``var_names``.

    Returns (present, missing, remapped) where ``remapped`` is a {requested:
    resolved} dict for any gene matched only through an alias. ``alias_map`` is
    an optional {old_symbol: new_symbol} dict (e.g. legacy HGNC aliases). Off by
    default — the user controls their own gene file — but available so a project
    can opt into alias resolution without silently dropping genes.
    """
    var = set(map(str, var_names))
    alias_map = alias_map or {}
    present, missing, remapped = [], [], {}
    for g in genes:
        if g in var:
            present.append(g)
        elif g in alias_map and alias_map[g] in var:
            present.append(alias_map[g])
            remapped[g] = alias_map[g]
        else:
            missing.append(g)
    return present, missing, remapped


###############################################################################
# PLOT STYLE
###############################################################################
PUBLICATION_RCPARAMS = {
    "pdf.fonttype": 42, "ps.fonttype": 42,          # editable text in vector output
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
}


def set_plot_style():
    """Apply the publication rcParams. Called once on import; safe to re-call."""
    mpl.rcParams.update(PUBLICATION_RCPARAMS)


set_plot_style()


###############################################################################
# SMALL HELPERS
###############################################################################
def niche_sort_key(x):
    """Sort niches numerically when possible, else lexically (numbers first)."""
    try:
        return (0, int(x))
    except Exception:
        return (1, str(x))


def pval_stars(p, ns=""):
    """Significance stars; returns ``ns`` for NaN/None/non-significant."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ns
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ns


def color_yticklabels(ax, labels, color_map):
    """Colour y tick labels by ``color_map`` keyed on str(label)."""
    for tick, lab in zip(ax.get_yticklabels(), labels):
        tick.set_color(color_map.get(str(lab), "black"))


def shorten_gp(name, maxlen=18):
    """Compact a NicheCompass gene-program name for axis labels."""
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
        base = base[:maxlen - 6] + "\u2026"
    return f"{base} {tag}".strip()


def hierarchical_order(mat_df):
    """Leaf order of an average-linkage tree on the rows (correlation distance)."""
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


# Backwards-compatible private alias (the scripts called it _hierarchical_order).
_hierarchical_order = hierarchical_order


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR that preserves NaNs (NaN p -> NaN q).

    Returns a float ndarray aligned to the input.
    """
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() > 0:
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def spearman_cols(M, y):
    """Spearman rho between each column of M (n x p) and vector y (n).

    Vectorised over columns via rank transform; NaN for degenerate columns.
    """
    from scipy.stats import rankdata
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


###############################################################################
# ARGPARSE INTEGRATION
###############################################################################
def add_common_args(parser):
    """Register the arguments shared by every Step-2.5 script.

    Call this on your script's ArgumentParser, add any script-specific args,
    parse, then call apply_common_args(args). Returns the parser for chaining.
    """
    parser.add_argument(
        "--genes", default=None, metavar="FILE",
        help="Comma/newline-separated gene file to score. Overrides "
             "CONDITION_GENES_FILE and the module default.")
    parser.add_argument(
        "--h5ad", default=None, metavar="FILE",
        help=f"Integrated .h5ad (default: {DEFAULT_H5AD}).")
    parser.add_argument(
        "--output_dir", "--output-dir", dest="output_dir",
        default=None, metavar="DIR", help="Output directory for this script.")
    parser.add_argument(
        "--exclude_samples", "--exclude-samples", dest="exclude_samples",
        nargs="*", default=None, metavar="ID",
        help="Sample IDs to exclude (space-separated).")
    parser.add_argument(
        "--condition_cntrl", "--condition-cntrl", dest="condition_cntrl",
        default=None, help="Control/baseline label in obs[condition].")
    parser.add_argument(
        "--condition_test", "--condition-test", dest="condition_test",
        default=None, help="Test/contrast label in obs[condition].")
    parser.add_argument(
        "--score_quantile", "--score-quantile", dest="score_quantile",
        type=float, default=None, help="Top-quantile threshold for the flag.")
    parser.add_argument(
        "--rebuild_score", "--rebuild-score", dest="rebuild_score",
        action="store_true", help="Ignore the cached condition score.")
    parser.add_argument(
        "--rebuild_style", "--rebuild-style", dest="rebuild_style",
        action="store_true", help="Ignore the cached niche colour map.")
    return parser


def apply_common_args(args):
    """Apply parsed common args to the module-level configuration.

    Mutates the module globals in place (so attribute access sees the new
    values) and rebuilds the derived condition objects. Per-run values that are
    not module config (--h5ad, --output_dir, --exclude_samples, --rebuild_*)
    are left on ``args`` for the caller to read.
    """
    global CONDITION_GENES_FILE, CONDITION_GENES
    global CONDITION_CNTRL, CONDITION_TEST, CONDITION_SCORE_QUANTILE

    if getattr(args, "genes", None):
        CONDITION_GENES_FILE = args.genes
        CONDITION_GENES = load_gene_list(CONDITION_GENES_FILE)
    if getattr(args, "condition_cntrl", None):
        CONDITION_CNTRL = args.condition_cntrl
    if getattr(args, "condition_test", None):
        CONDITION_TEST = args.condition_test
    if getattr(args, "score_quantile", None) is not None:
        CONDITION_SCORE_QUANTILE = float(args.score_quantile)
    _rebuild_condition_derived()


def resolve_exclude(args=None):
    """Resolve the exclude-samples list (CLI > env > []).

    Centralised so the ordering bug in the original scripts (env clobbering the
    CLI value) cannot recur.
    """
    if args is not None and getattr(args, "exclude_samples", None) is not None:
        return list(args.exclude_samples)
    return [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]


###############################################################################
# COLOUR MAPS
###############################################################################
def build_or_load_niche_style(adata, niche_key=None, celltype_key=None,
                              cache_csv=None, rebuild=False):
    """Return (color_dict, ordered_niches).

    Compositionally similar niches get adjacent hues; the result is cached so
    every figure and every script colours the same niche identically.
    """
    niche_key = niche_key or NICHE_KEY
    celltype_key = celltype_key or CELLTYPE_KEY
    cache_csv = cache_csv or NICHE_STYLE_CSV

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
    ordered = hierarchical_order(comp)
    hues = sns.color_palette("husl", n_colors=max(len(ordered), 3))
    color = {n: mpl.colors.to_hex(hues[i]) for i, n in enumerate(ordered)}
    pd.DataFrame({"niche": ordered, "order": range(len(ordered)),
                  "color": [color[n] for n in ordered]}).to_csv(cache_csv, index=False)
    return color, ordered


def build_or_load_celltype_colors(adata, celltype_key=None,
                                  candidates=None, shared_dir=None):
    """Return {cell_type: hex}. Prefers an R-pipeline colour map if present,
    else assigns a stable tab20 palette and writes it to the shared cache."""
    celltype_key = celltype_key or CELLTYPE_KEY
    candidates = candidates if candidates is not None else CELLTYPE_STYLE_CANDIDATES
    shared_dir = shared_dir or SHARED_DIR

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
            os.path.join(shared_dir, "celltype_style.csv"), index=False)
    return {c: color[c] for c in cts}


###############################################################################
# CONDITION SCORE  (cached + provenance-checked)
###############################################################################
def _gene_signature(genes):
    """Stable hash of a gene set (order-independent), for cache provenance."""
    payload = ",".join(sorted(set(map(str, genes))))
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def compute_condition_score(adata, gene_set=None, score_key=None,
                            cache_csv=None, counts_layer="counts", rebuild=False):
    """LogNorm -> sc.tl.score_genes(gene_set), cached per-cell so every script
    shares identical values. Returns the matched gene list.

    The cache is invalidated automatically when the gene set changes (a sidecar
    .meta.json records the matched genes). Numerics are kept identical to the
    original scripts (target_sum 1e4, log1p, ctrl_size 50, random_state 0) so
    pre-existing caches remain valid.
    """
    import scanpy as sc

    gene_set = list(CONDITION_GENES if gene_set is None else gene_set)
    score_key = score_key or CONDITION_SCORE_KEY
    cache_csv = cache_csv or CONDITION_SCORE_CSV
    meta_path = os.path.splitext(cache_csv)[0] + ".meta.json"

    var = set(map(str, adata.var_names))
    present = [g for g in gene_set if g in var]
    missing = [g for g in gene_set if g not in var]
    print(f"  [SCORE] {len(present)}/{len(gene_set)} signature genes present; "
          f"{len(missing)} missing")
    if missing:
        print(f"  [SCORE] missing (curate aliases if needed): "
              f"{', '.join(missing[:25])}" + (" ..." if len(missing) > 25 else ""))

    sig = _gene_signature(present)

    # Try cache (must match cells AND the gene signature)
    if (not rebuild) and os.path.exists(cache_csv):
        try:
            cached_sig = None
            if os.path.exists(meta_path):
                with open(meta_path) as fh:
                    cached_sig = json.load(fh).get("signature")
            s = pd.read_csv(cache_csv, index_col=0)
            cells_ok = score_key in s.columns and adata.obs_names.isin(s.index).all()
            sig_ok = (cached_sig is None) or (cached_sig == sig)
            if cells_ok and sig_ok:
                adata.obs[score_key] = s.loc[adata.obs_names, score_key].values
                tag = "" if cached_sig is not None else " (unverified: no meta)"
                print(f"  [SCORE] loaded cached score from {cache_csv}{tag}")
                return present
            if cells_ok and not sig_ok:
                print("  [SCORE] cached gene set differs from current — recomputing.")
        except Exception as exc:
            print(f"  [SCORE][WARN] cache read failed ({exc}); recomputing.")

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
    with open(meta_path, "w") as fh:
        json.dump({"signature": sig, "n_present": len(present),
                   "genes_present": present, "genes_missing": missing}, fh, indent=2)
    with open(os.path.join(SHARED_DIR, "condition_genes_used.txt"), "w") as fh:
        fh.write("# matched (%d)\n%s\n\n# missing (%d)\n%s\n" %
                 (len(present), "\n".join(present), len(missing), "\n".join(missing)))
    print(f"  [SCORE] score computed and cached -> {cache_csv}")
    return present


def add_condition_flag(adata, score_key=None, flag_key=None, quantile=None):
    """Flag the top-quantile cells (high condition score). Returns the threshold."""
    score_key = score_key or CONDITION_SCORE_KEY
    flag_key = flag_key or CONDITION_FLAG_KEY
    quantile = CONDITION_SCORE_QUANTILE if quantile is None else quantile
    if score_key not in adata.obs or adata.obs[score_key].isna().all():
        adata.obs[flag_key] = False
        return np.nan
    thr = float(np.nanquantile(adata.obs[score_key].values, quantile))
    adata.obs[flag_key] = (adata.obs[score_key].values >= thr)
    return thr


###############################################################################
# SAMPLE-AWARE LINEAR MIXED MODEL
###############################################################################
def mixedlm_condition(df, value_col, sample_col=None, condition_col=None,
                      extra_fixed=None, hi=None, lo=None, max_rows=30000, seed=0):
    """Fit ``value ~ C(condition)[+extra] + (1|sample)``.

    Reference level = ``lo`` (test), so the condition coefficient is the
    Ctrl-vs-Test effect (positive => higher in control). Returns
    dict(effect_CNTRL_vs_TEST, p, n, converged). NaN effect when the random
    effect is singular/unidentifiable (honest rather than spurious).
    """
    sample_col = sample_col or SAMPLE_KEY
    condition_col = condition_col or CONDITION_KEY
    hi = hi or CONDITION_CNTRL
    lo = lo or CONDITION_TEST

    out = {"effect_CNTRL_vs_TEST": np.nan, "p": np.nan, "n": int(len(df)),
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
    # deterministic order so the same cells give an identical fit
    d = d.sort_values(["_grp", condition_col, value_col]).reset_index(drop=True)

    rhs = f"C({condition_col})"
    for c in (extra_fixed or []):
        rhs += f" + C({c})" if str(d[c].dtype) in ("object", "category") else f" + {c}"

    try:
        md = smf.mixedlm(f"{value_col} ~ {rhs}", d, groups=d["_grp"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = md.fit(reml=False, method="lbfgs", maxiter=200, disp=False)
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
            out.update(effect_CNTRL_vs_TEST=float(r.params[cn[0]]),
                       p=float(r.pvalues[cn[0]]))
    except Exception:
        pass
    return out


###############################################################################
# ANNDATA LOAD / BOOKKEEPING
###############################################################################
_REQUIRE_TO_KEY = {
    "spatial":   ("obsm", lambda: SPATIAL_KEY),
    "celltype":  ("obs",  lambda: CELLTYPE_KEY),
    "niche":     ("obs",  lambda: NICHE_KEY),
    "sample":    ("obs",  lambda: SAMPLE_KEY),
    "condition": ("obs",  lambda: CONDITION_KEY),
    "latent":    ("obsm", lambda: LATENT_KEY),
}


def load_adata(h5ad=None, require=("spatial", "celltype", "niche", "sample",
                                   "condition"),
               exclude_samples=None, coerce_str=True, verbose=True):
    """Read the integrated h5ad, validate required keys, exclude samples, and
    coerce the categorical obs columns to str. Returns the AnnData.
    """
    import scanpy as sc
    path = h5ad or DEFAULT_H5AD
    if verbose:
        print(f"[LOAD] {path}")
    adata = sc.read_h5ad(path)

    for name in require:
        where, key_fn = _REQUIRE_TO_KEY[name]
        key = key_fn()
        store = adata.obsm if where == "obsm" else adata.obs
        present = (key in store) if where == "obsm" else (key in store.columns)
        if not present:
            raise KeyError(f"Missing {where}['{key}'] (required='{name}'). "
                           f"Available {where}: {list(store)[:12]}...")

    if exclude_samples:
        keep = ~adata.obs[SAMPLE_KEY].astype(str).isin([str(s) for s in exclude_samples])
        n_drop = int((~keep).sum())
        if verbose:
            print(f"  Excluding {list(exclude_samples)}: dropping {n_drop} cells")
        adata = adata[keep].copy()

    if coerce_str:
        for k in (NICHE_KEY, CELLTYPE_KEY, CONDITION_KEY, SAMPLE_KEY):
            if k in adata.obs.columns:
                adata.obs[k] = adata.obs[k].astype(str)

    if verbose:
        n_niche = adata.obs[NICHE_KEY].nunique() if NICHE_KEY in adata.obs else 0
        n_ct = adata.obs[CELLTYPE_KEY].nunique() if CELLTYPE_KEY in adata.obs else 0
        print(f"  {adata.n_obs} cells | {n_niche} niches | {n_ct} cell types")
    return adata


def sample_condition_map(adata, sample_key=None, condition_key=None):
    """{sample_id: condition_label} using the first observed condition per sample."""
    sample_key = sample_key or SAMPLE_KEY
    condition_key = condition_key or CONDITION_KEY
    return (adata.obs.groupby(sample_key)[condition_key].first().to_dict())


def condition_counts(adata, sample_key=None, condition_key=None):
    """Return (n_cntrl, n_test, n_note) at the SAMPLE level.

    ``n_note`` is the standard caption string, e.g. "n=3 Ctrl vs n=3 Test slides".
    """
    cmap = sample_condition_map(adata, sample_key, condition_key)
    n_cntrl = sum(1 for v in cmap.values() if v == CONDITION_CNTRL)
    n_test = sum(1 for v in cmap.values() if v == CONDITION_TEST)
    note = (f"n={n_cntrl} {COND_SHORT.get(CONDITION_CNTRL, 'Ctrl')} vs "
            f"n={n_test} {COND_SHORT.get(CONDITION_TEST, 'Test')} slides")
    return n_cntrl, n_test, note


###############################################################################
# GENE-PROGRAM (LATENT) ACTIVITY
###############################################################################
def get_gp_activity(adata, latent_key=None, active_gp_key=None, signed=True):
    """Return a per-cell gene-program activity DataFrame (cells x GPs).

    The NicheCompass latent space is one feature per active GP. When the number
    of active GP names matches the latent dimension we label columns by GP name;
    otherwise we fall back to generic latent_i labels and warn. ``signed=True``
    keeps the natural sign (correct for Spearman vs a score); set False for
    magnitude-only use.

    Returns None if no latent embedding is present.
    """
    latent_key = latent_key or LATENT_KEY
    active_gp_key = active_gp_key or ACTIVE_GP_KEY
    if latent_key not in adata.obsm:
        return None
    L = np.asarray(adata.obsm[latent_key])
    if not signed:
        L = np.abs(L)
    names = list(adata.uns.get(active_gp_key, []))
    if len(names) == L.shape[1]:
        cols = names
    else:
        cols = [f"latent_{i}" for i in range(L.shape[1])]
        if names:
            print(f"  [GP][WARN] active GP names ({len(names)}) != latent dim "
                  f"({L.shape[1]}); using generic latent labels")
    return pd.DataFrame(L, index=adata.obs_names, columns=cols)


###############################################################################
# RELABEL WITHOUT RETRAINING
###############################################################################
def _read_obs_only(path):
    """Read just .obs from an h5ad as cheaply as possible (backed, fallback)."""
    import anndata as ad
    try:
        a = ad.read_h5ad(path, backed="r")
        obs = a.obs.copy()
        try:
            a.file.close()
        except Exception:
            pass
        return obs
    except Exception:
        return ad.read_h5ad(path).obs.copy()


def relabel_from_source(adata, source, *, celltype_key=None, sample_key=None,
                        label_col=None, barcode_col=None, mode="auto",
                        index_prefix_sep="_", min_overlap=0.99,
                        out_path=None, keep_previous=True, verbose=True):
    """Overwrite cell-type labels on an *already trained* integrated object.

    NicheCompass niches depend only on expression + spatial graph + GP priors +
    batch covariate, never on cell-type labels. So to adopt a new reference /
    relabelling you do NOT retrain: load the trained object, swap the cell-type
    column, and re-run the analysis scripts. This helper does the swap safely.

    Parameters
    ----------
    adata : AnnData | str
        The trained integrated object (or path to its .h5ad).
    source : str
        Either (a) a CSV with a barcode column/index and a label column, or
        (b) a directory of per-sample ``*_annotated.h5ad`` files (the layout
        produced upstream). In the directory case, the integrated barcodes are
        assumed to be ``f"{sample_id}{sep}{original_barcode}"`` (matching the
        NicheCompass concat convention), where ``sample_id`` is the parent
        folder name.
    mode : {"auto","csv","annotated_dir"}
    label_col : str
        Column in the source holding the new label (default: celltype_key).
    barcode_col : str | None
        Barcode column in a CSV source; if None, the CSV index is used.
    min_overlap : float
        Require at least this fraction of integrated cells to be matched, else
        raise (prevents silently corrupting labels through a key mismatch).

    Returns
    -------
    AnnData with obs[celltype_key] replaced (and obs[celltype_key + "_prev"]
    preserving the old labels if keep_previous). Writes to ``out_path`` if given.
    """
    import anndata as ad

    celltype_key = celltype_key or CELLTYPE_KEY
    sample_key = sample_key or SAMPLE_KEY
    label_col = label_col or celltype_key

    if isinstance(adata, str):
        adata = ad.read_h5ad(adata)

    if mode == "auto":
        mode = "csv" if str(source).lower().endswith(".csv") else "annotated_dir"

    # Build {integrated_barcode: new_label}
    mapping = {}
    if mode == "csv":
        df = pd.read_csv(source)
        if label_col not in df.columns:
            raise KeyError(f"label_col '{label_col}' not in CSV columns {list(df.columns)}")
        # Raw 10x barcodes repeat across samples, so a bare-barcode match would
        # collide. Prefer (sample, barcode) reconstruction when both columns are
        # present (builds 'sample{sep}barcode' keys to match the integrated
        # obs_names). Otherwise treat a single id column as the FULL integrated
        # obs_name (which already carries the sample prefix) and match directly.
        lc = {c.lower(): c for c in df.columns}
        scol = lc.get(sample_key.lower(), lc.get("sample", lc.get("sample_id")))
        bcol = barcode_col or lc.get("barcode", lc.get("cell_id", lc.get("cellid")))
        if scol in df.columns and (bcol in df.columns if bcol else False):
            keys = df[scol].astype(str) + index_prefix_sep + df[bcol].astype(str)
            mapping = dict(zip(keys, df[label_col].astype(str)))
        else:
            id_col = barcode_col or df.columns[0]
            mapping = dict(zip(df[id_col].astype(str), df[label_col].astype(str)))

    elif mode == "annotated_dir":
        files = sorted(glob.glob(os.path.join(source, "*", "*_annotated.h5ad")))
        if not files:
            files = sorted(glob.glob(os.path.join(source, "*_annotated.h5ad")))
        if not files:
            raise FileNotFoundError(f"No *_annotated.h5ad under {source}")
        for f in files:
            sid = os.path.basename(os.path.dirname(f))
            obs = _read_obs_only(f)
            if label_col not in obs.columns:
                if verbose:
                    print(f"  [RELABEL][WARN] {os.path.basename(f)} lacks "
                          f"'{label_col}'; skipping")
                continue
            for bc, lab in zip(obs.index.astype(str), obs[label_col].astype(str)):
                mapping[f"{sid}{index_prefix_sep}{bc}"] = lab
    else:
        raise ValueError(f"unknown mode '{mode}'")

    # Match against integrated obs_names (direct only; no bare-barcode
    #      stripping, which would collide across samples)
    names = adata.obs_names.astype(str)
    new = names.map(mapping)
    overlap = float(new.notna().mean())

    if overlap < min_overlap:
        unmatched = [n for n, v in zip(names, new) if pd.isna(v)][:8]
        raise ValueError(
            f"relabel matched only {overlap:.1%} of cells (need {min_overlap:.0%}). "
            f"Keys must equal the integrated obs_names (with the sample prefix), "
            f"or supply sample + barcode columns so they can be reconstructed. "
            f"Example integrated names: {list(names[:4])}; "
            f"example source keys: {list(mapping)[:4]}; "
            f"example unmatched: {unmatched}")

    if verbose:
        print(f"  [RELABEL] matched {overlap:.2%} of {adata.n_obs} cells; "
              f"{new.nunique()} label categories")

    if keep_previous and celltype_key in adata.obs.columns:
        adata.obs[celltype_key + "_prev"] = adata.obs[celltype_key].astype(str).values
    adata.obs[celltype_key] = pd.Categorical(new.astype(str).values)

    if out_path:
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        adata.write_h5ad(out_path)
        if verbose:
            print(f"  [RELABEL] wrote {out_path}")
    return adata


###############################################################################
# GENERIC PATH RESOLUTION + COUNTS HANDLING
###############################################################################
def resolve_input_path(path, *, required=True, what="file"):
    """Resolve a user-supplied path to an existing absolute path.

    Honoured exactly as given: '~' is expanded and relative paths are taken
    relative to the current working directory. Nothing is searched next to the
    script. Returns the absolute path, or None when path is empty. If a path is
    given but missing: raises FileNotFoundError (required=True) or warns and
    returns None (required=False), naming the absolute path tried and the CWD.
    """
    if not path:
        return None
    resolved = os.path.abspath(os.path.expanduser(str(path)))
    if os.path.exists(resolved):
        return resolved
    msg = (f"Could not find {what}: '{path}'\n"
           f"  looked for: {resolved}\n"
           f"  working dir: {os.getcwd()}\n"
           f"  pass a path relative to your working directory or an absolute path.")
    if required:
        raise FileNotFoundError(msg)
    warnings.warn(msg)
    return None


def resolve_counts_layer(adata, counts_layer="counts", verbose=True):
    """Return the name of a layer holding raw integer-like counts.

    Prefers the named layer; falls back to .X if it looks like raw counts.
    Returns None if neither is integer-like (a normalised object), so callers can
    decide whether to proceed.
    """
    import scipy.sparse as sp

    def _looks_integer(M):
        if M is None:
            return False
        sample = M[:200] if M.shape[0] > 200 else M
        sample = sample.toarray() if sp.issparse(sample) else np.asarray(sample)
        if sample.size == 0:
            return False
        finite = sample[np.isfinite(sample)]
        if finite.size == 0:
            return False
        return bool(np.allclose(finite, np.round(finite))) and float(finite.min()) >= 0.0

    if counts_layer and counts_layer in getattr(adata, "layers", {}):
        if verbose:
            print(f"  [COUNTS] using layer '{counts_layer}'")
        return counts_layer
    if _looks_integer(adata.X):
        if verbose:
            print("  [COUNTS] using .X (looks like raw counts)")
        return None  # signal: use .X
    raise ValueError(
        "No raw-count layer found and .X is not integer-like. Pass --counts-layer "
        "with the name of the layer holding raw counts.")


def make_lognorm_view(adata, counts_layer="counts", target_sum=1e4):
    """Return a copy of adata with .X log-normalised from the raw counts, leaving
    the raw counts untouched in the original object. Used for L-R inference."""
    import scanpy as sc
    tmp = adata.copy()
    if counts_layer and counts_layer in tmp.layers:
        tmp.X = tmp.layers[counts_layer].copy()
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    sc.pp.log1p(tmp)
    return tmp


###############################################################################
# LIGAND-RECEPTOR RESOURCE
###############################################################################
def liana_resource(resource_name="consensus"):
    """Return the LIANA ligand-receptor resource as a DataFrame with standardised
    'ligand' / 'receptor' columns."""
    import liana as li
    res = li.rs.select_resource(resource_name)
    cols = {c.lower(): c for c in res.columns}
    return res.rename(columns={cols.get("ligand", "ligand"): "ligand",
                               cols.get("receptor", "receptor"): "receptor"})


def split_complex(token):
    """A LIANA complex like 'IL6_IL6R' or 'ITGA1&ITGB1' -> list of subunit symbols."""
    if not isinstance(token, str):
        return []
    for sep in ("_", "&", "+"):
        if sep in token:
            return [t for t in token.split(sep) if t]
    return [token]


def ligand_receptor_universe(resource_name="consensus"):
    """Return (ligand_set, receptor_set) of individual gene symbols in the resource."""
    res = liana_resource(resource_name)
    ligands, receptors = set(), set()
    for v in res["ligand"].astype(str):
        ligands.update(split_complex(v))
    for v in res["receptor"].astype(str):
        receptors.update(split_complex(v))
    return ligands, receptors


###############################################################################
# SPATIAL GRAPH + CO-LOCALISATION (the spatial constraint)
###############################################################################
def build_spatial_graph(adata, spatial_key=None, n_neighs=6, coord_type="generic"):
    """Build a kNN spatial neighbour graph with squidpy (in place)."""
    import squidpy as sq
    spatial_key = spatial_key or SPATIAL_KEY
    if "spatial" not in adata.obsm and spatial_key in adata.obsm:
        adata.obsm["spatial"] = np.asarray(adata.obsm[spatial_key])
    sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)
    return adata


def colocalization_table(adata, sample_key=None, celltype_key=None,
                         spatial_key=None, n_neighs=6, metric="contact_fraction"):
    """Per (sample, source_type, target_type) spatial co-occurrence.

    For metric='contact_fraction' the value is the fraction of target-type cells
    that have at least one source-type cell in their spatial neighbourhood
    (computed per sample, then stacked). 'nhood_zscore' uses squidpy's
    neighbourhood-enrichment z-score instead. Returns a long DataFrame with
    columns [sample, source, target, value, metric].
    """
    import scipy.sparse as sp
    import squidpy as sq
    sample_key = sample_key or SAMPLE_KEY
    celltype_key = celltype_key or CELLTYPE_KEY
    spatial_key = spatial_key or SPATIAL_KEY

    rows = []
    for samp in adata.obs[sample_key].astype(str).unique():
        sub = adata[adata.obs[sample_key].astype(str) == samp].copy()
        if sub.n_obs < 10:
            continue
        build_spatial_graph(sub, spatial_key=spatial_key, n_neighs=n_neighs)
        labels = sub.obs[celltype_key].astype(str).values
        types = sorted(set(labels))
        if metric == "nhood_zscore":
            sub.obs["_ct"] = pd.Categorical(labels, categories=types)
            sq.gr.nhood_enrichment(sub, cluster_key="_ct", show_progress_bar=False)
            z = sub.uns["_ct_nhood_enrichment"]["zscore"]
            for i, s in enumerate(types):
                for j, t in enumerate(types):
                    rows.append((samp, s, t, float(z[i, j]), metric))
        else:
            A = sub.obsp["spatial_connectivities"]
            A = (A > 0).astype(int)
            onehot = pd.get_dummies(pd.Series(labels, index=sub.obs_names))
            has_src = pd.DataFrame((A @ onehot.values) > 0,
                                   index=sub.obs_names, columns=onehot.columns).astype(int)
            for t in types:
                tgt_mask = labels == t
                if tgt_mask.sum() == 0:
                    continue
                for s in types:
                    if s not in has_src.columns:
                        frac = 0.0
                    else:
                        frac = float(has_src.loc[tgt_mask, s].mean())
                    rows.append((samp, s, t, frac, metric))
    return pd.DataFrame(rows, columns=["sample", "source", "target", "value", "metric"])


def stable_colocalized_pairs(coloc_df, threshold=0.05, min_samples_frac=1 / 3.0):
    """From a colocalization_table, return the set of (source, target) pairs that
    are co-localised (value > threshold) in at least min_samples_frac of samples."""
    if coloc_df is None or coloc_df.empty:
        return set()
    n_samples = coloc_df["sample"].nunique()
    need = max(1, int(np.ceil(min_samples_frac * n_samples)))
    hit = coloc_df[coloc_df["value"] > threshold]
    counts = hit.groupby(["source", "target"])["sample"].nunique()
    return set(counts[counts >= need].index)


def condition_split_grouping(adata, group_key, flag_key=None, out_key=None,
                             hi="cond_high", lo="cond_low"):
    """Build a grouping column that splits each `group_key` level by the condition
    flag, e.g. 'Fibroblast|cond_high' vs 'Fibroblast|cond_low'.

    Generalises 'treat condition-high cells of a type as a distinct sender'. Uses
    the boolean CONDITION_FLAG_KEY (see add_condition_flag). Returns the new
    column name written to obs.
    """
    flag_key = flag_key or CONDITION_FLAG_KEY
    out_key = out_key or f"{group_key}__condsplit"
    base = adata.obs[group_key].astype(str)
    if flag_key in adata.obs:
        tag = np.where(adata.obs[flag_key].astype(bool).values, hi, lo)
    else:
        tag = np.array([lo] * adata.n_obs)
    adata.obs[out_key] = (base.values.astype(object) + "|" + tag).astype(str)
    return out_key


###############################################################################
# EFFECT SIZE + NON-PARAMETRIC TEST  (slide-as-unit; report effect over p)
###############################################################################
def cliffs_delta(a, b):
    """Cliff's delta in [-1, 1]; +1 means group a stochastically dominates b.

    For n=3 vs n=3 designs prefer this (and log2FC) over a p-value, which floors
    near 0.1 for a two-sided rank test at that sample size.
    """
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return np.nan
    gt = int(sum((x > b).sum() for x in a))
    lt = int(sum((x < b).sum() for x in a))
    return (gt - lt) / (a.size * b.size)


def mannwhitney_test(values_cntrl, values_test):
    """Two-sided Mann-Whitney U with Cliff's delta (TEST vs CNTRL).

    Returns dict(pvalue, cliffs_delta, median_cntrl, median_test, n_cntrl, n_test).
    Positive delta => TEST tends to exceed CNTRL.
    """
    from scipy.stats import mannwhitneyu
    a = np.asarray(values_test, dtype=float); a = a[np.isfinite(a)]
    b = np.asarray(values_cntrl, dtype=float); b = b[np.isfinite(b)]
    out = {"pvalue": np.nan, "cliffs_delta": np.nan,
           "median_cntrl": float(np.median(b)) if b.size else np.nan,
           "median_test": float(np.median(a)) if a.size else np.nan,
           "n_cntrl": int(b.size), "n_test": int(a.size)}
    if a.size and b.size:
        out["cliffs_delta"] = cliffs_delta(a, b)
        try:
            out["pvalue"] = float(mannwhitneyu(a, b, alternative="two-sided")[1])
        except ValueError:
            pass
    return out


###############################################################################
# PSEUDOBULK AGGREGATION
###############################################################################
def pseudobulk_matrix(adata, group_keys, counts_layer="counts", min_cells=10):
    """Sum raw counts within each combination of `group_keys` (e.g.
    [sample, cell_type]).

    Returns (counts_df, coldata) where counts_df is genes x pseudobulk-sample and
    coldata is indexed by the same pseudobulk ids ('a||b' joined group values)
    with the group columns plus a '_n_cells' column. Groups with < min_cells are
    dropped.
    """
    import scipy.sparse as sp
    layer = counts_layer if (counts_layer and counts_layer in adata.layers) else None
    X = adata.layers[layer] if layer else adata.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    obs = adata.obs[list(group_keys)].astype(str)
    ids, mats, meta = [], [], []
    for key_vals, idx in obs.groupby(list(group_keys)).groups.items():
        if not isinstance(key_vals, tuple):
            key_vals = (key_vals,)
        pos = [adata.obs_names.get_loc(i) for i in idx]
        if len(pos) < min_cells:
            continue
        ids.append("||".join(map(str, key_vals)))
        mats.append(np.asarray(X[pos].sum(axis=0)).ravel())
        meta.append(list(key_vals) + [len(pos)])
    if not ids:
        return (pd.DataFrame(index=adata.var_names),
                pd.DataFrame(columns=list(group_keys) + ["_n_cells"]))
    counts_df = pd.DataFrame(np.vstack(mats), index=ids, columns=adata.var_names).T
    coldata = pd.DataFrame(meta, index=ids, columns=list(group_keys) + ["_n_cells"])
    return counts_df, coldata


###############################################################################
# OUTPUT / REPORTING
###############################################################################
def ensure_out(output_dir):
    """Create output_dir (and a figures/ subdir) and return its path."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "figures"), exist_ok=True)
    return output_dir


def banner(title):
    """Print a solid-# section banner to stdout."""
    print("\n" + "#" * 79 + f"\n# {title}\n" + "#" * 79)


def write_report(output_dir, payload, filename):
    """Write a JSON run-report into output_dir."""
    path = os.path.join(output_dir, filename)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    print(f"  [REPORT] {path}")
    return path


###############################################################################
# SHARED INTERPRETATION CAVEATS  (carry into figure legends)
###############################################################################
CAVEATS = (
    "Interpretation caveats:\n"
    "  1. Circularity: if the CNTRL/TEST split was defined from the same gene set "
    "being scored/tested here, 'TEST has more signal' is partly circular - report "
    "WHERE (which cells/niches) the signal concentrates, not merely that TEST "
    "exceeds CNTRL, unless the split came from an independent source.\n"
    "  2. n is small: the slide is the unit. A two-sided Mann-Whitney floors near "
    "p approx 0.1 for 3 vs 3, so lead with effect size (Cliff's delta, log2FC, "
    "factor-loading separation); treat p as supporting, not gating.\n"
    "  3. Donor confound: with few donors per arm, condition is confounded with "
    "any donor-level covariate; per-sample pseudobulk does not remove a group-level "
    "confound.\n"
    "  4. Contact is not communication: the spatial constraint removes pairs that "
    "never touch; co-localisation is necessary, not sufficient, for signalling.\n"
    "  5. Direction convention: every log2FC and 'enriched in' statement is TEST "
    "relative to CNTRL (positive = higher in TEST)."
)


###############################################################################
# SELF-TEST  (python common_py_functions.py -> prints the resolved configuration)
###############################################################################
if __name__ == "__main__":
    print("common_py_functions resolved configuration")
    print("-" * 60)
    print(f"  SHARED_DIR        = {SHARED_DIR}")
    print(f"  DEFAULT_H5AD      = {DEFAULT_H5AD}")
    print(f"  keys: spatial={SPATIAL_KEY} celltype={CELLTYPE_KEY} "
          f"niche={NICHE_KEY} sample={SAMPLE_KEY} condition={CONDITION_KEY}")
    print(f"        latent={LATENT_KEY} active_gp={ACTIVE_GP_KEY}")
    print(f"  score_key={CONDITION_SCORE_KEY} flag_key={CONDITION_FLAG_KEY} "
          f"quantile={CONDITION_SCORE_QUANTILE}")
    print(f"  conditions: CNTRL={CONDITION_CNTRL!r} TEST={CONDITION_TEST!r}")
    print(f"  COND_ORDER={CONDITION_ORDER} COND_SHORT={COND_SHORT}")
    print(f"  COND_COLORS={COND_COLORS}")
    print(f"  gene file = {CONDITION_GENES_FILE}")
    print(f"  genes loaded = {len(CONDITION_GENES)}")
