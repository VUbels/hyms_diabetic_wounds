"""
niche_common.py
===============================================================================
Shared utilities for the niche analysis pipeline (Step 2.5a–2.5d).

This module is the single source of truth for everything that was previously
copy-pasted across analysis_niche_composition.py, analysis_niche_network.py,
analysis_niche_spatial.py and analysis_niche_dossiers.py:

  * obs / obsm / uns key names
  * condition labels, colours, ordering
  * the condition gene signature and its (cached, provenance-checked) scoring
  * shared niche / cell-type colour maps
  * the sample-aware linear mixed model
  * small plotting / ordering / FDR helpers
  * AnnData load + validation + sample/condition bookkeeping
  * gene-program (NicheCompass latent) activity extraction

It also adds two helpers that did not exist before:

  * get_gp_activity()      — one consistent, signed GP-activity matrix
  * relabel_from_source()  — overwrite cell-type labels on an *already trained*
                             NicheCompass object WITHOUT retraining the model.
                             (NicheCompass niches depend only on expression +
                             spatial graph + GP priors + batch covariate, never
                             on cell-type labels, so relabelling is free.)

DESIGN CONTRACT
---------------
Nothing here is hard-wired to a particular biology. The two condition labels,
the gene list, the score quantile and all key names are resolved at IMPORT TIME
from environment variables (with sensible defaults), and may be overridden at
RUNTIME via add_common_args() / apply_common_args().

Because functions resolve their label/key/gene defaults from the module globals
*at call time* (default arg = None → read global inside the body), env vars and
CLI overrides both work regardless of how you import. The only caveat: if you
override the condition LABELS at runtime (--condition-cntrl/--condition-test),
reference the derived objects by attribute (nc.COND_COLORS, nc.CONDITION_ORDER,
nc.CONDITION_CNTRL, ...) rather than `from niche_common import COND_COLORS`, so
you see the rebuilt values. Pure functions (pval_stars, mixedlm_condition, ...)
are safe to `from`-import.

Per-project setup is therefore one of:
    export NICHE_COND_CNTRL=lesion NICHE_COND_TEST=diabetic
    export CONDITION_GENES_FILE=/path/to/genes.txt
or pass --condition-cntrl/--condition-test/--genes on the command line.
===============================================================================
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
    "mixedlm_condition",
    # io / bookkeeping
    "load_adata", "sample_condition_map", "condition_counts",
    "get_gp_activity", "relabel_from_source",
]


# =============================================================================
# PATHS / SHARED CACHE LOCATION
# =============================================================================
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


# =============================================================================
# OBS / OBSM / UNS KEYS  (override via env if your converter names them differently)
# =============================================================================
SPATIAL_KEY   = os.environ.get("NICHE_SPATIAL_KEY",   "X_spatial")
CELLTYPE_KEY  = os.environ.get("NICHE_CELLTYPE_KEY",  "predicted_cell_type")
NICHE_KEY     = os.environ.get("NICHE_NICHE_KEY",     "nichecompass_niche")
SAMPLE_KEY    = os.environ.get("NICHE_SAMPLE_KEY",    "sample_id")
CONDITION_KEY = os.environ.get("NICHE_CONDITION_KEY", "condition")
LATENT_KEY    = os.environ.get("NICHE_LATENT_KEY",    "nichecompass_latent")
ACTIVE_GP_KEY = os.environ.get("NICHE_ACTIVE_GP_KEY", "nichecompass_active_gp_names")


# =============================================================================
# SCORE KEYS / PARAMS
# =============================================================================
CONDITION_SCORE_KEY = os.environ.get("NICHE_SCORE_KEY", "condition_score")
# Canonical flag-column name. (The old scripts disagreed: composition used
# "is_control", the others "is_high_condition". This is the single value now;
# it marks the top-quantile cells of the score = "high condition".)
CONDITION_FLAG_KEY  = os.environ.get("NICHE_FLAG_KEY", "is_high_condition")
CONDITION_SCORE_QUANTILE = float(os.environ.get("CONDITION_SCORE_QUANTILE", "0.75"))


# =============================================================================
# CONDITION LABELS  (+ derived colours / ordering)
# =============================================================================
CONDITION_CNTRL = os.environ.get("NICHE_COND_CNTRL", "control")   # reference / baseline
CONDITION_TEST  = os.environ.get("NICHE_COND_TEST",  "test")      # contrast group

# Colours are configurable; defaults match the existing scripts.
_CNTRL_COLOR = os.environ.get("NICHE_COND_CNTRL_COLOR", "#C94A4A")
_TEST_COLOR  = os.environ.get("NICHE_COND_TEST_COLOR",  "#4A7EB8")


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


# =============================================================================
# GENE SIGNATURE
# =============================================================================
# Default file lives next to this module unless overridden by env / CLI.
CONDITION_GENES_FILE = os.environ.get(
    "CONDITION_GENES_FILE",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "condition_genes.txt"),
)


def load_gene_list(path: str | None) -> list[str]:
    """Load gene symbols from a text file (comma- and/or newline-separated).

    Returns [] (with a warning) if the path is missing or None, so that a score
    is simply skipped downstream rather than crashing the run.
    """
    if not path:
        print("  [SCORE][WARN] No gene file configured — condition score disabled.")
        return []
    try:
        with open(path) as fh:
            raw = fh.read()
    except FileNotFoundError:
        print(f"  [SCORE][WARN] Gene file not found: {path} — score disabled.")
        return []
    genes = [g.strip() for g in raw.replace("\n", ",").split(",") if g.strip()]
    # de-duplicate, preserve order
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    print(f"  [SCORE] Loaded {len(out)} genes from {path}")
    return out


CONDITION_GENES: list[str] = load_gene_list(CONDITION_GENES_FILE)


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


# =============================================================================
# PLOT STYLE
# =============================================================================
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


# =============================================================================
# SMALL HELPERS
# =============================================================================
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


# =============================================================================
# ARGPARSE INTEGRATION
# =============================================================================
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


# =============================================================================
# COLOUR MAPS
# =============================================================================
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


# =============================================================================
# CONDITION SCORE  (cached + provenance-checked)
# =============================================================================
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

    # ---- try cache (must match cells AND the gene signature) ----
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


# =============================================================================
# SAMPLE-AWARE LINEAR MIXED MODEL
# =============================================================================
def mixedlm_condition(df, value_col, sample_col=None, condition_col=None,
                      extra_fixed=None, hi=None, lo=None, max_rows=30000, seed=0):
    """Fit ``value ~ C(condition)[+extra] + (1|sample)``.

    Reference level = ``lo`` (test), so the condition coefficient is the
    Ctrl-vs-Test effect (positive => higher in control). Returns
    dict(effect_Ctrl_vs_Test, p, n, converged). NaN effect when the random
    effect is singular/unidentifiable (honest rather than spurious).
    """
    sample_col = sample_col or SAMPLE_KEY
    condition_col = condition_col or CONDITION_KEY
    hi = hi or CONDITION_CNTRL
    lo = lo or CONDITION_TEST

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
            out.update(effect_Ctrl_vs_Test=float(r.params[cn[0]]),
                       p=float(r.pvalues[cn[0]]))
    except Exception:
        pass
    return out


# =============================================================================
# ANNDATA LOAD / BOOKKEEPING
# =============================================================================
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


# =============================================================================
# GENE-PROGRAM (LATENT) ACTIVITY
# =============================================================================
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


# =============================================================================
# RELABEL WITHOUT RETRAINING
# =============================================================================
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

    # ---- build {integrated_barcode: new_label} ----
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

    # ---- match against integrated obs_names (direct only; no bare-barcode
    #      stripping, which would collide across samples) ----
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


# =============================================================================
# SELF-TEST  (python niche_common.py  -> prints the resolved configuration)
# =============================================================================
if __name__ == "__main__":
    print("niche_common resolved configuration")
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