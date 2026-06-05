"""
==============================================================================
SPATIAL TRANSCRIPTOMICS PIPELINE — STEP 2.5c  (revised)
Per-slide spatial maps: cell types, niches, and condition score
==============================================================================

What changed vs the previous version
------------------------------------
  * Shared cell-type colour map (reused from the R pipeline / shared cache) and
    shared niche colour map, so cell-type and niche spatial maps match the
    colours used in 2.5a / 2.5b exactly.
  * Condition score added to the spatial view:
        - continuous condition score grid (per-slide heatmap on the tissue, with a
          common colour scale across slides for comparability), and
        - high-condition-score cell highlight grid (cells above the shared top-quartile
          threshold lit up), Test arranged above CNTRL
  * A compact condition-score-by-niche summary (Ctrl vs Test, shared niche colours,
    mixed-model significance markers) ties the maps back to the niche analysis.
  * Conditions are ordered Test-above-CNTRL in every grid.

Inputs
------
    nichecompass_results/objects/nichecompass_integrated.h5ad   (or $NICHE_H5AD)

Outputs (spatial_maps/)
-----------------------
    {sample_id}_panel.pdf            per-slide cell-type breakouts
    condition_grid.pdf               cell types, all slides by condition
    niche_grid.pdf                   niches, shared colours
    condition_score_grid.pdf        continuous condition score score on tissue
    high_condition_cells_grid.pdf         high-condition-score cell highlight on tissue
    condition_score_by_niche_summary.pdf  per-niche Ctrl vs Test (mixed-model stars)
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
    "--output_dir", default=None, metavar="DIR",
    help="Output directory (default: spatial_maps)")
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
OUTPUT_DIR  = "spatial_maps"

# CLI overrides (parser runs before this block)
if _args.h5ad:           INPUT_H5AD  = _args.h5ad
if _args.output_dir:     OUTPUT_DIR  = _args.output_dir
if _args.exclude_samples is not None: EXCLUDE_SAMPLES = _args.exclude_samples

SPATIAL_KEY   = "X_spatial"
CELLTYPE_KEY  = "predicted_cell_type"
NICHE_KEY     = "nichecompass_niche"
SAMPLE_KEY    = "sample_id"
CONDITION_KEY = "condition"

EXCLUDE_SAMPLES = [s for s in os.environ.get("EXCLUDE_SAMPLES", "").split(",") if s]

# Cell types to highlight in the breakout panels (None -> top-N most abundant).
HIGHLIGHT_CELLTYPES = None
AUTO_SELECT_TOP_N = 6

POINT_SIZE     = 0.8
POINT_ALPHA    = 0.85
BG_COLOR       = "black"
BG_POINT_COLOR = "#303030"
FIG_DPI        = 150
SCORE_CMAP     = "magma"
SCORE_HIGH_COLOR = "#FFD23F"     # high-condition-score cell highlight

os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# 1. LOAD + CONDITION SCORE + SHARED STYLES
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

samples = sorted(adata.obs[SAMPLE_KEY].unique())
print(f"  {adata.n_obs} cells, {len(samples)} samples")

print("\n[SCORE] Scoring condition score condition score signature...")
compute_condition_score(adata)
score_thr = add_condition_flag(adata)
HAS_SCORE = (CONDITION_SCORE_KEY in adata.obs
           and not adata.obs[CONDITION_SCORE_KEY].isna().all())

niche_colors, _ = build_or_load_niche_style(adata, NICHE_KEY, CELLTYPE_KEY)
ct_colors_all = build_or_load_celltype_colors(adata, CELLTYPE_KEY)
all_niches = sorted(adata.obs[NICHE_KEY].unique(), key=niche_sort_key)

if HIGHLIGHT_CELLTYPES is None:
    HIGHLIGHT_CELLTYPES = adata.obs[CELLTYPE_KEY].value_counts() \
        .head(AUTO_SELECT_TOP_N).index.tolist()
    print(f"  Highlight cell types: {HIGHLIGHT_CELLTYPES}")
ct_colors = {ct: ct_colors_all.get(ct, "#BBBBBB") for ct in HIGHLIGHT_CELLTYPES}

# conditions present, ordered HS-above-LS
conds_present = [c for c in CONDITION_ORDER if c in set(adata.obs[CONDITION_KEY])]
samples_by_cond = {c: sorted(adata.obs.loc[adata.obs[CONDITION_KEY] == c, SAMPLE_KEY]
                             .unique()) for c in conds_present}
n_conds = len(conds_present)
max_per_cond = max((len(v) for v in samples_by_cond.values()), default=1)
n_ctrl = len(samples_by_cond.get(CONDITION_CNTRL, []))
n_test = len(samples_by_cond.get(CONDITION_TEST, []))
N_NOTE = f"n={n_ctrl} HS vs n={n_test} LS slides"

# global score colour scale (robust percentiles) for cross-slide comparability
if HAS_SCORE:
    sv = adata.obs[CONDITION_SCORE_KEY].values
    SCORE_VMIN = float(np.nanpercentile(sv, 2))
    SCORE_VMAX = float(np.nanpercentile(sv, 98))

# =============================================================================
# 2. GENERIC CONDITION GRID
# =============================================================================
def condition_grid(draw_fn, out_path, suptitle, facecolor="white",
                   legend_handles=None, colorbar=None, colorbar_label="",
                   label_color=None):
    """draw_fn(ax, sub) renders one slide. Rows = conditions (HS top), cols =
    slides; shared layout so panels are comparable."""
    lc = label_color or ("white" if facecolor == "black" else "black")
    fig, axes = plt.subplots(n_conds, max_per_cond,
                             figsize=(4 * max_per_cond, 4 * n_conds),
                             facecolor=facecolor, squeeze=False)
    for ri, cond in enumerate(conds_present):
        for ci in range(max_per_cond):
            ax = axes[ri, ci]
            ax.set_facecolor(facecolor)
            cs = samples_by_cond[cond]
            if ci >= len(cs):
                ax.axis("off"); continue
            sub = adata[adata.obs[SAMPLE_KEY] == cs[ci]]
            draw_fn(ax, sub)
            ax.set_title(cs[ci], color=lc, fontsize=9)
            ax.set_aspect("equal"); ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color(lc)
        axes[ri, 0].set_ylabel(COND_SHORT[cond], color=lc, fontsize=12,
                               fontweight="bold")
    if legend_handles:
        fig.legend(handles=legend_handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.02),
                   ncol=min(10, len(legend_handles)),
                   facecolor=facecolor, labelcolor=lc, fontsize=8, frameon=False)
    if colorbar is not None:
        fig.subplots_adjust(right=0.90)
        cax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
        cb = fig.colorbar(colorbar, cax=cax)
        cb.set_label(colorbar_label, color=lc)
        cb.ax.yaxis.set_tick_params(color=lc)
        plt.setp(cb.ax.get_yticklabels(), color=lc)
        fig.suptitle(suptitle, color=lc, fontsize=13, y=1.0)
        plt.savefig(out_path, bbox_inches="tight", dpi=FIG_DPI, facecolor=facecolor)
    else:
        fig.suptitle(suptitle, color=lc, fontsize=13, y=1.0)
        plt.tight_layout()
        plt.savefig(out_path, bbox_inches="tight", dpi=FIG_DPI, facecolor=facecolor)
    plt.close()


# =============================================================================
# 3. PER-SLIDE CELL-TYPE BREAKOUT PANELS
# =============================================================================
def plot_slide_panel(sub, sample_id, condition, out_path):
    coords = np.asarray(sub.obsm[SPATIAL_KEY])
    cts = sub.obs[CELLTYPE_KEY].values
    n_panels = 1 + len(HIGHLIGHT_CELLTYPES)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 4),
                             facecolor=BG_COLOR, squeeze=False)
    axes = axes[0]
    ax = axes[0]; ax.set_facecolor(BG_COLOR)
    bg = ~np.isin(cts, HIGHLIGHT_CELLTYPES)
    if bg.any():
        ax.scatter(coords[bg, 0], coords[bg, 1], c=BG_POINT_COLOR, s=POINT_SIZE,
                   alpha=POINT_ALPHA * 0.4, rasterized=True, linewidths=0)
    for ct in HIGHLIGHT_CELLTYPES:
        m = cts == ct
        if m.any():
            ax.scatter(coords[m, 0], coords[m, 1], c=[ct_colors[ct]], s=POINT_SIZE,
                       alpha=POINT_ALPHA, rasterized=True, linewidths=0)
    ax.set_title(f"{sample_id}\n({COND_SHORT.get(condition, condition)})",
                 color="white", fontsize=10)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("white")
    for ax, ct in zip(axes[1:], HIGHLIGHT_CELLTYPES):
        ax.set_facecolor(BG_COLOR)
        ax.scatter(coords[:, 0], coords[:, 1], c=BG_POINT_COLOR, s=POINT_SIZE,
                   alpha=POINT_ALPHA * 0.3, rasterized=True, linewidths=0)
        m = cts == ct
        if m.any():
            ax.scatter(coords[m, 0], coords[m, 1], c=[ct_colors[ct]],
                       s=POINT_SIZE * 1.4, alpha=POINT_ALPHA, rasterized=True,
                       linewidths=0)
        ax.set_title(ct, color=ct_colors[ct], fontsize=9, fontweight="bold")
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("white")
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight", dpi=FIG_DPI, facecolor=BG_COLOR)
    plt.close()


print("\n[PLOT] Per-slide cell-type panels...")
for sample in samples:
    sub = adata[adata.obs[SAMPLE_KEY] == sample]
    cond = sub.obs[CONDITION_KEY].iloc[0]
    plot_slide_panel(sub, sample, cond,
                     os.path.join(OUTPUT_DIR, f"{sample}_panel.pdf"))

# =============================================================================
# 4. CELL-TYPE CONDITION GRID
# =============================================================================
print("[PLOT] Cell-type condition grid...")


def _draw_celltype(ax, sub):
    coords = np.asarray(sub.obsm[SPATIAL_KEY]); cts = sub.obs[CELLTYPE_KEY].values
    bg = ~np.isin(cts, HIGHLIGHT_CELLTYPES)
    if bg.any():
        ax.scatter(coords[bg, 0], coords[bg, 1], c=BG_POINT_COLOR,
                   s=POINT_SIZE * 0.6, alpha=POINT_ALPHA * 0.35, rasterized=True,
                   linewidths=0)
    for ct in HIGHLIGHT_CELLTYPES:
        m = cts == ct
        if m.any():
            ax.scatter(coords[m, 0], coords[m, 1], c=[ct_colors[ct]], s=POINT_SIZE,
                       alpha=POINT_ALPHA, rasterized=True, linewidths=0)


condition_grid(_draw_celltype,
               os.path.join(OUTPUT_DIR, "condition_grid.pdf"),
               f"Cell-type distributions by condition  ·  {N_NOTE}",
               facecolor=BG_COLOR,
               legend_handles=[plt.Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=c, markersize=8, label=ct)
                               for ct, c in ct_colors.items()])

# =============================================================================
# 5. NICHE GRID (shared colours)
# =============================================================================
print("[PLOT] Niche grid...")


def _draw_niche(ax, sub):
    coords = np.asarray(sub.obsm[SPATIAL_KEY]); niches = sub.obs[NICHE_KEY].values
    for n in all_niches:
        m = niches == n
        if m.any():
            ax.scatter(coords[m, 0], coords[m, 1], c=[niche_colors.get(n, "#999999")],
                       s=POINT_SIZE, alpha=POINT_ALPHA, rasterized=True, linewidths=0)


condition_grid(_draw_niche,
               os.path.join(OUTPUT_DIR, "niche_grid.pdf"),
               f"Niche assignments by condition (shared colours)  ·  {N_NOTE}",
               facecolor="white",
               legend_handles=[plt.Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=niche_colors.get(n, "#999999"),
                               markersize=8, label=f"N{n}") for n in all_niches])

# =============================================================================
# 6. CONDITION SCORE SPATIAL GRIDS
# =============================================================================
if HAS_SCORE:
    print("[PLOT] Condition score grid...")
    norm = mpl.colors.Normalize(vmin=SCORE_VMIN, vmax=SCORE_VMAX)

    def _draw_score(ax, sub):
        coords = np.asarray(sub.obsm[SPATIAL_KEY])
        s = sub.obs[CONDITION_SCORE_KEY].values
        order = np.argsort(s)                       # plot high scores on top
        ax.scatter(coords[order, 0], coords[order, 1], c=s[order], cmap=SCORE_CMAP,
                   norm=norm, s=POINT_SIZE, alpha=POINT_ALPHA, rasterized=True,
                   linewidths=0)

    sm_for_bar = mpl.cm.ScalarMappable(cmap=SCORE_CMAP, norm=norm)
    sm_for_bar.set_array([])

    condition_grid(_draw_score,
                   os.path.join(OUTPUT_DIR, "condition_score_grid.pdf"),
                   f"condition score on tissue (shared scale)  ·  {N_NOTE}",
                   facecolor=BG_COLOR, colorbar=sm_for_bar,
                   colorbar_label="condition score")

    print("[PLOT] Senescent-cell highlight grid...")

    def _draw_senescent(ax, sub):
        coords = np.asarray(sub.obsm[SPATIAL_KEY])
        flag = sub.obs[CONDITION_FLAG_KEY].values.astype(bool)
        if (~flag).any():
            ax.scatter(coords[~flag, 0], coords[~flag, 1], c=BG_POINT_COLOR,
                       s=POINT_SIZE * 0.6, alpha=POINT_ALPHA * 0.3, rasterized=True,
                       linewidths=0)
        if flag.any():
            ax.scatter(coords[flag, 0], coords[flag, 1], c=SCORE_HIGH_COLOR,
                       s=POINT_SIZE * 1.3, alpha=POINT_ALPHA, rasterized=True,
                       linewidths=0)

    q = int(CONDITION_SCORE_QUANTILE * 100)
    condition_grid(_draw_senescent,
                   os.path.join(OUTPUT_DIR, "high_condition_cells_grid.pdf"),
                   f"High-condition-score cells (top {100-q}% condition score, shared threshold)  ·  "
                   f"{N_NOTE}", facecolor=BG_COLOR,
                   legend_handles=[plt.Line2D([0], [0], marker="o", color="w",
                                   markerfacecolor=SCORE_HIGH_COLOR, markersize=8,
                                   label=f"senescent (≥ q{q})")])

    # ---- 6b. condition-score-by-niche summary (mixed-model stars) ----
    print("[PLOT] Condition-score-by-niche summary...")
    cell_df = adata.obs[[CONDITION_SCORE_KEY, NICHE_KEY, SAMPLE_KEY,
                         CONDITION_KEY]].copy()
    cell_df.columns = ["score", "niche", "sample", "condition"]
    rows = []
    for niche in all_niches:
        d = cell_df[cell_df["niche"] == niche]
        res = mixedlm_condition(d, "score", sample_col="sample")
        rows.append({"niche": niche, "mean_score": d["score"].mean(),
                     "mean_Ctrl": d[d.condition == CONDITION_CNTRL]["score"].mean(),
                     "mean_Test": d[d.condition == CONDITION_TEST]["score"].mean(),
                     **res})
    score_niche = pd.DataFrame(rows)
    vp = score_niche["p"].dropna()
    score_niche["padj"] = np.nan
    if len(vp):
        score_niche.loc[vp.index, "padj"] = multipletests(vp, method="fdr_bh")[1]
    score_niche = score_niche.sort_values("mean_score", ascending=False)
    score_niche.to_csv(os.path.join(OUTPUT_DIR, "condition_score_by_niche_summary.csv"),
                     index=False)

    order = score_niche["niche"].tolist()
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(max(7, len(order) * 0.7), 4.2))
    for gi, niche in enumerate(order):
        gd = cell_df[cell_df["niche"] == niche]
        for ci, cond in enumerate(CONDITION_ORDER):
            cd = gd[gd["condition"] == cond]["score"].values
            if len(cd) == 0:
                continue
            xc = gi + (ci - 0.5) * 0.34
            bp = ax.boxplot([cd], positions=[xc], widths=0.28, showfliers=False,
                            patch_artist=True, manage_ticks=False)
            for b in bp["boxes"]:
                b.set(facecolor=COND_COLORS[cond], alpha=0.55, edgecolor="black",
                      linewidth=0.6)
            for med in bp["medians"]:
                med.set(color="black", linewidth=1.0)
        star = pval_stars(score_niche.set_index("niche").loc[niche, "padj"])
        if star:
            ytop = gd["score"].quantile(0.97)
            ax.plot([gi - 0.17, gi + 0.17], [ytop, ytop], color="black", lw=0.8)
            ax.text(gi, ytop, star, ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"N{n}" for n in order])
    for tick, n in zip(ax.get_xticklabels(), order):
        tick.set_color(niche_colors.get(n, "black"))
    ax.set_ylabel("condition score")
    ax.set_title(f"Condition score by niche (Ctrl vs Test)\nstars = mixed-model condition "
                 f"effect, FDR; {N_NOTE}", fontsize=9)
    handles = [plt.Line2D([0], [0], marker="s", color="w",
                          markerfacecolor=COND_COLORS[c], markersize=9,
                          label=COND_SHORT[c]) for c in CONDITION_ORDER]
    ax.legend(handles=handles, frameon=False, loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "condition_score_by_niche_summary.pdf"))
    plt.close()

# =============================================================================
# 7. DONE
# =============================================================================
print("\n[DONE] Step 2.5c complete. Outputs in spatial_maps/")
print(f"  Cell-type and niche colours are shared with 2.5a / 2.5b. Condition score")
print(f"  grids use a shared colour scale / threshold across slides. {N_NOTE}.")