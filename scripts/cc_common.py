"""
cc_common.py
============
Shared utilities for the two communication-layer analyses:

  1. analysis_cellcell_communication.py   (spatially-constrained LIANA+ -> Tensor-cell2cell)
  2. analysis_condition_factor_roles.py    (per-factor pseudobulk DE + dual-role assignment)

Design principles (deliberate, so these scripts are reusable across datasets):

  * Everything is phrased as CNTRL vs TEST, never as a hard-coded biological label.
    For the senescence/wound instantiation, set --cntrl LS --test HS, but the code
    never assumes that. CNTRL is the reference/baseline level; TEST is the perturbed
    level. log2FC and "enriched in" are always reported as TEST-relative-to-CNTRL.

  * The signature is supplied as a flat gene list (condition_genes.txt), never
    hard-coded as "senescence". The factor-role script decomposes THAT list into
    individual signalling factors. Swap the file to repurpose the pipeline.

  * Cell-type / sample / condition / niche obs keys, the spatial obsm key and the
    raw-counts layer are all resolved in one place, with this precedence:
        explicit CLI argument  >  attribute on the project's niche_common module
        >  auto-detection against a candidate list  >  documented default.
    This means if the project's niche_common.py defines e.g. CELLTYPE_KEY, that
    value wins automatically and the two scripts stay consistent with the rest of
    the pipeline without importing anything fragile.

  * Heavy third-party imports (scanpy, squidpy, liana, cell2cell, pydeseq2,
    decoupler) are done lazily inside the functions that need them, so `--help`
    and unit-level use never require the full stack to be installed.

Author: analysis pipeline (communication layer)
"""

from __future__ import annotations

import json
import os
import sys
import textwrap
from dataclasses import dataclass, field, asdict
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Optionally pull authoritative keys/colour maps from the project's shared
# module. We never *require* it and we never overwrite it; we only read from it
# so the new scripts agree with the existing four analysis scripts.
# ---------------------------------------------------------------------------
try:
    import niche_common as _nc  # type: ignore
except Exception:  # pragma: no cover - the module is project-local and optional
    _nc = None


def _from_niche_common(name: str, default):
    """Return getattr(niche_common, name) if it exists and is truthy, else default."""
    if _nc is not None:
        val = getattr(_nc, name, None)
        if val is not None:
            return val
    return default


# ---------------------------------------------------------------------------
# Candidate names used for auto-detection when neither the CLI nor niche_common
# pins a key. Order = priority.
# ---------------------------------------------------------------------------
_CELLTYPE_CANDIDATES = [
    "predicted_cell_type", "consensus_label", "cell_type", "celltype",
    "predicted.id", "CellType", "annotation",
]
_SAMPLE_CANDIDATES = [
    "sample", "sample_id", "Sample", "library_id", "batch", "orig.ident", "donor",
]
_CONDITION_CANDIDATES = [
    "condition", "Condition", "group", "Group", "status", "senescence",
]
_NICHE_CANDIDATES = [
    "niche", "nichecompass_niche", "niche_label", "leiden_nichecompass",
    "nichecompass_latent_cluster", "leiden",
]
_SPATIAL_OBSM_CANDIDATES = ["X_spatial", "spatial", "X_umap_spatial"]
_COUNTS_LAYER_CANDIDATES = ["counts", "raw_counts", "X_counts", "spliced"]


# ---------------------------------------------------------------------------
# Configuration object
# ---------------------------------------------------------------------------
@dataclass
class CCConfig:
    """Resolved keys + run parameters shared by both scripts."""
    # obs / obsm / layer keys
    celltype_key: str = field(default_factory=lambda: _from_niche_common("CELLTYPE_KEY", "predicted_cell_type"))
    groupby_key: Optional[str] = None          # what LIANA groups on; defaults to celltype_key
    sample_key: str = field(default_factory=lambda: _from_niche_common("SAMPLE_KEY", "sample"))
    condition_key: str = field(default_factory=lambda: _from_niche_common("CONDITION_KEY", "condition"))
    niche_key: Optional[str] = field(default_factory=lambda: _from_niche_common("NICHE_KEY", None))
    spatial_key: str = field(default_factory=lambda: _from_niche_common("SPATIAL_KEY", "X_spatial"))
    counts_layer: Optional[str] = None         # None -> auto-resolve

    # condition levels (CNTRL = reference/baseline, TEST = perturbed)
    cntrl: str = "CNTRL"
    test: str = "TEST"

    # resource / signature
    resource_name: str = "consensus"           # LIANA L-R resource
    genes_path: Optional[str] = None           # condition_genes.txt

    # spatial constraint
    spatial_mode: str = "colocalization"        # one of: colocalization | within_niche | none
    n_neighs: int = 6                          # k for spatial kNN graph
    coloc_metric: str = "contact_fraction"      # contact_fraction | nhood_zscore
    coloc_threshold: float = 0.05               # min contact fraction (or 0.0 for z>0)
    coloc_min_samples_frac: float = 1 / 3.0     # keep a pair if co-localized in >= this frac of samples

    # focus sets (interpretation only; the tensor is always built on the full matrix)
    sender_celltypes: list = field(default_factory=list)
    receiver_celltypes: list = field(default_factory=list)

    # optional: split a cell type into condition-signature +/- compartments
    split_by_signature: bool = False
    signature_score_key: str = "condition_score"
    signature_pos_quantile: float = 0.75        # cells above this quantile -> "|pos"

    # io
    out_dir: str = "communication_results"
    random_state: int = 1337

    def resolve_groupby(self) -> str:
        return self.groupby_key or self.celltype_key

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, default=str)


# ---------------------------------------------------------------------------
# Argument plumbing shared by both scripts
# ---------------------------------------------------------------------------
def add_common_args(parser):
    """Attach the arguments shared by both communication scripts."""
    g = parser.add_argument_group("data + keys")
    g.add_argument("--h5ad", required=True,
                   help="Path to the integrated AnnData (e.g. "
                        "nichecompass_results/objects/nichecompass_integrated.h5ad).")
    g.add_argument("--genes", dest="genes_path", default=None,
                   help="Path to condition_genes.txt (one gene symbol per line). "
                        "Generic signature; for senescence point it at your SASP/senescence list.")
    g.add_argument("--celltype-key", default=None, help="obs column with cell-type labels.")
    g.add_argument("--groupby", dest="groupby_key", default=None,
                   help="obs column LIANA groups on (defaults to --celltype-key; "
                        "point at a fine sub-state column if you have one).")
    g.add_argument("--sample-key", default=None, help="obs column identifying each slide/sample.")
    g.add_argument("--condition-key", default=None, help="obs column with the CNTRL/TEST grouping.")
    g.add_argument("--niche-key", default=None, help="obs column with niche labels (for within-niche mode).")
    g.add_argument("--spatial-key", default=None, help="obsm key with spatial coords (default X_spatial).")
    g.add_argument("--counts-layer", default=None,
                   help="layer holding raw integer counts. Auto-resolved if omitted.")

    g2 = parser.add_argument_group("condition levels")
    g2.add_argument("--cntrl", required=True,
                    help="Value in --condition-key used as the CONTROL/baseline level "
                         "(e.g. LS).")
    g2.add_argument("--test", required=True,
                    help="Value in --condition-key used as the TEST/perturbed level "
                         "(e.g. HS). log2FC is reported as TEST vs CNTRL.")

    g3 = parser.add_argument_group("output")
    g3.add_argument("--out-dir", default="communication_results",
                    help="Directory for tables and figures.")
    g3.add_argument("--random-state", type=int, default=1337)


def config_from_args(args) -> CCConfig:
    cfg = CCConfig()
    # CLI overrides take precedence over the dataclass defaults (which already
    # consulted niche_common).
    for attr in ["celltype_key", "groupby_key", "sample_key", "condition_key",
                 "niche_key", "spatial_key", "counts_layer", "out_dir", "random_state"]:
        v = getattr(args, attr, None)
        if v is not None:
            setattr(cfg, attr, v)
    if getattr(args, "genes_path", None):
        cfg.genes_path = args.genes_path
    cfg.cntrl = args.cntrl
    cfg.test = args.test
    return cfg


# ---------------------------------------------------------------------------
# AnnData loading + key resolution against the actual object
# ---------------------------------------------------------------------------
def load_adata(path: str):
    import anndata as ad
    if not os.path.exists(path):
        raise FileNotFoundError(f"AnnData not found: {path}")
    adata = ad.read_h5ad(path)
    return adata


def _resolve_obs_key(adata, current: Optional[str], candidates: Sequence[str], label: str,
                     required: bool = True) -> Optional[str]:
    cols = list(adata.obs.columns)
    # 1) explicit / niche_common value, if present on the object
    if current is not None and current in cols:
        return current
    # 2) auto-detect
    for c in candidates:
        if c in cols:
            return c
    if required:
        raise KeyError(
            f"Could not resolve the {label} obs key. Tried '{current}' and "
            f"{candidates}. Available obs columns: {cols}. "
            f"Pass it explicitly with the corresponding --*-key argument."
        )
    return None


def resolve_keys(adata, cfg: CCConfig, require_niche: bool = False) -> CCConfig:
    """Validate/auto-detect every key against the actual AnnData and update cfg in place."""
    cfg.celltype_key = _resolve_obs_key(adata, cfg.celltype_key, _CELLTYPE_CANDIDATES, "cell-type")
    if cfg.groupby_key is None:
        cfg.groupby_key = cfg.celltype_key
    else:
        cfg.groupby_key = _resolve_obs_key(adata, cfg.groupby_key, [cfg.groupby_key], "groupby")
    cfg.sample_key = _resolve_obs_key(adata, cfg.sample_key, _SAMPLE_CANDIDATES, "sample")
    cfg.condition_key = _resolve_obs_key(adata, cfg.condition_key, _CONDITION_CANDIDATES, "condition")
    cfg.niche_key = _resolve_obs_key(adata, cfg.niche_key, _NICHE_CANDIDATES, "niche",
                                     required=require_niche)

    # spatial obsm
    if cfg.spatial_key not in adata.obsm:
        found = next((k for k in _SPATIAL_OBSM_CANDIDATES if k in adata.obsm), None)
        if found is None:
            raise KeyError(
                f"No spatial coordinates in obsm. Tried '{cfg.spatial_key}' and "
                f"{_SPATIAL_OBSM_CANDIDATES}. Available obsm keys: {list(adata.obsm)}."
            )
        cfg.spatial_key = found

    # condition levels present?
    cond_vals = set(map(str, pd.unique(adata.obs[cfg.condition_key])))
    missing = [lvl for lvl in (cfg.cntrl, cfg.test) if str(lvl) not in cond_vals]
    if missing:
        raise ValueError(
            f"Condition level(s) {missing} not found in obs['{cfg.condition_key}']. "
            f"Observed levels: {sorted(cond_vals)}."
        )
    return cfg


# ---------------------------------------------------------------------------
# Counts handling: LIANA wants log-normalised expression; DE wants raw integers.
# ---------------------------------------------------------------------------
def _looks_like_counts(mat) -> bool:
    import scipy.sparse as sp
    sub = mat[:2000] if mat.shape[0] > 2000 else mat
    if sp.issparse(sub):
        data = sub.data
    else:
        data = np.asarray(sub).ravel()
    if data.size == 0:
        return False
    return bool(np.allclose(data, np.round(data)) and data.min() >= 0)


def resolve_counts_layer(adata, cfg: CCConfig) -> str:
    """Return the name of a layer holding raw integer counts, creating one in '.layers'
    from a counts-like .X or adata.raw if needed. Never mutates expression values."""
    # explicit / candidate layers
    cand = [cfg.counts_layer] if cfg.counts_layer else []
    cand += [c for c in _COUNTS_LAYER_CANDIDATES if c not in cand]
    for name in cand:
        if name and name in adata.layers and _looks_like_counts(adata.layers[name]):
            return name
    # .X integer?
    if _looks_like_counts(adata.X):
        adata.layers["counts"] = adata.X.copy()
        return "counts"
    # adata.raw integer?
    if adata.raw is not None and _looks_like_counts(adata.raw.X):
        raw = adata.raw.to_adata()
        # align to current var order where possible
        common = [g for g in adata.var_names if g in set(raw.var_names)]
        adata.layers["counts"] = raw[:, common].X.copy() if len(common) == adata.n_vars else raw.X.copy()
        return "counts"
    raise ValueError(
        "Could not find raw integer counts in a layer, .X, or .raw. "
        "Pseudobulk DE and LIANA normalisation both need them. "
        "Pass --counts-layer or ensure the conversion wrote integer counts."
    )


def make_lognorm_view(adata, cfg: CCConfig, counts_layer: str):
    """Return a *copy* with log1p(CPM-like) normalised expression in .X, suitable for LIANA.
    Raw counts are preserved in .layers['counts']."""
    import scanpy as sc
    work = adata.copy()
    work.X = work.layers[counts_layer].copy()
    work.layers["counts"] = work.layers[counts_layer].copy()
    sc.pp.normalize_total(work)      # default target = median library size
    sc.pp.log1p(work)
    return work


# ---------------------------------------------------------------------------
# Signature (condition_genes) + ligand/receptor universe
# ---------------------------------------------------------------------------
def load_condition_genes(path: Optional[str]) -> list:
    if not path:
        return []
    with open(path) as fh:
        genes = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    # de-dup, preserve order
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def liana_resource(resource_name: str = "consensus") -> pd.DataFrame:
    """Return the LIANA ligand-receptor resource as a DataFrame
    (columns include ligand / receptor complexes)."""
    import liana as li
    res = li.rs.select_resource(resource_name)
    # standardise column names across liana versions
    cols = {c.lower(): c for c in res.columns}
    return res.rename(columns={cols.get("ligand", "ligand"): "ligand",
                               cols.get("receptor", "receptor"): "receptor"})


def split_complex(token: str) -> list:
    """A LIANA complex like 'IL6_IL6R' or 'ITGA1&ITGB1' -> its subunits."""
    if not isinstance(token, str):
        return []
    for sep in ("_", "&", "+"):
        if sep in token:
            return [t for t in token.split(sep) if t]
    return [token]


def ligand_receptor_universe(resource_name: str = "consensus"):
    """Return (ligand_set, receptor_set) of individual gene symbols from the resource."""
    res = liana_resource(resource_name)
    ligands, receptors = set(), set()
    for v in res["ligand"].astype(str):
        ligands.update(split_complex(v))
    for v in res["receptor"].astype(str):
        receptors.update(split_complex(v))
    return ligands, receptors


# ---------------------------------------------------------------------------
# Spatial graph + co-localisation (the spatial constraint)
# ---------------------------------------------------------------------------
def build_spatial_graph(adata_sample, cfg: CCConfig):
    """Build a per-sample spatial kNN graph in-place (squidpy)."""
    import squidpy as sq
    # ensure squidpy finds coordinates under the conventional key
    if "spatial" not in adata_sample.obsm:
        adata_sample.obsm["spatial"] = np.asarray(adata_sample.obsm[cfg.spatial_key])
    sq.gr.spatial_neighbors(adata_sample, coord_type="generic", n_neighs=cfg.n_neighs)
    return adata_sample


def colocalization_table(adata, cfg: CCConfig) -> pd.DataFrame:
    """For every (sample, source, target) cell-type pair, compute a spatial
    co-localisation statistic and a boolean 'colocalized' flag.

    contact_fraction: fraction of target cells that have >=1 source-type cell in
                      their spatial kNN neighbourhood. Directional, density-aware,
                      no permutation cost.
    nhood_zscore:     squidpy neighbourhood-enrichment z-score (symmetric).
    """
    import scipy.sparse as sp
    groupby = cfg.resolve_groupby()
    rows = []
    for samp, idx in adata.obs.groupby(cfg.sample_key).groups.items():
        sub = adata[list(idx)].copy()
        if sub.n_obs < 10:
            continue
        build_spatial_graph(sub, cfg)
        A = sub.obsp["spatial_connectivities"]
        if not sp.issparse(A):
            A = sp.csr_matrix(A)
        A = (A > 0).astype(int)
        labels = sub.obs[groupby].astype(str).values
        cats = pd.unique(labels)
        # one-hot of source membership: neighbours_of_source = A @ onehot
        onehot = pd.get_dummies(pd.Series(labels, index=sub.obs_names))
        # neighbour count of each cell to each source type
        nbr_counts = A @ onehot.values            # (n_cells, n_types)
        has_src = (nbr_counts > 0).astype(int)    # cell has >=1 neighbour of that type
        nbr_df = pd.DataFrame(has_src, index=sub.obs_names, columns=onehot.columns)

        if cfg.coloc_metric == "nhood_zscore":
            import squidpy as sq
            sub.obs["_grp"] = pd.Categorical(labels)
            sq.gr.nhood_enrichment(sub, cluster_key="_grp", seed=cfg.random_state)
            z = sub.uns["_grp_nhood_enrichment"]["zscore"]
            zcats = list(sub.obs["_grp"].cat.categories)
            zdf = pd.DataFrame(z, index=zcats, columns=zcats)

        for tgt in cats:
            tgt_mask = labels == tgt
            n_tgt = int(tgt_mask.sum())
            if n_tgt == 0:
                continue
            for src in cats:
                if src not in nbr_df.columns:
                    frac = 0.0
                else:
                    frac = float(nbr_df.loc[tgt_mask, src].mean())
                stat = frac
                if cfg.coloc_metric == "nhood_zscore":
                    stat = float(zdf.loc[src, tgt]) if (src in zdf.index and tgt in zdf.columns) else np.nan
                rows.append({
                    "sample": samp, "source": src, "target": tgt,
                    "contact_fraction": frac, "nhood_zscore": (stat if cfg.coloc_metric == "nhood_zscore" else np.nan),
                    "n_target": n_tgt,
                })
    tab = pd.DataFrame(rows)
    if tab.empty:
        return tab
    if cfg.coloc_metric == "nhood_zscore":
        tab["colocalized"] = tab["nhood_zscore"] > (cfg.coloc_threshold if cfg.coloc_threshold else 0.0)
    else:
        tab["colocalized"] = tab["contact_fraction"] >= cfg.coloc_threshold
    return tab


def stable_colocalized_pairs(coloc: pd.DataFrame, cfg: CCConfig) -> set:
    """Return the set of (source, target) pairs co-localised in >= min_samples_frac of samples."""
    if coloc.empty:
        return set()
    n_samp = coloc["sample"].nunique()
    by_pair = (coloc.groupby(["source", "target"])["colocalized"].sum() / n_samp)
    keep = by_pair[by_pair >= cfg.coloc_min_samples_frac].index
    return set(map(tuple, keep))


# ---------------------------------------------------------------------------
# Optional: split a cell type into signature-positive / -negative compartments
# (generalises "senescent cells as senders")
# ---------------------------------------------------------------------------
def add_signature_split(adata, cfg: CCConfig, genes: Sequence[str]) -> str:
    """Score `genes` per cell, flag the top quantile as '|pos', append to the
    groupby labels, and return the name of the new grouping column."""
    import scanpy as sc
    present = [g for g in genes if g in adata.var_names]
    if not present:
        raise ValueError("None of the condition genes are present; cannot split by signature.")
    sc.tl.score_genes(adata, gene_list=present, score_name=cfg.signature_score_key, use_raw=False)
    groupby = cfg.resolve_groupby()
    out_col = f"{groupby}__sig"
    thr = adata.obs[cfg.signature_score_key].quantile(cfg.signature_pos_quantile)
    pos = adata.obs[cfg.signature_score_key] >= thr
    suffix = np.where(pos, "|pos", "|neg")
    adata.obs[out_col] = (adata.obs[groupby].astype(str) + suffix).astype("category")
    return out_col


# ---------------------------------------------------------------------------
# Pseudobulk (the correct n=3 substrate for DE and for slide-as-unit summaries)
# ---------------------------------------------------------------------------
def pseudobulk_matrix(adata, cfg: CCConfig, counts_layer: str,
                      group_keys: Sequence[str]) -> tuple:
    """Sum raw counts within each combination of group_keys.

    Returns (counts_df [genes x pseudobulk], coldata [pseudobulk x metadata]).
    Pseudobulk id = '||'.join(group values).
    """
    import scipy.sparse as sp
    X = adata.layers[counts_layer]
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    obs = adata.obs[list(group_keys)].astype(str)
    obs["_pb"] = obs[list(group_keys)].agg("||".join, axis=1)
    groups = obs.groupby("_pb").indices

    mats, ids, meta = [], [], []
    for pb, idx in groups.items():
        mats.append(np.asarray(X[idx].sum(axis=0)).ravel())
        ids.append(pb)
        row = {k: obs.iloc[idx[0]][k] for k in group_keys}
        row["_n_cells"] = len(idx)
        meta.append(row)
    counts_df = pd.DataFrame(np.vstack(mats), index=ids, columns=adata.var_names).T  # genes x pb
    coldata = pd.DataFrame(meta, index=ids)
    return counts_df, coldata


# ---------------------------------------------------------------------------
# Statistics for n=3 vs n=3 (honest, slide-as-unit)
# ---------------------------------------------------------------------------
def cliffs_delta(a: Sequence[float], b: Sequence[float]) -> float:
    """Cliff's delta effect size (a vs b); +1 means a stochastically dominates b."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return np.nan
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (a.size * b.size)


def mannwhitney_test(values_cntrl: Sequence[float], values_test: Sequence[float]) -> dict:
    """Two-sided exact Mann-Whitney, TEST vs CNTRL, with effect sizes.
    With n=3 vs n=3 the minimum two-sided p is ~0.1; effect size carries the signal."""
    from scipy.stats import mannwhitneyu
    a = np.asarray(values_test, float)
    b = np.asarray(values_cntrl, float)
    out = {"n_test": a.size, "n_cntrl": b.size, "median_test": np.nanmedian(a),
           "median_cntrl": np.nanmedian(b), "cliffs_delta": cliffs_delta(a, b)}
    if a.size < 1 or b.size < 1 or (np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]):
        out["U"], out["pvalue"] = np.nan, np.nan
        return out
    try:
        U, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        U, p = np.nan, np.nan
    out["U"], out["pvalue"] = U, p
    return out


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaNs pass through."""
    p = np.asarray(pvals, float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    n = ranked.size
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def pval_stars(p) -> str:
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "ns"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Plotting theme (kept minimal; both scripts share it)
# ---------------------------------------------------------------------------
def set_plot_theme():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    return plt


def condition_palette(cfg: CCConfig) -> dict:
    """A muted, consistent CNTRL/TEST palette (TEST warm, CNTRL cool)."""
    pal = _from_niche_common("CONDITION_PALETTE", None)
    if isinstance(pal, dict) and cfg.cntrl in pal and cfg.test in pal:
        return pal
    return {cfg.cntrl: "#4C72B0", cfg.test: "#C44E52"}


# ---------------------------------------------------------------------------
# IO + reporting
# ---------------------------------------------------------------------------
def ensure_out(cfg: CCConfig, *subdirs) -> str:
    path = os.path.join(cfg.out_dir, *subdirs) if subdirs else cfg.out_dir
    os.makedirs(path, exist_ok=True)
    return path


CAVEATS = textwrap.dedent("""\
    Interpretation caveats (carry these into any figure legend / talk):

    1. PROVENANCE / CIRCULARITY. If the CNTRL/TEST split was defined from the same
       signature you are now testing (e.g. a senescence score on this dataset),
       then 'TEST has more signal' is circular. Report where/in which cells the
       signal concentrates, not that TEST exceeds CNTRL, unless the split came
       from an independent assay.
    2. n = 3 vs 3. The slide is the unit. Two-sided Mann-Whitney bottoms out at
       p ~= 0.1 for 3 vs 3, so nothing survives FDR on p alone. Lead with effect
       size (Cliff's delta, log2FC, factor-loading separation); treat p as
       supporting, not gating.
    3. DONOR CONFOUND. With 3 donors per arm, the condition is confounded with any
       donor-level covariate (age, comorbidity, wound chronicity, slide source).
       Per-sample pseudobulk does not remove a group-level confound.
    4. CONTACT != COMMUNICATION. The spatial constraint removes pairs that never
       touch; it does not prove signalling. Co-localisation is necessary, not
       sufficient.
    5. DIRECTION CONVENTION. log2FC and 'enriched in' are always TEST relative to
       CNTRL. Positive log2FC = higher in TEST.
    """)


def write_report(cfg: CCConfig, extra: dict, filename: str):
    out = ensure_out(cfg)
    payload = {"config": json.loads(cfg.to_json()), **extra}
    with open(os.path.join(out, filename), "w") as fh:
        json.dump(payload, fh, indent=2, default=str)


def banner(msg: str):
    print("\n" + "=" * 78 + f"\n{msg}\n" + "=" * 78, flush=True)


# ---------------------------------------------------------------------------
# Self-test: print the resolved configuration (no heavy work).
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    cfg = CCConfig(cntrl="LS", test="HS")
    print("niche_common importable:", _nc is not None)
    print(cfg.to_json())
    print(CAVEATS)
