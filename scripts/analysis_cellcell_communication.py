#!/usr/bin/env python
"""
analysis_cellcell_communication.py
==================================
Cell-cell communication layer (the directional sender -> receiver ligand-receptor
analysis; the Figure-7 analog of the reference paper).

WHAT IT DOES
------------
1. Loads the integrated AnnData, log-normalises a working copy for LIANA while
   keeping raw counts in a layer.
2. Runs LIANA+ consensus (rank_aggregate) ligand-receptor scoring PER SAMPLE,
   grouped by cell type (or a finer sub-state column you point it at).
3. SPATIALLY CONSTRAINS the inference: builds a per-sample spatial kNN graph and
   keeps only ligand-receptor calls between sender/receiver cell types that are
   actually spatially co-localised (contact-fraction or neighbourhood-enrichment),
   so you never infer contact between cells that never touch. Three modes:
       colocalization (default) | within_niche | none
4. Builds a 4D communication tensor (samples x LR x sender x receiver) from the
   spatially-filtered LIANA output and decomposes it with Tensor-cell2cell to get
   context-driven communication factors.
5. Compares the factor CONTEXT loadings between CNTRL and TEST (slide-as-unit
   Mann-Whitney + Cliff's delta + BH-FDR), i.e. which communication programmes
   separate the two conditions.
6. Produces a focused sender -> receiver view (you choose which cell types are
   senders and which are receivers) and writes the long-form ligand-receptor
   table that the factor-role script consumes downstream.

Everything is phrased as CNTRL vs TEST and driven by a generic gene list
(condition_genes.txt) so the script is reusable across datasets/conditions.

EXAMPLE
-------
    python analysis_cellcell_communication.py \
        --h5ad nichecompass_results/objects/nichecompass_integrated.h5ad \
        --genes condition_genes.txt \
        --cntrl L --test D \
        --sender-celltypes "Fibroblast" "Macrophage" "Endothelial" \
        --receiver-celltypes "Keratinocyte" "Fibroblast" "Endothelial" \
        --spatial-mode colocalization \
        --split-by-signature \
        --out-dir communication_results

DEPENDENCIES
------------
    scanpy anndata squidpy liana cell2cell  (+ numpy pandas matplotlib scipy)
Tensor-cell2cell factorisation can use a GPU (--device cuda) but runs on CPU.
APIs of liana / cell2cell are version-sensitive; pinned calls are noted inline.
"""

from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd

import cc_common as cc


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    cc.add_common_args(p)

    g = p.add_argument_group("LIANA")
    g.add_argument("--resource", default="consensus",
                   help="LIANA ligand-receptor resource (default: consensus). Human symbols.")
    g.add_argument("--score-key", default="magnitude_rank",
                   help="LIANA score used to build the tensor (default: magnitude_rank).")
    g.add_argument("--n-perms", type=int, default=100,
                   help="LIANA permutations for specificity (0/None to skip for speed).")
    g.add_argument("--min-cells-per-group", type=int, default=10,
                   help="Drop sender/receiver groups with fewer cells than this, per sample.")

    g2 = p.add_argument_group("spatial constraint")
    g2.add_argument("--spatial-mode", default="colocalization",
                    choices=["colocalization", "within_niche", "none"],
                    help="How to spatially constrain L-R inference.")
    g2.add_argument("--n-neighs", type=int, default=6, help="k for the spatial kNN graph.")
    g2.add_argument("--coloc-metric", default="contact_fraction",
                    choices=["contact_fraction", "nhood_zscore"])
    g2.add_argument("--coloc-threshold", type=float, default=0.05,
                    help="Min contact fraction (or z>this) to call a pair co-localised.")
    g2.add_argument("--coloc-min-samples-frac", type=float, default=1 / 3.0,
                    help="Keep a pair if co-localised in at least this fraction of samples.")

    g3 = p.add_argument_group("focus + signature split")
    g3.add_argument("--sender-celltypes", nargs="*", default=[],
                    help="Cell-type labels to treat as senders in the focused view.")
    g3.add_argument("--receiver-celltypes", nargs="*", default=[],
                    help="Cell-type labels to treat as receivers in the focused view.")
    g3.add_argument("--split-by-signature", action="store_true",
                    help="Also split each group into signature +/- compartments "
                         "(generalises 'senescent cells as senders').")
    g3.add_argument("--signature-pos-quantile", type=float, default=0.75)

    g4 = p.add_argument_group("tensor")
    g4.add_argument("--tensor-rank", type=int, default=None,
                    help="Number of factors. If omitted, an elbow rank-selection is run.")
    g4.add_argument("--max-rank", type=int, default=15, help="Upper rank for elbow search.")
    g4.add_argument("--elbow-runs", type=int, default=10,
                    help="Decomposition repeats per rank in the elbow search "
                         "(cell2cell 'regular'=10, 'robust'=20; higher is more stable, slower).")
    g4.add_argument("--tf-init", default="svd", choices=["svd", "random"])
    g4.add_argument("--device", default="cpu", help="cpu or cuda (GPU) for factorisation.")
    return p.parse_args()


# ---------------------------------------------------------------------------
def run_liana_by_sample(adata_ln, cfg: cc.CCConfig, args, groupby: str) -> pd.DataFrame:
    """Run LIANA rank_aggregate per sample; return the long-form result frame."""
    import liana as li
    cc.banner(f"LIANA rank_aggregate (resource={args.resource}, groupby={groupby})")

    # drop tiny groups per sample so LIANA does not score ghost cell types
    keep = np.ones(adata_ln.n_obs, dtype=bool)
    for samp, idx in adata_ln.obs.groupby(cfg.sample_key).groups.items():
        sub = adata_ln.obs.loc[idx]
        counts = sub[groupby].value_counts()
        small = set(counts[counts < args.min_cells_per_group].index)
        if small:
            mask = adata_ln.obs[cfg.sample_key].astype(str).isin([str(samp)]) & \
                   adata_ln.obs[groupby].astype(str).isin(map(str, small))
            keep &= ~mask.values
    adata_use = adata_ln[keep].copy()

    n_perms = None if (args.n_perms in (0, None)) else args.n_perms
    li.mt.rank_aggregate.by_sample(
        adata_use,
        groupby=groupby,
        resource_name=args.resource,   # human gene symbols
        sample_key=cfg.sample_key,
        use_raw=False,                 # expression is log-normalised in .X
        n_perms=n_perms,
        return_all_lrs=True,
        verbose=True,
    )
    res = adata_use.uns["liana_res"].copy()
    # normalise column names across liana versions
    rename = {}
    for a, b in [("ligand.complex", "ligand_complex"), ("receptor.complex", "receptor_complex")]:
        if a in res.columns and b not in res.columns:
            rename[a] = b
    res = res.rename(columns=rename)
    return res


def apply_spatial_constraint(liana_res: pd.DataFrame, adata, cfg: cc.CCConfig, args):
    """Filter the long-form LIANA frame to spatially-plausible (sample, source, target)
    triples. Returns (filtered_res, coloc_table)."""
    if args.spatial_mode == "none":
        cc.banner("Spatial constraint: NONE (using all sender/receiver pairs)")
        return liana_res, pd.DataFrame()

    if args.spatial_mode == "within_niche":
        cc.banner("Spatial constraint: WITHIN-NICHE")
        # Restrict each cell to its niche, recompute LIANA per (sample, niche) would be
        # ideal but expensive; instead we keep L-R calls only where sender & receiver
        # co-occur within the same niche in that sample (cheap, faithful proxy).
        groupby = cfg.resolve_groupby()
        allowed = set()
        for samp, idx in adata.obs.groupby(cfg.sample_key).groups.items():
            sub = adata.obs.loc[idx]
            for niche, nidx in sub.groupby(cfg.niche_key).groups.items():
                present = set(map(str, sub.loc[nidx, groupby].unique()))
                for s in present:
                    for t in present:
                        allowed.add((str(samp), s, t))
        before = len(liana_res)
        samp_col = cfg.sample_key if cfg.sample_key in liana_res.columns else "sample"
        keys = zip(liana_res[samp_col].astype(str),
                   liana_res["source"].astype(str),
                   liana_res["target"].astype(str))
        mask = pd.Series([k in allowed for k in keys], index=liana_res.index)
        filt = liana_res[mask]
        print(f"  kept {len(filt)}/{before} L-R rows co-occurring within a niche")
        return filt, pd.DataFrame()

    # default: colocalization
    cc.banner(f"Spatial constraint: CO-LOCALISATION ({args.coloc_metric}, "
              f"thr={args.coloc_threshold}, in >= {args.coloc_min_samples_frac:.2f} of samples)")
    coloc = cc.colocalization_table(adata, cfg)
    if coloc.empty:
        warnings.warn("Co-localisation table is empty; skipping spatial filter.")
        return liana_res, coloc

    # per-(sample, source, target) flag
    samp_col = cfg.sample_key if cfg.sample_key in liana_res.columns else "sample"
    coloc_key = coloc.set_index(["sample", "source", "target"])["colocalized"]
    keys = list(zip(liana_res[samp_col].astype(str),
                    liana_res["source"].astype(str),
                    liana_res["target"].astype(str)))
    flags = np.array([bool(coloc_key.get((s, src, tgt), False)) for s, src, tgt in keys])
    before = len(liana_res)
    filt = liana_res[flags]
    print(f"  kept {len(filt)}/{before} per-sample L-R rows between co-localised pairs")

    # also report which pairs are *stably* co-localised (for the focused view)
    stable = cc.stable_colocalized_pairs(coloc, cfg)
    print(f"  {len(stable)} sender->receiver pairs co-localised in "
          f">= {cfg.coloc_min_samples_frac:.2f} of samples")
    return filt, coloc


def build_and_decompose_tensor(liana_res: pd.DataFrame, cfg: cc.CCConfig, args,
                               sample_to_condition: dict):
    """Build the 4D tensor with liana.multi.to_tensor_c2c and decompose it with
    Tensor-cell2cell.
    Returns (tensor, factors_dict, context_df, tensor_metadata, rank)."""
    import liana as li
    import cell2cell as c2c

    cc.banner("Building 4D communication tensor (samples x LR x sender x receiver)")
    samp_col = cfg.sample_key if cfg.sample_key in liana_res.columns else "sample"
    build_kwargs = dict(
        liana_res=liana_res,
        sample_key=samp_col,
        source_key="source",
        target_key="target",
        ligand_key="ligand_complex",
        receptor_key="receptor_complex",
        score_key=args.score_key,          # e.g. magnitude_rank
        non_negative=True,                 # decomposition needs >= 0
        inverse_fun=lambda x: 1.0 - x,     # rank: lower is stronger -> 1-rank so higher is stronger
        how="outer",                       # union of cells & LRs across samples
        outer_fraction=cfg.coloc_min_samples_frac,
        lr_fill=np.nan, cell_fill=np.nan,
    )
    if args.device and args.device != "cpu":
        build_kwargs["device"] = args.device
    try:
        tensor = li.multi.to_tensor_c2c(**build_kwargs)
    except TypeError:
        # older/newer signature without a device kwarg -> build on CPU
        build_kwargs.pop("device", None)
        tensor = li.multi.to_tensor_c2c(**build_kwargs)
    print(f"  tensor shape (contexts, LR, senders, receivers): {tensor.tensor.shape}")

    # context metadata: map each sample/context -> condition for plotting/testing
    context_order = list(tensor.order_names[0])
    ctx_meta = {s: sample_to_condition.get(str(s), "NA") for s in context_order}
    meta = c2c.tensor.generate_tensor_metadata(
        interaction_tensor=tensor,
        metadata_dicts=[ctx_meta, None, None, None],
        fill_with_order_elements=True,
    )

    cc.banner("Tensor-cell2cell decomposition")
    if args.tensor_rank is None:
        print("  running elbow rank selection ...")
        # NB: tf_optimization is an R-wrapper concept, NOT a cell2cell kwarg.
        # The Python API controls robustness via `runs` (and tol/n_iter_max in
        # compute_tensor_factorization). automatic_elbow=True stores the chosen
        # rank on the tensor as `tensor.rank`.
        tensor.elbow_rank_selection(
            upper_rank=args.max_rank, runs=args.elbow_runs, init=args.tf_init,
            automatic_elbow=True, random_state=cfg.random_state,
        )
        rank = int(getattr(tensor, "rank", 0) or 0)
        if rank < 1:
            rank = max(2, min(8, args.max_rank))
            warnings.warn(f"Elbow analysis did not return a usable rank "
                          f"(tensor may be sparse); defaulting to rank={rank}. "
                          f"Consider re-running LIANA with a smaller expr_prop or "
                          f"raising --coloc-min-samples-frac.")
        print(f"  selected rank = {rank}")
    else:
        rank = int(args.tensor_rank)

    tensor.compute_tensor_factorization(
        rank=rank, init=args.tf_init, random_state=cfg.random_state,
    )

    # factors: dict with keys 'Contexts','Ligand-Receptor Pairs','Sender Cells','Receiver Cells'
    factors = {k: pd.DataFrame(v) for k, v in tensor.factors.items()}
    # context loadings -> attach condition
    ctx_key = next(k for k in factors if "ontext" in k or "Context" in k)
    context_df = factors[ctx_key].copy()
    context_df.index = [str(i) for i in context_df.index]
    context_df[cfg.condition_key] = [sample_to_condition.get(i, "NA") for i in context_df.index]
    return tensor, factors, context_df, meta, rank


def test_factors_by_condition(context_df: pd.DataFrame, cfg: cc.CCConfig) -> pd.DataFrame:
    """Mann-Whitney (TEST vs CNTRL) on each factor's context loadings."""
    cc.banner("Comparing factor context-loadings: TEST vs CNTRL")
    factor_cols = [c for c in context_df.columns if c != cfg.condition_key]
    cond = context_df[cfg.condition_key].astype(str)
    rows = []
    for fc in factor_cols:
        b = context_df.loc[cond == cfg.cntrl, fc].values
        a = context_df.loc[cond == cfg.test, fc].values
        res = cc.mannwhitney_test(b, a)
        res["factor"] = fc
        res["direction"] = ("higher_in_TEST" if res["median_test"] >= res["median_cntrl"]
                            else "higher_in_CNTRL")
        rows.append(res)
    df = pd.DataFrame(rows).set_index("factor")
    df["padj_BH"] = cc.bh_fdr(df["pvalue"].values)
    df["stars"] = df["padj_BH"].map(cc.pval_stars)
    df = df.sort_values("cliffs_delta", key=lambda s: s.abs(), ascending=False)
    return df[["n_cntrl", "n_test", "median_cntrl", "median_test", "cliffs_delta",
               "U", "pvalue", "padj_BH", "stars", "direction"]]


def focused_sender_receiver(liana_res: pd.DataFrame, cfg: cc.CCConfig, args,
                            sample_to_condition: dict) -> pd.DataFrame:
    """Aggregate the focused sender->receiver interactions and compare CNTRL vs TEST,
    slide-as-unit. Returns one row per (ligand, receptor, source, target)."""
    cc.banner("Focused sender -> receiver differential (slide-as-unit)")
    senders = set(map(str, args.sender_celltypes)) or None
    receivers = set(map(str, args.receiver_celltypes)) or None
    samp_col = cfg.sample_key if cfg.sample_key in liana_res.columns else "sample"

    df = liana_res.copy()
    # focused view matches either the raw label or the '|pos'/'|neg' split label prefix
    def _match(series, allowed):
        if allowed is None:
            return np.ones(len(series), dtype=bool)
        base = series.astype(str).str.split("|").str[0]
        return series.astype(str).isin(allowed) | base.isin(allowed)
    df = df[_match(df["source"], senders) & _match(df["target"], receivers)]
    if df.empty:
        warnings.warn("No focused sender->receiver interactions found.")
        return df

    # higher score = stronger (we already use magnitude_rank; convert to strength)
    score = args.score_key
    df["_strength"] = 1.0 - df[score] if "rank" in score else df[score]
    df["_cond"] = df[samp_col].astype(str).map(lambda s: sample_to_condition.get(s, "NA"))

    # per-(LR, source, target) per-sample mean strength, then MW across slides
    grp_cols = ["ligand_complex", "receptor_complex", "source", "target"]
    per_sample = (df.groupby(grp_cols + [samp_col, "_cond"])["_strength"]
                    .mean().reset_index())
    rows = []
    for keys, sub in per_sample.groupby(grp_cols):
        b = sub.loc[sub["_cond"] == cfg.cntrl, "_strength"].values
        a = sub.loc[sub["_cond"] == cfg.test, "_strength"].values
        res = cc.mannwhitney_test(b, a)
        rec = dict(zip(grp_cols, keys))
        rec.update({
            "mean_strength_cntrl": np.nanmean(b) if b.size else np.nan,
            "mean_strength_test": np.nanmean(a) if a.size else np.nan,
            "delta_test_minus_cntrl": (np.nanmean(a) if a.size else np.nan) -
                                       (np.nanmean(b) if b.size else np.nan),
            "cliffs_delta": res["cliffs_delta"], "pvalue": res["pvalue"],
        })
        rows.append(rec)
    out = pd.DataFrame(rows)
    if not out.empty:
        out["padj_BH"] = cc.bh_fdr(out["pvalue"].values)
        out["enriched_in"] = np.where(out["delta_test_minus_cntrl"] >= 0, cfg.test, cfg.cntrl)
        out = out.sort_values("delta_test_minus_cntrl", key=lambda s: s.abs(), ascending=False)
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_factor_heatmaps(factors: dict, out_dir: str, top_n: int = 15):
    plt = cc.set_plot_theme()
    import matplotlib.pyplot as _plt  # noqa
    for name, df in factors.items():
        if df.shape[0] > 60:
            # keep the top loaders per factor for readability
            keep = set()
            for col in df.columns:
                keep.update(df[col].abs().sort_values(ascending=False).head(top_n).index)
            d = df.loc[sorted(keep)]
        else:
            d = df
        fig, ax = plt.subplots(figsize=(1.2 + 0.5 * d.shape[1], 2 + 0.18 * d.shape[0]))
        im = ax.imshow(d.values, aspect="auto", cmap="magma")
        ax.set_xticks(range(d.shape[1])); ax.set_xticklabels(d.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(d.shape[0])); ax.set_yticklabels(d.index, fontsize=5)
        ax.set_title(f"Tensor factors: {name}")
        fig.colorbar(im, ax=ax, shrink=0.6, label="loading")
        safe = name.replace(" ", "_").replace("/", "-")
        fig.savefig(os.path.join(out_dir, f"tensor_factor_{safe}.pdf"))
        plt.close(fig)


def plot_context_boxplots(context_df: pd.DataFrame, factor_stats: pd.DataFrame,
                          cfg: cc.CCConfig, out_dir: str):
    plt = cc.set_plot_theme()
    pal = cc.condition_palette(cfg)
    factor_cols = [c for c in context_df.columns if c != cfg.condition_key]
    n = len(factor_cols)
    ncol = min(4, n); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.4 * ncol, 2.4 * nrow), squeeze=False)
    cond = context_df[cfg.condition_key].astype(str)
    for i, fc in enumerate(factor_cols):
        ax = axes[i // ncol][i % ncol]
        for j, lvl in enumerate([cfg.cntrl, cfg.test]):
            vals = context_df.loc[cond == lvl, fc].values
            ax.scatter(np.full(vals.size, j) + np.random.uniform(-0.07, 0.07, vals.size),
                       vals, color=pal.get(lvl, "gray"), s=18, zorder=3)
            ax.hlines(np.median(vals) if vals.size else 0, j - 0.2, j + 0.2,
                      color="k", lw=1.5, zorder=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels([cfg.cntrl, cfg.test])
        star = factor_stats.loc[fc, "stars"] if fc in factor_stats.index else ""
        cd = factor_stats.loc[fc, "cliffs_delta"] if fc in factor_stats.index else np.nan
        ax.set_title(f"{fc}\n{star}  (delta={cd:.2f})", fontsize=8)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle("Communication-factor context loadings by condition", y=1.0)
    fig.savefig(os.path.join(out_dir, "tensor_context_loadings_by_condition.pdf"))
    plt.close(fig)


def plot_focused_dotplot(focused: pd.DataFrame, cfg: cc.CCConfig, out_dir: str, top: int = 30):
    if focused.empty:
        return
    plt = cc.set_plot_theme()
    d = focused.head(top).copy()
    d["lr"] = d["ligand_complex"].astype(str) + " -> " + d["receptor_complex"].astype(str)
    d["sr"] = d["source"].astype(str) + "\n->" + d["target"].astype(str)
    lrs = list(dict.fromkeys(d["lr"]))
    srs = list(dict.fromkeys(d["sr"]))
    fig, ax = plt.subplots(figsize=(1.5 + 0.7 * len(srs), 1.5 + 0.32 * len(lrs)))
    for _, r in d.iterrows():
        x = srs.index(r["sr"]); y = lrs.index(r["lr"])
        size = 20 + 260 * (abs(r["delta_test_minus_cntrl"]) /
                           (d["delta_test_minus_cntrl"].abs().max() + 1e-9))
        color = cc.condition_palette(cfg).get(r["enriched_in"], "gray")
        ax.scatter(x, y, s=size, color=color, edgecolor="k", lw=0.3)
    ax.set_xticks(range(len(srs))); ax.set_xticklabels(srs, rotation=0, fontsize=6)
    ax.set_yticks(range(len(lrs))); ax.set_yticklabels(lrs, fontsize=6)
    ax.set_title("Focused sender -> receiver L-R\n(size=|delta|, colour=enriched-in)")
    fig.savefig(os.path.join(out_dir, "focused_sender_receiver_dotplot.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    cfg = cc.config_from_args(args)
    # pull spatial/tensor params into cfg
    cfg.resource_name = args.resource
    cfg.spatial_mode = args.spatial_mode
    cfg.n_neighs = args.n_neighs
    cfg.coloc_metric = args.coloc_metric
    cfg.coloc_threshold = args.coloc_threshold
    cfg.coloc_min_samples_frac = args.coloc_min_samples_frac
    cfg.sender_celltypes = args.sender_celltypes
    cfg.receiver_celltypes = args.receiver_celltypes
    cfg.split_by_signature = args.split_by_signature
    cfg.signature_pos_quantile = args.signature_pos_quantile

    out_dir = cc.ensure_out(cfg)
    cc.banner("LOADING DATA")
    adata = cc.load_adata(args.h5ad)
    cfg = cc.resolve_keys(adata, cfg, require_niche=(args.spatial_mode == "within_niche"))
    counts_layer = cc.resolve_counts_layer(adata, cfg)
    print(f"  cells={adata.n_obs}  genes={adata.n_vars}")
    print(f"  keys: celltype={cfg.celltype_key} sample={cfg.sample_key} "
          f"condition={cfg.condition_key} niche={cfg.niche_key} spatial={cfg.spatial_key}")
    print(f"  counts layer = {counts_layer};  CNTRL={cfg.cntrl}  TEST={cfg.test}")

    genes = cc.load_condition_genes(args.genes_path)
    print(f"  condition genes loaded: {len(genes)}")

    # log-normalised working copy for LIANA (raw counts preserved in layer)
    adata_ln = cc.make_lognorm_view(adata, cfg, counts_layer)

    # optionally split groups by the signature (e.g. senescent vs not)
    groupby = cfg.resolve_groupby()
    if cfg.split_by_signature and genes:
        groupby = cc.add_signature_split(adata_ln, cfg, genes)
        # mirror the split column onto the raw object too (used by spatial table)
        adata.obs[groupby] = adata_ln.obs[groupby].values
        print(f"  signature split -> grouping on '{groupby}'")

    # sample -> condition map
    s2c = (adata.obs[[cfg.sample_key, cfg.condition_key]].astype(str)
           .drop_duplicates().set_index(cfg.sample_key)[cfg.condition_key].to_dict())

    # 1) LIANA per sample
    liana_res = run_liana_by_sample(adata_ln, cfg, args, groupby)
    liana_res.to_csv(os.path.join(out_dir, "liana_res_per_sample_raw.csv"), index=False)

    # 2) spatial constraint
    liana_filt, coloc = apply_spatial_constraint(liana_res, adata, cfg, args)
    liana_filt.to_csv(os.path.join(out_dir, "liana_res_per_sample_spatial.csv"), index=False)
    if not coloc.empty:
        coloc.to_csv(os.path.join(out_dir, "spatial_colocalization_table.csv"), index=False)

    # 3) tensor + decomposition
    tensor, factors, context_df, meta, rank = build_and_decompose_tensor(
        liana_filt, cfg, args, s2c)
    for name, df in factors.items():
        safe = name.replace(" ", "_").replace("/", "-")
        df.to_csv(os.path.join(out_dir, f"tensor_factor_{safe}.csv"))
    context_df.to_csv(os.path.join(out_dir, "tensor_context_loadings.csv"))

    # 4) factor differential by condition
    factor_stats = test_factors_by_condition(context_df, cfg)
    factor_stats.to_csv(os.path.join(out_dir, "tensor_factor_condition_tests.csv"))
    print(factor_stats.to_string())

    # 5) focused sender -> receiver view (bridge to the factor-role script)
    focused = focused_sender_receiver(liana_filt, cfg, args, s2c)
    if not focused.empty:
        focused.to_csv(os.path.join(out_dir, "focused_sender_receiver_interactions.csv"),
                       index=False)

    # 6) figures
    cc.banner("WRITING FIGURES")
    fig_dir = cc.ensure_out(cfg, "figures")
    try:
        plot_factor_heatmaps(factors, fig_dir)
        plot_context_boxplots(context_df, factor_stats, cfg, fig_dir)
        plot_focused_dotplot(focused, cfg, fig_dir)
    except Exception as e:  # plotting must never sink the numerical outputs
        warnings.warn(f"Figure generation issue: {e}")

    # report
    cc.write_report(cfg, {
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "tensor_rank": int(rank), "tensor_shape": list(tensor.tensor.shape),
        "n_liana_rows_raw": int(len(liana_res)),
        "n_liana_rows_spatial": int(len(liana_filt)),
        "groupby": groupby,
        "n_focused_interactions": int(len(focused)),
        "caveats": cc.CAVEATS,
    }, "cellcell_communication_report.json")

    print(cc.CAVEATS)
    cc.banner("DONE: cell-cell communication")
    print(f"Outputs in: {out_dir}")
    print("Hand 'focused_sender_receiver_interactions.csv' and 'liana_res_per_sample_spatial.csv' "
          "to analysis_condition_factor_roles.py.")


if __name__ == "__main__":
    main()
