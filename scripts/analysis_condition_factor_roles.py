#!/usr/bin/env python
"""
analysis_condition_factor_roles.py
==================================
Signalling-factor-level dual-role analysis.

SASP-factor-level analysis:
it decomposes the composite signature (condition_genes.txt = your SASP, this can be metal/senescence
list) into individual secreted signalling factors and asks, per factor:

    (i)   is it enriched in CNTRL or TEST?         -> pseudobulk DE by SAMPLE
                                                       (pydeseq2 / DESeq2; the
                                                       correct n=3 engine)
    (ii)  which cell type secretes it?             -> per-cell-type expression +
                                                       where the DE is strongest
    (iii) which receptor on which receiver does    -> joined from the cell-cell
          it target?                                  communication step's
                                                       focused L-R table
    (iv)  does the literature class it pro- or      -> from an editable role table
          anti-(healing / the process of interest)?    (factor_roles_template.csv)

Deliverable: "pro-process factors enriched in CNTRL  vs  anti-process factors enriched in TEST"
i.e. a dual-role map where x = log2FC(TEST vs CNTRL) and colour = pro/anti role.

It is deliberately CNTRL/TEST + condition_genes driven so it is reusable for any
two-level contrast and any signature.

EXAMPLE
-------
    python analysis_condition_factor_roles.py \
        --h5ad nichecompass_results/objects/nichecompass_integrated.h5ad \
        --genes condition_genes.txt \
        --cntrl LS --test HS \
        --liana-csv communication_results/liana_res_per_sample_spatial.csv \
        --focused-csv communication_results/focused_sender_receiver_interactions.csv \
        --roles-csv factor_roles_template.csv \
        --out-dir factor_roles_results

DEPENDENCIES
------------
    scanpy anndata pydeseq2 liana  (+ numpy pandas matplotlib scipy)
If pydeseq2 is unavailable, pass --de-backend export to write the pseudobulk
matrices + a ready-to-run R DESeq2/edgeR script (run_pseudobulk_DE.R) instead.
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
    # the shared default out-dir is for the communication step; override it here
    p.set_defaults(out_dir="factor_roles_results")

    g = p.add_argument_group("inputs from the communication step")
    g.add_argument("--liana-csv", default=None,
                   help="liana_res_per_sample_spatial.csv from the communication step.")
    g.add_argument("--focused-csv", default=None,
                   help="focused_sender_receiver_interactions.csv from the communication step.")
    g.add_argument("--roles-csv", default=None,
                   help="Editable factor-role table (cols: factor, role[, note]). "
                        "role in {pro, anti, context}. See factor_roles_template.csv.")

    g2 = p.add_argument_group("DE (pseudobulk by sample)")
    g2.add_argument("--de-backend", default="pydeseq2",
                    choices=["pydeseq2", "export"],
                    help="pydeseq2 = run DESeq2 in Python; export = write matrices + R script.")
    g2.add_argument("--resource", default="consensus",
                    help="LIANA resource used to decide which signature genes are ligands.")
    g2.add_argument("--min-cells-pb", type=int, default=10,
                    help="Drop a (sample x cell-type) pseudobulk with fewer cells than this.")
    g2.add_argument("--min-samples-per-group", type=int, default=2,
                    help="A cell type needs >= this many samples in BOTH arms to be tested.")
    g2.add_argument("--padj-threshold", type=float, default=0.1,
                    help="FDR threshold for calling a factor differential.")
    g2.add_argument("--lfc-threshold", type=float, default=0.0,
                    help="Optional |log2FC| floor for calling a factor differential.")
    g2.add_argument("--celltypes-de", nargs="*", default=[],
                    help="Restrict DE to these cell types (default: all with enough samples).")

    g3 = p.add_argument_group("scope")
    g3.add_argument("--ligands-only", action="store_true",
                    help="Restrict the dual-role map to signature genes that are ligands "
                         "in the resource (the secreted-factor view). Recommended.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Pseudobulk differential expression per cell type
# ---------------------------------------------------------------------------
def _make_dds(counts_ct: pd.DataFrame, coldata_ct: pd.DataFrame, cond_col: str,
              cntrl: str, test: str, n_cpus: int = 4):
    """Construct a pydeseq2 DeseqDataSet across the supported signatures."""
    from pydeseq2.dds import DeseqDataSet
    counts_samples = counts_ct.T  # samples x genes, integer
    counts_samples = counts_samples.round().astype(int)
    md = coldata_ct.copy()
    md[cond_col] = pd.Categorical(md[cond_col], categories=[cntrl, test])
    # try the newer `design=` formula signature, then older variants
    try:
        return DeseqDataSet(counts=counts_samples, metadata=md,
                            design=f"~{cond_col}", n_cpus=n_cpus, quiet=True)
    except TypeError:
        pass
    try:
        return DeseqDataSet(counts=counts_samples, metadata=md,
                            design_factors=cond_col, ref_level=[cond_col, cntrl],
                            n_cpus=n_cpus, quiet=True)
    except TypeError:
        return DeseqDataSet(counts=counts_samples, clinical=md,
                            design_factors=cond_col, reference_level=cntrl,
                            n_cpus=n_cpus)


def _deseq_stats(dds, cond_col: str, cntrl: str, test: str):
    """Run Wald test for contrast TEST vs CNTRL and return results_df."""
    from pydeseq2.ds import DeseqStats
    dds.deseq2()
    try:
        from pydeseq2.default_inference import DefaultInference
        inference = DefaultInference(n_cpus=4)
        ds = DeseqStats(dds, contrast=[cond_col, test, cntrl], inference=inference, quiet=True)
    except Exception:
        ds = DeseqStats(dds, contrast=[cond_col, test, cntrl])
    ds.summary()
    res = ds.results_df.copy()
    res = res.rename(columns={"log2FoldChange": "log2FC", "pvalue": "pval", "padj": "padj"})
    return res


def pseudobulk_de_by_celltype(adata, cfg: cc.CCConfig, counts_layer: str, args) -> pd.DataFrame:
    """For each cell type, build a (sample x cell-type) pseudobulk and run DESeq2
    TEST vs CNTRL. Returns long-form DE table (gene, cell_type, log2FC, padj, ...)."""
    cc.banner("Pseudobulk DE per cell type (TEST vs CNTRL)")
    counts_df, coldata = cc.pseudobulk_matrix(
        adata, cfg, counts_layer,
        group_keys=[cfg.sample_key, cfg.celltype_key])
    # drop tiny pseudobulks
    coldata = coldata[coldata["_n_cells"] >= args.min_cells_pb]
    counts_df = counts_df[coldata.index]
    coldata[cfg.condition_key] = (
        adata.obs[[cfg.sample_key, cfg.condition_key]].astype(str)
        .drop_duplicates().set_index(cfg.sample_key)[cfg.condition_key]
        .reindex(coldata[cfg.sample_key].values).values)

    celltypes = (args.celltypes_de or
                 sorted(coldata[cfg.celltype_key].astype(str).unique()))
    all_res = []
    for ct in celltypes:
        sel = coldata[coldata[cfg.celltype_key].astype(str) == str(ct)]
        n_c = (sel[cfg.condition_key] == cfg.cntrl).sum()
        n_t = (sel[cfg.condition_key] == cfg.test).sum()
        if n_c < args.min_samples_per_group or n_t < args.min_samples_per_group:
            print(f"  skip {ct}: CNTRL={n_c} TEST={n_t} (need >= {args.min_samples_per_group} each)")
            continue
        c_ct = counts_df[sel.index]
        # filter genes with all-zero across this cell type's pseudobulks
        c_ct = c_ct.loc[c_ct.sum(axis=1) > 0]
        md = sel[[cfg.condition_key]].copy()
        print(f"  {ct}: CNTRL={n_c} TEST={n_t}, genes={c_ct.shape[0]}")
        try:
            dds = _make_dds(c_ct, md, cfg.condition_key, cfg.cntrl, cfg.test)
            res = _deseq_stats(dds, cfg.condition_key, cfg.cntrl, cfg.test)
        except Exception as e:
            warnings.warn(f"DESeq2 failed for {ct}: {e}")
            continue
        res["cell_type"] = str(ct)
        res["gene"] = res.index
        all_res.append(res.reset_index(drop=True))
    if not all_res:
        return pd.DataFrame()
    return pd.concat(all_res, ignore_index=True)


def export_pseudobulk_for_R(adata, cfg: cc.CCConfig, counts_layer: str, args, out_dir: str):
    """Write per-cell-type pseudobulk count matrices + colData for an R DESeq2/edgeR run."""
    cc.banner("Exporting pseudobulk matrices for R (DESeq2/edgeR)")
    counts_df, coldata = cc.pseudobulk_matrix(
        adata, cfg, counts_layer, group_keys=[cfg.sample_key, cfg.celltype_key])
    coldata = coldata[coldata["_n_cells"] >= args.min_cells_pb]
    counts_df = counts_df[coldata.index]
    coldata[cfg.condition_key] = (
        adata.obs[[cfg.sample_key, cfg.condition_key]].astype(str)
        .drop_duplicates().set_index(cfg.sample_key)[cfg.condition_key]
        .reindex(coldata[cfg.sample_key].values).values)
    pb_dir = cc.ensure_out(cfg, "pseudobulk")
    counts_df.to_csv(os.path.join(pb_dir, "pseudobulk_counts_genes_x_pb.csv"))
    coldata.to_csv(os.path.join(pb_dir, "pseudobulk_coldata.csv"))
    meta = {"sample_key": cfg.sample_key, "celltype_key": cfg.celltype_key,
            "condition_key": cfg.condition_key, "cntrl": cfg.cntrl, "test": cfg.test}
    with open(os.path.join(pb_dir, "pseudobulk_meta.json"), "w") as fh:
        import json
        json.dump(meta, fh, indent=2)
    print(f"  wrote pseudobulk matrices to {pb_dir}")
    print("  run:  Rscript run_pseudobulk_DE.R "
          f"{pb_dir}/pseudobulk_counts_genes_x_pb.csv {pb_dir}/pseudobulk_coldata.csv "
          f"{cfg.condition_key} {cfg.cntrl} {cfg.test} {pb_dir}/de_results.csv")


# ---------------------------------------------------------------------------
# Per-cell-type expression of each factor (who secretes it)
# ---------------------------------------------------------------------------
def expression_by_celltype(adata_ln, cfg: cc.CCConfig, genes) -> pd.DataFrame:
    """Mean log-normalised expression + fraction expressing, per cell type, for each gene."""
    import scipy.sparse as sp
    present = [g for g in genes if g in adata_ln.var_names]
    if not present:
        return pd.DataFrame()
    X = adata_ln[:, present].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    df = pd.DataFrame(X, columns=present, index=adata_ln.obs_names)
    df[cfg.celltype_key] = adata_ln.obs[cfg.celltype_key].astype(str).values
    mean_expr = df.groupby(cfg.celltype_key).mean()
    frac_expr = df.drop(columns=[cfg.celltype_key]).gt(0).groupby(
        adata_ln.obs[cfg.celltype_key].astype(str).values).mean()
    long = (mean_expr.reset_index().melt(id_vars=cfg.celltype_key,
            var_name="gene", value_name="mean_expr"))
    fr = (frac_expr.reset_index().rename(columns={"index": cfg.celltype_key})
          .melt(id_vars=cfg.celltype_key, var_name="gene", value_name="frac_expr"))
    return long.merge(fr, on=[cfg.celltype_key, "gene"], how="left")


def top_secretors(expr_long: pd.DataFrame, cfg: cc.CCConfig, n: int = 2) -> pd.DataFrame:
    """For each gene, the top-n cell types by mean expression."""
    if expr_long.empty:
        return pd.DataFrame()
    out = []
    for gene, sub in expr_long.groupby("gene"):
        top = sub.sort_values("mean_expr", ascending=False).head(n)
        out.append({"gene": gene,
                    "top_secretors": "; ".join(
                        f"{r[cfg.celltype_key]}({r['mean_expr']:.2f})"
                        for _, r in top.iterrows())})
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Assemble the dual-role table
# ---------------------------------------------------------------------------
def receptor_targets_for_ligand(focused: pd.DataFrame, ligand: str,
                                fallback: pd.DataFrame = None) -> str:
    """From the focused L-R table, summarise receptor(s) and receiver(s) for a ligand.
    If the ligand is absent from the focused view, fall back to the full spatial
    LIANA frame (so factors whose secretor sits outside the focused sender set still
    get a receptor/receiver annotation)."""
    def _summarise(tab):
        if tab is None or tab.empty or "ligand_complex" not in tab.columns:
            return ""
        sub = tab[tab["ligand_complex"].astype(str).str.split("_").str[0].eq(ligand) |
                  tab["ligand_complex"].astype(str).eq(ligand)]
        if sub.empty:
            return ""
        pairs = (sub.assign(pair=lambda d: d["receptor_complex"].astype(str) +
                            " on " + d["target"].astype(str))
                 ["pair"].value_counts().head(4).index.tolist())
        return "; ".join(pairs)

    hit = _summarise(focused)
    if hit:
        return hit
    return _summarise(fallback)


def build_dual_role_table(de: pd.DataFrame, expr_long: pd.DataFrame, focused: pd.DataFrame,
                          roles: pd.DataFrame, genes, ligand_set: set,
                          cfg: cc.CCConfig, args, full_liana: pd.DataFrame = None) -> pd.DataFrame:
    """One row per signature factor: direction, secretor, receptor/receiver, role."""
    cc.banner("Assembling the dual-role factor table")
    # universe of factors = signature genes (optionally restricted to ligands)
    factors = [g for g in genes if (not args.ligands_only or g in ligand_set)]
    if not factors:
        warnings.warn("No signature factors after filtering; check --genes / --ligands-only.")
        return pd.DataFrame()

    sec = top_secretors(expr_long, cfg)
    role_map = {}
    if roles is not None and not roles.empty:
        rc = {c.lower(): c for c in roles.columns}
        fcol = rc.get("factor", "factor"); rcol = rc.get("role", "role")
        role_map = dict(zip(roles[fcol].astype(str), roles[rcol].astype(str).str.lower()))

    rows = []
    for g in factors:
        sub = de[de["gene"].astype(str) == g] if not de.empty else pd.DataFrame()
        # pick the cell type with the most significant / largest-magnitude change
        if not sub.empty:
            sub = sub.assign(_score=lambda d: -np.log10(d["padj"].fillna(1) + 1e-300) *
                             d["log2FC"].abs())
            best = sub.sort_values("_score", ascending=False).iloc[0]
            l2fc = float(best["log2FC"]); padj = float(best["padj"]) if pd.notna(best["padj"]) else np.nan
            de_ct = str(best["cell_type"])
            n_sig_ct = int((sub["padj"] < args.padj_threshold).sum())
            # overall direction by sign-consensus across cell types where significant
            sig = sub[sub["padj"] < args.padj_threshold]
            if not sig.empty:
                direction = cfg.test if sig["log2FC"].mean() > 0 else cfg.cntrl
            else:
                direction = cfg.test if l2fc > 0 else cfg.cntrl
        else:
            l2fc, padj, de_ct, n_sig_ct, direction = np.nan, np.nan, "", 0, "NA"

        differential = (pd.notna(padj) and padj < args.padj_threshold and
                        abs(l2fc) >= args.lfc_threshold)
        secretor = (sec.loc[sec["gene"] == g, "top_secretors"].iloc[0]
                    if (not sec.empty and (sec["gene"] == g).any()) else "")
        receptors = receptor_targets_for_ligand(focused, g, fallback=full_liana)
        role = role_map.get(g, "unassigned")
        rows.append({
            "factor": g,
            "is_ligand": g in ligand_set,
            "log2FC_TEST_vs_CNTRL": l2fc,
            "padj": padj,
            "differential": bool(differential),
            "enriched_in": (direction if differential else "ns"),
            "DE_cell_type": de_ct,
            "n_significant_celltypes": n_sig_ct,
            "top_secretors": secretor,
            "receptor_on_receiver": receptors,
            "literature_role": role,
        })
    out = pd.DataFrame(rows)

    # the headline dual-role cross-tab: role x enrichment direction
    if not out.empty:
        out["dual_role_call"] = out.apply(
            lambda r: _dual_role_label(r, cfg), axis=1)
    return out.sort_values(["differential", "log2FC_TEST_vs_CNTRL"],
                           ascending=[False, False])


def _dual_role_label(r, cfg: cc.CCConfig) -> str:
    """Translate (role, direction) into an interpretable call without hard-coding biology."""
    role = str(r["literature_role"])
    enr = str(r["enriched_in"])
    if enr == "ns":
        return "not differential"
    if role == "pro":
        return f"pro-process & enriched in {enr}"
    if role == "anti":
        return f"anti-process & enriched in {enr}"
    return f"role-unassigned & enriched in {enr}"


# ---------------------------------------------------------------------------
# Figures: the Fig-7D/E analog
# ---------------------------------------------------------------------------
def plot_dual_role_map(table: pd.DataFrame, cfg: cc.CCConfig, out_dir: str):
    if table.empty:
        return
    plt = cc.set_plot_theme()
    role_colors = {"pro": "#2E7D32", "anti": "#C62828", "context": "#6A1B9A",
                   "unassigned": "#9E9E9E"}
    d = table.dropna(subset=["log2FC_TEST_vs_CNTRL"]).copy()
    if d.empty:
        return
    d = d.sort_values("log2FC_TEST_vs_CNTRL")
    y = np.arange(len(d))
    colors = [role_colors.get(str(r), "#9E9E9E") for r in d["literature_role"]]
    fig, ax = plt.subplots(figsize=(6, 2 + 0.22 * len(d)))
    ax.barh(y, d["log2FC_TEST_vs_CNTRL"].values, color=colors,
            edgecolor="k", lw=0.3)
    # mark significance
    for yi, (_, r) in zip(y, d.iterrows()):
        if r["differential"]:
            ax.text(r["log2FC_TEST_vs_CNTRL"] +
                    (0.05 if r["log2FC_TEST_vs_CNTRL"] >= 0 else -0.05),
                    yi, cc.pval_stars(r["padj"]), va="center",
                    ha="left" if r["log2FC_TEST_vs_CNTRL"] >= 0 else "right",
                    fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(d["factor"], fontsize=6)
    ax.set_xlabel(f"log2FC  (TEST={cfg.test}  vs  CNTRL={cfg.cntrl})")
    ax.set_title("Signalling-factor dual-role map\n"
                 f"left = enriched in {cfg.cntrl}   |   right = enriched in {cfg.test}")
    # legend
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=v, label=k) for k, v in role_colors.items()],
              title="literature role", fontsize=6, title_fontsize=7,
              loc="lower right", frameon=False)
    fig.savefig(os.path.join(out_dir, "dual_role_factor_map.pdf"))
    plt.close(fig)


def plot_role_direction_crosstab(table: pd.DataFrame, cfg: cc.CCConfig, out_dir: str):
    """The summary cross-tab: do anti-process factors concentrate in TEST and pro- in CNTRL?"""
    if table.empty:
        return
    plt = cc.set_plot_theme()
    d = table[table["differential"]].copy()
    if d.empty:
        return
    ct = pd.crosstab(d["literature_role"], d["enriched_in"])
    fig, ax = plt.subplots(figsize=(2.2 + 0.5 * ct.shape[1], 2 + 0.4 * ct.shape[0]))
    im = ax.imshow(ct.values, cmap="Blues", aspect="auto")
    ax.set_xticks(range(ct.shape[1])); ax.set_xticklabels(ct.columns)
    ax.set_yticks(range(ct.shape[0])); ax.set_yticklabels(ct.index)
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            ax.text(j, i, int(ct.values[i, j]), ha="center", va="center", fontsize=9)
    ax.set_xlabel("enriched in"); ax.set_ylabel("literature role")
    ax.set_title("Differential factors:\nrole x direction")
    fig.colorbar(im, ax=ax, shrink=0.6, label="n factors")
    fig.savefig(os.path.join(out_dir, "role_by_direction_crosstab.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    cfg = cc.config_from_args(args)
    cfg.resource_name = args.resource
    out_dir = cc.ensure_out(cfg)

    cc.banner("LOADING DATA")
    adata = cc.load_adata(args.h5ad)
    cfg = cc.resolve_keys(adata, cfg, require_niche=False)
    counts_layer = cc.resolve_counts_layer(adata, cfg)
    genes = cc.load_condition_genes(args.genes_path)
    if not genes:
        raise SystemExit("This analysis needs a signature: pass --genes condition_genes.txt")
    print(f"  cells={adata.n_obs} genes={adata.n_vars}; signature factors={len(genes)}")
    print(f"  CNTRL={cfg.cntrl} TEST={cfg.test}; counts layer={counts_layer}")

    # which signature genes are secreted ligands?
    try:
        ligand_set, _ = cc.ligand_receptor_universe(args.resource)
    except Exception as e:
        warnings.warn(f"Could not load LIANA resource ({e}); treating all factors as ligands.")
        ligand_set = set(genes)
    n_lig = len([g for g in genes if g in ligand_set])
    print(f"  signature factors that are ligands in '{args.resource}': {n_lig}/{len(genes)}")

    # log-normalised view for expression-by-cell-type
    adata_ln = cc.make_lognorm_view(adata, cfg, counts_layer)

    # (i) DE by sample (pseudobulk) ---------------------------------------
    if args.de_backend == "export":
        export_pseudobulk_for_R(adata, cfg, counts_layer, args, out_dir)
        de = pd.DataFrame()
    else:
        try:
            de = pseudobulk_de_by_celltype(adata, cfg, counts_layer, args)
        except Exception as e:
            warnings.warn(f"pydeseq2 unavailable/failed ({e}); exporting for R instead.")
            export_pseudobulk_for_R(adata, cfg, counts_layer, args, out_dir)
            de = pd.DataFrame()
    if not de.empty:
        de.to_csv(os.path.join(out_dir, "pseudobulk_DE_all_genes_by_celltype.csv"), index=False)
        de[de["gene"].isin(genes)].to_csv(
            os.path.join(out_dir, "pseudobulk_DE_signature_factors.csv"), index=False)

    # (ii) who secretes ----------------------------------------------------
    expr_long = expression_by_celltype(adata_ln, cfg, genes)
    if not expr_long.empty:
        expr_long.to_csv(os.path.join(out_dir, "factor_expression_by_celltype.csv"), index=False)

    # (iii) receptor/receiver from the communication step ------------------
    focused = None
    if args.focused_csv and os.path.exists(args.focused_csv):
        focused = pd.read_csv(args.focused_csv)
        print(f"  loaded focused L-R interactions: {len(focused)} rows")
    else:
        warnings.warn("No --focused-csv given; receptor/receiver column will be empty. "
                      "Run analysis_cellcell_communication.py first for the full picture.")
    # full spatial LIANA frame -> fallback receptor/receiver source for factors
    # whose secretor falls outside the focused sender set.
    liana_full = None
    if args.liana_csv and os.path.exists(args.liana_csv):
        liana_full = pd.read_csv(args.liana_csv)
        print(f"  loaded full spatial L-R table (receptor-join fallback): {len(liana_full)} rows")

    # (iv) literature roles ------------------------------------------------
    roles = None
    if args.roles_csv and os.path.exists(args.roles_csv):
        # the template carries a leading '#'-comment block documenting the schema;
        # comment='#' skips it (and any inline '#' notes) instead of mis-parsing it.
        roles = pd.read_csv(args.roles_csv, comment="#")
        roles = roles.dropna(how="all")
        print(f"  loaded role annotations for {len(roles)} factors")
    else:
        warnings.warn("No --roles-csv; literature_role will be 'unassigned'. "
                      "Edit factor_roles_template.csv and pass it in.")

    # assemble + plot ------------------------------------------------------
    table = build_dual_role_table(de, expr_long, focused, roles, genes,
                                  ligand_set, cfg, args, full_liana=liana_full)
    if not table.empty:
        table.to_csv(os.path.join(out_dir, "factor_dual_role_table.csv"), index=False)
        print("\nDual-role table (head):")
        print(table.head(25).to_string(index=False))

    fig_dir = cc.ensure_out(cfg, "figures")
    try:
        plot_dual_role_map(table, cfg, fig_dir)
        plot_role_direction_crosstab(table, cfg, fig_dir)
    except Exception as e:
        warnings.warn(f"Figure generation issue: {e}")

    cc.write_report(cfg, {
        "n_signature_factors": len(genes),
        "n_ligand_factors": int(n_lig),
        "de_backend": args.de_backend,
        "n_de_rows": int(len(de)),
        "n_factors_tabulated": int(len(table)),
        "n_differential": int(table["differential"].sum()) if not table.empty else 0,
        "padj_threshold": args.padj_threshold,
        "caveats": cc.CAVEATS,
    }, "condition_factor_roles_report.json")

    print(cc.CAVEATS)
    cc.banner("DONE: signalling-factor dual-role analysis")
    print(f"Outputs in: {out_dir}")


if __name__ == "__main__":
    main()
