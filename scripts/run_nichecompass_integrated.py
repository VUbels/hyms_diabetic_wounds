"""
run_nichecompass_integrated.py

NicheCompass niche identification pipeline for spatial transcriptomics data.

Loads converted h5ad files from annotated_data/, integrates all 6 samples,
trains NicheCompass, identifies niches via Leiden clustering, and saves results.

Usage:
    python run_nichecompass_integrated.py [--base_dir ./annotated_data] [--output_dir ./nichecompass_results]

Requires:
    pip install nichecompass scanpy squidpy anndata
    (with PyTorch + PyG already installed for GPU support)
"""

import os
import sys
import glob
import argparse
import warnings
from datetime import datetime

# Force non-interactive backend BEFORE importing pyplot anywhere
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import scipy.sparse as sp

# Suppress scanpy/squidpy interactive plot attempts
sc.settings.autoshow = False

from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    create_new_color_dict,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
    filter_and_combine_gp_dict_gps_v2,
)

warnings.filterwarnings("ignore")


# =============================================================================
# CONFIG — edit these if needed
# =============================================================================

# Data keys (matching your h5ad structure)
SPATIAL_KEY = "X_spatial"          # obsm key for proseg centroid coordinates
CELL_TYPE_KEY = "predicted_cell_type"  # obs column for cell annotations
COUNTS_KEY = "counts"              # layer key NicheCompass expects for raw counts

# Condition mapping — token -> condition label.
# Matched against sample_id using the mode set by CONDITION_MAP_MODE:
#   "prefix" — match the part before the first underscore (e.g. "Cntrl" in "Cntrl_Region_2")
#   "suffix" — match the part after the last hyphen or underscore (e.g. "L" in "XE765-L")
CONDITION_MAP = {
    "L": "lesion",
    "D": "diabetic",
}
CONDITION_MAP_MODE = "suffix"   # "prefix" or "suffix"

# These two values must exactly match what CONDITION_MAP produces.
# Used wherever condition labels are referenced in analysis and plots.
_cond_labels      = list(CONDITION_MAP.values())
CONDITION_CNTRL   = _cond_labels[0]
CONDITION_TEST    = _cond_labels[1] 

# Spatial graph
N_NEIGHBORS = 8  # 8 is recommended for single-cell resolution spatial data

# Gene programs
SPECIES = "mouse"

# NicheCompass AnnData keys
ADJ_KEY = "spatial_connectivities"
GP_NAMES_KEY = "nichecompass_gp_names"
ACTIVE_GP_NAMES_KEY = "nichecompass_active_gp_names"
GP_TARGETS_MASK_KEY = "nichecompass_gp_targets"
GP_TARGETS_CATEGORIES_MASK_KEY = "nichecompass_gp_targets_categories"
GP_SOURCES_MASK_KEY = "nichecompass_gp_sources"
GP_SOURCES_CATEGORIES_MASK_KEY = "nichecompass_gp_sources_categories"
LATENT_KEY = "nichecompass_latent"

# Architecture — GATv2 recommended for single-cell resolution data
CONV_LAYER_ENCODER = "gcnconv"
ACTIVE_GP_THRESH_RATIO = 0.01

# Training
N_EPOCHS = 300
N_EPOCHS_ALL_GPS = 25
LR = 0.0001
LAMBDA_EDGE_RECON = 5000.0         # Reduced from 50000 to prevent spatial-graph
                                    # dominance over gene expression reconstruction
LAMBDA_GENE_EXPR_RECON = 300.0
LAMBDA_L1_MASKED = 50.0            # Enforce sparsity on masked GP weights for
                                    # interpretable latent dimensions
LAMBDA_L1_ADDON = 30.0
EDGE_BATCH_SIZE = 512             
N_SAMPLED_NEIGHBORS = 8

# Clustering
LATENT_LEIDEN_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.2, 1.5]
LATENT_LEIDEN_RESOLUTION = 0.8     # Fallback default
LATENT_CLUSTER_KEY = "nichecompass_niche"


# =============================================================================
# STEP 1: Load and prepare data
# =============================================================================

def discover_h5ad_files(base_dir):
    """Find all _annotated.h5ad files and assign sample IDs and conditions."""
    pattern = os.path.join(base_dir, "*", "*_annotated.h5ad")
    files = sorted(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No *_annotated.h5ad files found under {base_dir}")

    samples = []
    for f in files:
        # Extract sample name from parent directory (e.g., "XE765-L" or "Cntrl_Region_2")
        sample_id = os.path.basename(os.path.dirname(f))

        if CONDITION_MAP_MODE == "suffix":
            # Split on the last hyphen or underscore to get the trailing token
            import re
            parts = re.split(r"[-_]", sample_id)
            token = parts[-1] if parts else sample_id
        else:  # "prefix"
            token = sample_id.split("_")[0]

        condition = CONDITION_MAP.get(token)
        if condition is None:
            print(f"[WARN] Unknown {CONDITION_MAP_MODE} token '{token}' "
                  f"for sample '{sample_id}', skipping")
            continue
        samples.append({
            "sample_id": sample_id,
            "condition": condition,
            "path": f,
        })

    print(f"Discovered {len(samples)} samples:")
    for s in samples:
        print(f"  {s['sample_id']} ({s['condition']}): {s['path']}")

    return samples


def load_and_prepare_sample(sample_info):
    """
    Load a single h5ad, add metadata, prepare counts layer and spatial graph.
    """
    sample_id = sample_info["sample_id"]
    condition = sample_info["condition"]
    path = sample_info["path"]

    print(f"\n[LOAD] {sample_id} from {path}")
    adata = sc.read_h5ad(path)

    # --- Add sample/condition metadata ---
    adata.obs["sample_id"] = sample_id
    adata.obs["condition"] = condition

    # --- Prepare counts layer ---
    # scCustomize put raw counts into .X; NicheCompass expects a named layer.
    # CRITICAL: cast to float32 — NicheCompass / PyTorch expects float32,
    # but h5ad files often store .X as float64 which causes dtype mismatch
    # errors in the encoder ("mat1 and mat2 must have the same dtype").
    if COUNTS_KEY not in adata.layers:
        if sp.issparse(adata.X):
            adata.layers[COUNTS_KEY] = adata.X.astype(np.float32).copy()
        else:
            adata.layers[COUNTS_KEY] = adata.X.astype(np.float32).copy()
        print(f"  Copied .X to .layers['{COUNTS_KEY}'] (cast to float32)")
    else:
        # Ensure existing counts layer is also float32
        if sp.issparse(adata.layers[COUNTS_KEY]):
            if adata.layers[COUNTS_KEY].dtype != np.float32:
                adata.layers[COUNTS_KEY] = adata.layers[COUNTS_KEY].astype(np.float32)
                print(f"  Cast existing .layers['{COUNTS_KEY}'] to float32")
        else:
            if adata.layers[COUNTS_KEY].dtype != np.float32:
                adata.layers[COUNTS_KEY] = adata.layers[COUNTS_KEY].astype(np.float32)
                print(f"  Cast existing .layers['{COUNTS_KEY}'] to float32")

    # Also ensure .X itself is float32 for consistency
    if sp.issparse(adata.X):
        if adata.X.dtype != np.float32:
            adata.X = adata.X.astype(np.float32)
    else:
        if adata.X.dtype != np.float32:
            adata.X = adata.X.astype(np.float32)

    # --- Ensure spatial coordinates exist ---
    if SPATIAL_KEY not in adata.obsm:
        raise KeyError(f"Spatial key '{SPATIAL_KEY}' not in obsm for {sample_id}")

    # --- Build per-sample spatial graph ---
    # Must be done BEFORE concatenation so samples remain disconnected
    sq.gr.spatial_neighbors(
        adata,
        spatial_key=SPATIAL_KEY,
        n_neighs=N_NEIGHBORS,
        coord_type="generic",
    )
    n_edges = adata.obsp["spatial_connectivities"].nnz
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes, {n_edges} spatial edges")

    # --- Ensure var_names are unique ---
    adata.var_names_make_unique()

    # --- Make obs index unique by prefixing sample_id ---
    adata.obs_names = [f"{sample_id}_{i}" for i in adata.obs_names]

    return adata


def concatenate_samples(adatas):
    """
    Concatenate per-sample AnnDatas with block-diagonal spatial graph.
    """
    from scipy.sparse import block_diag

    print(f"\n[CONCAT] Merging {len(adatas)} samples...")

    # Collect spatial graphs before concat
    connectivities = [a.obsp["spatial_connectivities"] for a in adatas]
    distances = [a.obsp["spatial_distances"] for a in adatas]

    # Concatenate — use inner join on var to handle slight gene differences
    adata = ad.concat(adatas, join="inner", merge="same")

    # Rebuild block-diagonal spatial graph (samples stay disconnected)
    adata.obsp["spatial_connectivities"] = block_diag(connectivities, format="csr")
    adata.obsp["spatial_distances"] = block_diag(distances, format="csr")

    # Ensure categorical types
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    adata.obs["condition"] = adata.obs["condition"].astype("category")
    adata.obs[CELL_TYPE_KEY] = adata.obs[CELL_TYPE_KEY].astype("category")

    print(f"  Combined: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"  Samples: {adata.obs['sample_id'].value_counts().to_dict()}")
    print(f"  Conditions: {adata.obs['condition'].value_counts().to_dict()}")

    return adata


# =============================================================================
# STEP 2: Build gene program prior mask
# =============================================================================

def build_gene_programs(adata):
    """
    Extract and combine gene program (GP) databases for the prior mask.
    NicheCompass uses these to make its latent space interpretable.

    MEBOCOST requires bundled data files that may not be present outside
    the NicheCompass repo. If it fails, we proceed with OmniPath + NicheNet
    which already provide strong ligand-receptor coverage.
    """
    print("\n[GP] Building gene program prior mask...")

    # Extract from three databases
    print("  Extracting OmniPath LR interactions...")
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=SPECIES,
        min_curation_effort=0,
    )

    print("  Extracting NicheNet LRT interactions...")
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=SPECIES,
        version="v2",
        keep_target_genes_ratio=1.,
        max_n_target_genes_per_gp=250,
    )

    mebocost_gp_dict = None
    try:
        print("  Extracting MEBOCOST MS interactions...")
        mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
            species=SPECIES,
        )
    except (FileNotFoundError, Exception) as e:
        print(f"  [WARN] MEBOCOST extraction failed: {e}")
        print("  Proceeding with OmniPath + NicheNet only.")

    # Combine and filter — function takes a positional list of GP dicts,
    # NOT named keyword arguments
    print("  Combining and filtering gene programs...")
    gp_dicts = [omnipath_gp_dict, nichenet_gp_dict]
    if mebocost_gp_dict is not None:
        gp_dicts.append(mebocost_gp_dict)

    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
        gp_dicts,
        verbose=True,
    )

    # Add GPs to adata — note: active_gp_names_key is NOT a parameter of
    # this function (it is only used at NicheCompass model init)
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_gp_dict,
        adata=adata,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        gp_names_key=GP_NAMES_KEY,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None,
    )

    n_gps = len(adata.uns.get(GP_NAMES_KEY, []))
    print(f"  {n_gps} total gene programs")

    return adata


# =============================================================================
# STEP 3: Train NicheCompass
# =============================================================================

def train_nichecompass(adata):
    """Initialize and train the NicheCompass model."""
    print("\n[TRAIN] Initializing NicheCompass model...")

    model = NicheCompass(
        adata,
        counts_key=COUNTS_KEY,
        adj_key=ADJ_KEY,
        cat_covariates_embeds_injection=["gene_expr_decoder"],
        cat_covariates_keys=["sample_id"],
        cat_covariates_no_edges=[True],
        cat_covariates_embeds_nums=[len(adata.obs["sample_id"].cat.categories)],
        gp_names_key=GP_NAMES_KEY,
        active_gp_names_key=ACTIVE_GP_NAMES_KEY,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        latent_key=LATENT_KEY,
        conv_layer_encoder=CONV_LAYER_ENCODER,
        active_gp_thresh_ratio=ACTIVE_GP_THRESH_RATIO,
    )

    print("[TRAIN] Starting training...")
    model.train(
        n_epochs=N_EPOCHS,
        n_epochs_all_gps=N_EPOCHS_ALL_GPS,
        lr=LR,
        lambda_edge_recon=LAMBDA_EDGE_RECON,
        lambda_gene_expr_recon=LAMBDA_GENE_EXPR_RECON,
        lambda_l1_masked=LAMBDA_L1_MASKED,
        lambda_l1_addon=LAMBDA_L1_ADDON,
        edge_batch_size=EDGE_BATCH_SIZE,
        n_sampled_neighbors=N_SAMPLED_NEIGHBORS,
        use_cuda_if_available=True,
        verbose=True,
    )

    print("[TRAIN] Training complete.")
    return model


# =============================================================================
# STEP 4: Cluster niches and compute embeddings
# =============================================================================

def cluster_niches(model, resolution=LATENT_LEIDEN_RESOLUTION):
    """
    Compute neighbor graph on latent space, sweep Leiden resolutions with
    silhouette scoring, run UMAP. Uses the resolution with the best
    silhouette score from LATENT_LEIDEN_RESOLUTIONS; falls back to the
    provided resolution if the sweep is empty.
    """
    from sklearn.metrics import silhouette_score

    print(f"\n[CLUSTER] Leiden clustering with resolution sweep...")

    adata = model.adata

    # Neighbor graph on NicheCompass latent space
    sc.pp.neighbors(
        adata,
        use_rep=LATENT_KEY,
        key_added=LATENT_KEY,
    )

    # Resolution sweep with silhouette scoring
    resolutions = LATENT_LEIDEN_RESOLUTIONS
    latent_rep = adata.obsm[LATENT_KEY]

    best_res = resolution
    best_sil = -1
    sweep_results = []

    for res in resolutions:
        temp_key = f"_leiden_sweep_{res}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=temp_key,
            neighbors_key=LATENT_KEY,
        )
        labels = adata.obs[temp_key].values
        n_clusters = len(np.unique(labels))

        if n_clusters < 2 or n_clusters >= adata.n_obs - 1:
            sil = -1
        else:
            # Subsample for speed if large dataset
            if adata.n_obs > 50000:
                rng = np.random.default_rng(42)
                idx = rng.choice(adata.n_obs, size=50000, replace=False)
                sil = silhouette_score(latent_rep[idx], labels[idx],
                                       metric="euclidean", sample_size=None)
            else:
                sil = silhouette_score(latent_rep, labels, metric="euclidean")

        sweep_results.append({
            "resolution": res,
            "n_clusters": n_clusters,
            "silhouette": sil,
        })
        print(f"  res={res:.2f}: {n_clusters} clusters, silhouette={sil:.4f}")

        if sil > best_sil:
            best_sil = sil
            best_res = res

        # Clean up temp key
        del adata.obs[temp_key]

    print(f"\n  Best resolution: {best_res} (silhouette={best_sil:.4f})")

    # Final clustering with best resolution
    sc.tl.leiden(
        adata,
        resolution=best_res,
        key_added=LATENT_CLUSTER_KEY,
        neighbors_key=LATENT_KEY,
    )

    # UMAP for visualization
    sc.tl.umap(adata, neighbors_key=LATENT_KEY)

    n_niches = adata.obs[LATENT_CLUSTER_KEY].nunique()
    print(f"  Identified {n_niches} niches at resolution {best_res}")

    # Store sweep results for inspection
    adata.uns["leiden_resolution_sweep"] = pd.DataFrame(sweep_results)
    adata.uns["leiden_best_resolution"] = best_res

    # Summarize niche composition
    ct_by_niche = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs[CELL_TYPE_KEY],
        normalize="index",
    )
    print("\n  Top cell types per niche:")
    for niche in ct_by_niche.index:
        top3 = ct_by_niche.loc[niche].nlargest(3)
        top_str = ", ".join([f"{ct} ({pct:.0%})" for ct, pct in top3.items()])
        print(f"    Niche {niche}: {top_str}")

    return adata


# =============================================================================
# STEP 4b: Gene program differential analysis
# =============================================================================

def analyze_gene_programs(adata, model, output_dir):
    """
    Extract and analyze NicheCompass gene program (GP) activity scores.

    NicheCompass encodes each cell's niche state as a vector of GP activities,
    where each GP corresponds to a known ligand-receptor or metabolite-sensor
    interaction. By comparing GP activities between conditions, we directly
    identify which communication axes change between CNTRL and TEST conditions.

    Produces:
        - Per-cell GP activity matrix (saved in adata.obsm)
        - Per-niche, per-condition mean GP activity
        - Differential GP activity between conditions (pseudobulk by sample)
        - Top differential GPs per niche
    """
    from scipy.stats import mannwhitneyu, ranksums

    print("\n[GP] Analyzing gene program activities...")

    gp_dir = os.path.join(output_dir, "gene_programs")
    os.makedirs(gp_dir, exist_ok=True)

    # Extract GP activity scores from the latent space
    # NicheCompass latent dimensions correspond to gene programs
    gp_names = adata.uns.get(GP_NAMES_KEY, [])
    active_gp_names = adata.uns.get(ACTIVE_GP_NAMES_KEY, [])

    latent = adata.obsm[LATENT_KEY]

    if len(active_gp_names) == 0:
        print("  [WARN] No active GP names found; using all GP names")
        active_gp_names = gp_names

    # The latent space has 2 * n_active_gps dimensions (source + target per GP)
    n_active = len(active_gp_names)
    n_latent = latent.shape[1]

    print(f"  {n_active} active GPs, {n_latent} latent dimensions")

    # If latent dims = 2 * n_active_gps, first half = source, second = target
    # Aggregate source + target as total GP activity per cell
    if n_latent == 2 * n_active:
        gp_source = latent[:, :n_active]
        gp_target = latent[:, n_active:]
        gp_activity = np.abs(gp_source) + np.abs(gp_target)
        gp_columns = list(active_gp_names)
        print("  Aggregated source + target GP activities")
    elif n_latent == n_active:
        gp_activity = np.abs(latent)
        gp_columns = list(active_gp_names)
    else:
        # Fallback: use raw latent dimensions
        gp_activity = np.abs(latent)
        gp_columns = [f"latent_{i}" for i in range(n_latent)]
        print(f"  [WARN] Latent dim ({n_latent}) != 2 * active GPs ({n_active}); "
              f"using raw latent dimensions")

    gp_df = pd.DataFrame(gp_activity, index=adata.obs_names, columns=gp_columns)

    # ---- Per-niche, per-condition mean GP activity ----
    meta = adata.obs[[LATENT_CLUSTER_KEY, "condition", "sample_id"]].copy()
    meta = meta.join(gp_df)

    niche_cond_gp = (
        meta.groupby([LATENT_CLUSTER_KEY, "condition"])[gp_columns]
        .mean()
    )
    niche_cond_gp.to_csv(os.path.join(gp_dir, "gp_activity_niche_condition.csv"))
    print(f"  Saved niche × condition GP activity means")

    # ---- Differential GP activity: sample-level pseudobulk ----
    # Aggregate GP activity per sample per niche, then compare conditions
    sample_niche_gp = (
        meta.groupby([LATENT_CLUSTER_KEY, "sample_id", "condition"])[gp_columns]
        .mean()
        .reset_index()
    )
    sample_niche_gp.to_csv(
        os.path.join(gp_dir, "gp_activity_sample_niche.csv"), index=False
    )

    niches = sorted(adata.obs[LATENT_CLUSTER_KEY].unique(),
                    key=lambda x: int(x))

    diff_gp_rows = []
    for niche in niches:
        niche_data = sample_niche_gp[
            sample_niche_gp[LATENT_CLUSTER_KEY] == niche
        ]
        hs_data = niche_data[niche_data["condition"] == CONDITION_CNTRL]
        ls_data = niche_data[niche_data["condition"] == CONDITION_TEST]

        if len(hs_data) < 2 or len(ls_data) < 2:
            continue

        for gp in gp_columns:
            hs_vals = hs_data[gp].values
            ls_vals = ls_data[gp].values

            mean_ctrl = hs_vals.mean()
            mean_test = ls_vals.mean()
            log2fc = np.log2((mean_ctrl + 1e-6) / (mean_test + 1e-6))

            # Wilcoxon rank-sum (n=3 vs 3; limited power, exploratory)
            try:
                _, p_val = ranksums(hs_vals, ls_vals)
            except Exception:
                p_val = np.nan

            diff_gp_rows.append({
                "niche": niche,
                "gene_program": gp,
                f"mean_{CONDITION_CNTRL}": mean_ctrl,
                f"mean_{CONDITION_TEST}": mean_test,
                f"log2FC_{CONDITION_CNTRL}_vs_{CONDITION_TEST}": log2fc,
                "ranksums_p": p_val,
                f"n_{CONDITION_CNTRL}_samples": len(hs_vals),
                f"n_{CONDITION_TEST}_samples": len(ls_vals),
            })

    diff_gp_df = pd.DataFrame(diff_gp_rows)

    if len(diff_gp_df) > 0:
        # FDR correction per niche (separate multiple testing domains)
        from statsmodels.stats.multitest import multipletests
        diff_gp_df["padj"] = np.nan
        for niche in diff_gp_df["niche"].unique():
            mask = (diff_gp_df["niche"] == niche) & diff_gp_df["ranksums_p"].notna()
            if mask.sum() > 0:
                diff_gp_df.loc[mask, "padj"] = multipletests(
                    diff_gp_df.loc[mask, "ranksums_p"], method="fdr_bh"
                )[1]

        lfc_col = f"log2FC_{CONDITION_CNTRL}_vs_{CONDITION_TEST}"
        diff_gp_df = diff_gp_df.sort_values(
            ["niche", lfc_col],
            key=lambda x: x.abs() if x.name == lfc_col else x,
            ascending=[True, False],
        )
        diff_gp_df.to_csv(
            os.path.join(gp_dir, f"gp_differential_{CONDITION_CNTRL}_vs_{CONDITION_TEST}.csv"),
            index=False,
        )

        # Top GPs per niche (by absolute log2FC)
        top_n = 10
        print(f"\n  Top {top_n} differential GPs per niche (by |log2FC|):")
        for niche in niches:
            niche_gps = diff_gp_df[diff_gp_df["niche"] == niche].head(top_n)
            if len(niche_gps) == 0:
                continue
            top_gp = niche_gps.iloc[0]
            print(f"    Niche {niche}: {top_gp['gene_program']} "
                  f"(log2FC={top_gp[lfc_col]:.2f}, "
                  f"p={top_gp['ranksums_p']:.3f})")

    # ---- GP activity heatmap: top variable GPs across niches ----
    import matplotlib.pyplot as plt
    import seaborn as sns

    niche_gp_mean = meta.groupby(LATENT_CLUSTER_KEY)[gp_columns].mean()
    gp_var = niche_gp_mean.var(axis=0)
    top_variable_gps = gp_var.nlargest(min(50, len(gp_var))).index.tolist()

    if len(top_variable_gps) > 0:
        plot_data = niche_gp_mean[top_variable_gps]
        fig, ax = plt.subplots(
            figsize=(max(12, len(top_variable_gps) * 0.25),
                     max(4, len(niches) * 0.35))
        )
        sns.heatmap(
            plot_data, cmap="viridis", ax=ax, xticklabels=True,
            yticklabels=True, linewidths=0.2,
            cbar_kws={"label": "Mean GP activity"},
        )
        ax.set_title("Top variable gene programs across niches", fontsize=11)
        ax.set_ylabel("Niche")
        ax.set_xlabel("Gene Program")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
        plt.setp(ax.get_yticklabels(), fontsize=8)
        plt.tight_layout()
        fig.savefig(os.path.join(gp_dir, "gp_activity_heatmap.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close(fig)

    # ---- Differential GP heatmap (log2FC across niches × top GPs) ----
    if len(diff_gp_df) > 0:
        # Select top differential GPs across all niches
        top_diff_gps = (
            diff_gp_df.groupby("gene_program")[lfc_col]
            .apply(lambda x: x.abs().max())
            .nlargest(min(40, len(diff_gp_df["gene_program"].unique())))
            .index.tolist()
        )
        lfc_pivot = diff_gp_df.pivot(
            index="niche", columns="gene_program", values=lfc_col
        )
        lfc_plot = lfc_pivot.reindex(columns=top_diff_gps).reindex(niches)

        fig, ax = plt.subplots(
            figsize=(max(10, len(top_diff_gps) * 0.25),
                     max(4, len(niches) * 0.35))
        )
        sns.heatmap(
            lfc_plot, cmap="RdBu_r", center=0, ax=ax,
            xticklabels=True, yticklabels=True, linewidths=0.2,
            cbar_kws={"label": f"log2FC ({CONDITION_CNTRL} / {CONDITION_TEST})"},
        )
        ax.set_title(
            f"Differential gene program activity: "
            f"{CONDITION_CNTRL} vs {CONDITION_TEST}\n"
            f"(sample-level pseudobulk, exploratory)", fontsize=10)
        ax.set_ylabel("Niche")
        ax.set_xlabel("Gene Program")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
        plt.setp(ax.get_yticklabels(), fontsize=8)
        plt.tight_layout()
        fig.savefig(os.path.join(gp_dir, "gp_differential_heatmap.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close(fig)

    print(f"\n[GP] Gene program analysis complete -> {gp_dir}")


# =============================================================================
# STEP 5: Save results
# =============================================================================

def save_results(adata, model, output_dir):
    """
    Save all outputs in organized subdirectories:
        output_dir/
            objects/         — h5ad and model
            tables/          — CSV summary tables
            plots/
                umap/        — UMAP embeddings
                spatial/     — per-sample spatial niche maps
                composition/ — niche composition heatmaps/bars
    """
    import matplotlib.pyplot as plt

    # --- Create directory structure ---
    dirs = {
        "objects":     os.path.join(output_dir, "objects"),
        "tables":      os.path.join(output_dir, "tables"),
        "plots_umap":  os.path.join(output_dir, "plots", "umap"),
        "plots_spatial": os.path.join(output_dir, "plots", "spatial"),
        "plots_composition": os.path.join(output_dir, "plots", "composition"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    samples = adata.obs["sample_id"].cat.categories.tolist()

    # ---- OBJECTS ----
    adata_path = os.path.join(dirs["objects"], "nichecompass_integrated.h5ad")
    print(f"\n[SAVE] Writing {adata_path}")
    adata.write_h5ad(adata_path)

    model_dir = os.path.join(dirs["objects"], f"model_{timestamp}")
    if model is not None:
        print(f"[SAVE] Saving model to {model_dir}")
        model.save(dir_path=model_dir, overwrite=True, save_adata=False)
    else:
        print(f"[SAVE] model=None (plot-only mode) — skipping model save")

    # ---- TABLES ----
    # Per-cell niche assignments
    niche_df = adata.obs[
        ["sample_id", "condition", CELL_TYPE_KEY, LATENT_CLUSTER_KEY]
    ].copy()
    niche_csv = os.path.join(dirs["tables"], "niche_assignments.csv")
    niche_df.to_csv(niche_csv)
    print(f"[SAVE] {niche_csv}")

    # Resolution sweep results
    if "leiden_resolution_sweep" in adata.uns:
        sweep_path = os.path.join(dirs["tables"], "leiden_resolution_sweep.csv")
        adata.uns["leiden_resolution_sweep"].to_csv(sweep_path, index=False)
        print(f"[SAVE] {sweep_path}")
        best_res = adata.uns.get("leiden_best_resolution", "unknown")
        print(f"  Best Leiden resolution: {best_res}")

    # Niche x cell type composition (row-normalized)
    ct_by_niche = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs[CELL_TYPE_KEY],
        normalize="index",
    )
    ct_by_niche.to_csv(os.path.join(dirs["tables"], "niche_composition.csv"))

    # Niche x sample frequency (column-normalized)
    niche_by_sample = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs["sample_id"],
        normalize="columns",
    )
    niche_by_sample.to_csv(os.path.join(dirs["tables"], "niche_sample_frequency.csv"))

    # Niche x condition frequency (column-normalized)
    niche_by_cond = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs["condition"],
        normalize="columns",
    )
    niche_by_cond.to_csv(os.path.join(dirs["tables"], "niche_condition_frequency.csv"))

    # Niche cell counts (absolute)
    niche_counts = pd.crosstab(
        adata.obs[LATENT_CLUSTER_KEY],
        adata.obs["sample_id"],
    )
    niche_counts.to_csv(os.path.join(dirs["tables"], "niche_cell_counts.csv"))

    print(f"[SAVE] All tables -> {dirs['tables']}")

    # ---- PLOTS: UMAP ----
    print("[PLOT] Generating UMAP plots...")

    # UMAP colored by niche
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=LATENT_CLUSTER_KEY, ax=ax, show=False,
               title="NicheCompass Niches", legend_loc="on data")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_niches.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # UMAP colored by sample
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="sample_id", ax=ax, show=False,
               title="Samples in Latent Space")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_samples.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # UMAP colored by condition
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="condition", ax=ax, show=False,
               title="Condition in Latent Space")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_condition.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # UMAP colored by cell type
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.umap(adata, color=CELL_TYPE_KEY, ax=ax, show=False,
               title="Cell Types in Latent Space")
    fig.savefig(os.path.join(dirs["plots_umap"], "umap_cell_types.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    print(f"  UMAP plots -> {dirs['plots_umap']}")

    # ---- PLOTS: SPATIAL per sample ----
    print("[PLOT] Generating spatial niche maps per sample...")

    for sample in samples:
        adata_sub = adata[adata.obs["sample_id"] == sample].copy()
        condition = adata_sub.obs["condition"].iloc[0]

        # Spatial scatter colored by niche
        fig, ax = plt.subplots(figsize=(10, 10))
        coords = adata_sub.obsm[SPATIAL_KEY]
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=adata_sub.obs[LATENT_CLUSTER_KEY].cat.codes,
            cmap="tab20", s=1, alpha=0.8, rasterized=True,
        )
        ax.set_title(f"{sample} ({condition}) — Niches")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        plt.colorbar(scatter, ax=ax, label="Niche")
        fig.savefig(
            os.path.join(dirs["plots_spatial"], f"{sample}_niches.pdf"),
            bbox_inches="tight", dpi=150,
        )
        plt.close(fig)

        # Spatial scatter colored by cell type
        fig, ax = plt.subplots(figsize=(10, 10))
        ct_codes = adata_sub.obs[CELL_TYPE_KEY].cat.codes
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=ct_codes, cmap="tab20", s=1, alpha=0.8, rasterized=True,
        )
        ax.set_title(f"{sample} ({condition}) — Cell Types")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        fig.savefig(
            os.path.join(dirs["plots_spatial"], f"{sample}_cell_types.pdf"),
            bbox_inches="tight", dpi=150,
        )
        plt.close(fig)

    print(f"  Spatial plots -> {dirs['plots_spatial']}")

    # ---- PLOTS: COMPOSITION ----
    print("[PLOT] Generating composition plots...")

    # Heatmap: niche x cell type composition
    fig, ax = plt.subplots(figsize=(max(12, ct_by_niche.shape[1] * 0.5),
                                    max(6, ct_by_niche.shape[0] * 0.4)))
    import seaborn as sns
    sns.heatmap(ct_by_niche, cmap="YlOrRd", ax=ax, linewidths=0.5,
                xticklabels=True, yticklabels=True)
    ax.set_title("Cell Type Composition per Niche")
    ax.set_ylabel("Niche")
    ax.set_xlabel("Cell Type")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    fig.savefig(
        os.path.join(dirs["plots_composition"], "niche_celltype_heatmap.pdf"),
        bbox_inches="tight", dpi=150,
    )
    plt.close(fig)

    # Stacked bar: niche proportions per sample, grouped by condition
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for idx, cond in enumerate([CONDITION_CNTRL, CONDITION_TEST]):
        cond_samples = [s for s in samples if adata.obs.loc[
            adata.obs["sample_id"] == s, "condition"].iloc[0] == cond]
        if not cond_samples:
            continue
        cond_freq = niche_by_sample[cond_samples]
        cond_freq.T.plot(kind="bar", stacked=True, ax=axes[idx],
                         colormap="tab20", legend=False)
        axes[idx].set_title(f"{cond.capitalize()} wounds")
        axes[idx].set_xlabel("Sample")
        axes[idx].set_ylabel("Niche proportion")
        axes[idx].tick_params(axis="x", rotation=45)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Niche", bbox_to_anchor=(1.02, 0.5),
               loc="center left", fontsize=8)
    fig.suptitle("Niche Proportions per Sample by Condition", y=1.02)
    fig.savefig(
        os.path.join(dirs["plots_composition"], "niche_proportions_by_condition.pdf"),
        bbox_inches="tight", dpi=150,
    )
    plt.close(fig)

    # Bar: niche proportions averaged per condition (side-by-side)
    fig, ax = plt.subplots(figsize=(max(8, niche_by_cond.shape[0] * 0.6), 5))
    niche_by_cond.plot(kind="bar", ax=ax)
    ax.set_title("Niche Proportions by Condition")
    ax.set_xlabel("Niche")
    ax.set_ylabel("Proportion of Cells")
    ax.legend(title="Condition")
    plt.xticks(rotation=0)
    fig.savefig(
        os.path.join(dirs["plots_composition"], "niche_proportions_condition_comparison.pdf"),
        bbox_inches="tight", dpi=150,
    )
    plt.close(fig)

    print(f"  Composition plots -> {dirs['plots_composition']}")
    print(f"\n[DONE] All outputs in {output_dir}")


# =============================================================================
# PLOT-ONLY: reload saved h5ad and re-run analysis + plotting
# =============================================================================

def patch_condition_labels(adata, condition_map, mode):
    """
    Re-derive obs['condition'] from obs['sample_id'] using the supplied map
    and mode ('prefix' or 'suffix'). Updates in-place and returns adata.
    Useful when the original run used wrong condition labels.
    """
    import re

    def extract_token(sample_id, mode):
        if mode == "suffix":
            parts = re.split(r"[-_]", sample_id)
            return parts[-1] if parts else sample_id
        else:  # prefix
            return sample_id.split("_")[0]

    new_labels = []
    for sid in adata.obs["sample_id"].astype(str):
        token = extract_token(sid, mode)
        label = condition_map.get(token)
        if label is None:
            raise ValueError(
                f"sample_id '{sid}' token '{token}' not found in "
                f"condition_map {condition_map} (mode='{mode}'). "
                f"Check --condition_map or CONDITION_MAP_MODE."
            )
        new_labels.append(label)

    adata.obs["condition"] = pd.Categorical(new_labels)
    counts = adata.obs["condition"].value_counts().to_dict()
    print(f"  [PATCH] condition labels updated: {counts}")
    return adata


def plot_only_from_h5ad(h5ad_path, output_dir, condition_map, mode):
    """
    Load a previously saved nichecompass_integrated.h5ad, optionally patch
    condition labels, then re-run analyze_gene_programs and save_results
    (tables + all plots). Model object is not needed for either.
    """
    print(f"\n[PLOT-ONLY] Loading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(adata)

    # Validate required keys are present
    for k in [SPATIAL_KEY]:
        if k not in adata.obsm:
            raise KeyError(f"Missing obsm['{k}'] in {h5ad_path}")
    for k in ["sample_id", "condition", CELL_TYPE_KEY, LATENT_CLUSTER_KEY]:
        if k not in adata.obs.columns:
            raise KeyError(f"Missing obs['{k}'] in {h5ad_path}")

    # Patch condition labels using the current map
    print(f"\n[PLOT-ONLY] Patching condition labels (mode='{mode}')...")
    adata = patch_condition_labels(adata, condition_map, mode)

    # Ensure categorical types expected downstream
    adata.obs["sample_id"]   = adata.obs["sample_id"].astype("category")
    adata.obs["condition"]   = adata.obs["condition"].astype("category")
    adata.obs[CELL_TYPE_KEY] = adata.obs[CELL_TYPE_KEY].astype("category")
    adata.obs[LATENT_CLUSTER_KEY] = adata.obs[LATENT_CLUSTER_KEY].astype("category")

    # UMAP must exist for save_results; recompute if missing
    if "X_umap" not in adata.obsm:
        print("  [PLOT-ONLY] X_umap not found — recomputing UMAP from latent space...")
        if LATENT_KEY not in adata.obsm:
            raise KeyError(f"Missing obsm['{LATENT_KEY}'] — cannot recompute UMAP")
        sc.pp.neighbors(adata, use_rep=LATENT_KEY, key_added=LATENT_KEY)
        sc.tl.umap(adata, neighbors_key=LATENT_KEY)

    # model=None is safe: save_results only calls model.save(), which we skip
    print("\n[PLOT-ONLY] Re-running gene program analysis...")
    analyze_gene_programs(adata, model=None, output_dir=output_dir)

    print("\n[PLOT-ONLY] Re-running save / plot pipeline...")
    save_results(adata, model=None, output_dir=output_dir)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="NicheCompass niche analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
# Full training run (default):
  python run_nichecompass_integrated.py

# Re-plot from a saved h5ad with corrected lesion/diabetic labels:
  python run_nichecompass_integrated.py --plot_only \\
      --h5ad nichecompass_results/objects/nichecompass_integrated.h5ad

# Override the condition map from the command line (suffix mode, L=lesion D=diabetic):
  python run_nichecompass_integrated.py --plot_only \\
      --condition_map L=lesion D=diabetic --condition_map_mode suffix

# Prefix mode example (Cntrl=lesion Test=diabetic):
  python run_nichecompass_integrated.py --plot_only \\
      --condition_map Cntrl=lesion Test=diabetic --condition_map_mode prefix
""",
    )

    # ── Paths ─────────────────────────────────────────────────────────────────
    parser.add_argument(
        "--base_dir", default="./annotated_data",
        help="Directory containing sample subfolders with .h5ad files "
             "(full run only)",
    )
    parser.add_argument(
        "--output_dir", default="./nichecompass_results",
        help="Output directory for results",
    )
    parser.add_argument(
        "--leiden_resolution", type=float, default=None,
        help="Override Leiden resolution (full run only; default: auto-select)",
    )

    # ── Plot-only mode ────────────────────────────────────────────────────────
    parser.add_argument(
        "--plot_only", action="store_true",
        help="Skip training: load an existing h5ad, patch condition labels, "
             "and re-run all analysis and plotting.",
    )
    parser.add_argument(
        "--h5ad",
        default="./nichecompass_results/objects/nichecompass_integrated.h5ad",
        help="Path to the saved h5ad to load in --plot_only mode",
    )

    # ── Condition map overrides ───────────────────────────────────────────────
    parser.add_argument(
        "--condition_map", nargs="+", metavar="TOKEN=LABEL",
        help="Override CONDITION_MAP. Pass space-separated TOKEN=LABEL pairs, "
             "e.g. --condition_map L=lesion D=diabetic",
    )
    parser.add_argument(
        "--condition_map_mode", choices=["prefix", "suffix"], default=None,
        help="Override CONDITION_MAP_MODE ('prefix' or 'suffix'). "
             "Prefix splits on first underscore; suffix splits on last hyphen/underscore.",
    )

    args = parser.parse_args()

    # ── Apply any CLI overrides to the module-level config ───────────────────
    # Declare globals first so CLI overrides can be written back for full runs
    global CONDITION_MAP, CONDITION_MAP_MODE

    condition_map = dict(CONDITION_MAP)   # start from script defaults
    if args.condition_map:
        try:
            condition_map = dict(pair.split("=", 1) for pair in args.condition_map)
        except ValueError:
            parser.error(
                "--condition_map pairs must be in TOKEN=LABEL format, "
                "e.g. --condition_map L=lesion D=diabetic"
            )

    mode = args.condition_map_mode or CONDITION_MAP_MODE

    print("=" * 70)
    print("NicheCompass Pipeline")
    print("=" * 70)
    print(f"  Condition map : {condition_map}  (mode: {mode})")

    # ── Branch: plot-only vs full training ───────────────────────────────────
    if args.plot_only:
        print(f"  Mode          : plot-only (loading {args.h5ad})")
        plot_only_from_h5ad(
            h5ad_path=args.h5ad,
            output_dir=args.output_dir,
            condition_map=condition_map,
            mode=mode,
        )
        return

    print(f"  Mode          : full training run")

    # Step 1: Discover and load samples
    # Write CLI overrides back to module globals so discover/load pick them up
    CONDITION_MAP      = condition_map
    CONDITION_MAP_MODE = mode

    sample_info_list = discover_h5ad_files(args.base_dir)
    adatas = [load_and_prepare_sample(s) for s in sample_info_list]
    adata  = concatenate_samples(adatas)
    del adatas

    # Step 2: Build gene programs
    adata = build_gene_programs(adata)

    # Step 3: Train model
    model = train_nichecompass(adata)

    # Step 4: Cluster niches
    if args.leiden_resolution is not None:
        adata = cluster_niches(model, resolution=args.leiden_resolution)
    else:
        adata = cluster_niches(model)

    # Step 4b: Gene program differential analysis
    analyze_gene_programs(adata, model, args.output_dir)

    # Step 5: Save results
    save_results(adata, model, args.output_dir)


if __name__ == "__main__":
    main()