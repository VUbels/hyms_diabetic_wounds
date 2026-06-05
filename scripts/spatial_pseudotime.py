"""
==============================================================================
WOUND HEALING XENIUM 5K PIPELINE — STEP 3 (updated for NicheCompass output)
Pseudotime Construction: Spatial + NicheCompass-latent Hybrid
==============================================================================

Inputs:
    nichecompass_results/objects/nichecompass_integrated.h5ad

Object structure assumed (from run_nichecompass_integrated.py):
    adata.obsm['X_spatial']              : per-cell spatial coordinates
    adata.obsm['nichecompass_latent']    : spatially-aware latent embedding
    adata.obs['predicted_cell_type']     : cell annotations
    adata.obs['nichecompass_niche']      : niche labels (numeric Leiden IDs)
    adata.obs['sample_id']               : slide identifier
    adata.obs['condition']               : 'HS' or 'LS'
    adata.obs['wound_region']            : pre-annotated wound edge/bed/intact
    adata.layers['counts']               : raw counts (float32)

Outputs:
    results_step3/adata_step3_pseudotime.h5ad
    results_step3/pseudotime_validation.pdf
    results_step3/pseudotime_spatial_per_slide.pdf
    results_step3/pseudotime_decoupling_scatter.pdf

Conceptual note on using NicheCompass latent for DPT:
    The NicheCompass embedding encodes both cell-intrinsic gene expression
    and spatial niche context. DPT on this latent therefore gives you a
    "context-aware" transcriptomic pseudotime rather than a purely cell-
    intrinsic one. This is biologically appropriate for wound healing where
    niche context drives fate, but frame findings accordingly.
==============================================================================
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

# -----------------------------------------------------------------------------
# CONFIG — matches NicheCompass script's key names
# -----------------------------------------------------------------------------
INPUT_H5AD    = "nichecompass_results/objects/nichecompass_integrated.h5ad"
OUTPUT_DIR    = "pseudotime"
SPATIAL_KEY   = "X_spatial"
LATENT_KEY    = "nichecompass_latent"
CELLTYPE_KEY  = "predicted_cell_type"
NICHE_KEY     = "nichecompass_niche"
SAMPLE_KEY    = "sample_id"
CONDITION_KEY = "condition"

# Candidate column names for wound region — script will try each
WOUND_REGION_CANDIDATES = [
    "wound_region", "region", "wound_annotation", "annotation",
    "manual_annotation", "histology", "zone",
]

# Expected labels for wound compartments — substring matching, case-insensitive
WOUND_EDGE_LABELS = ["wound_edge", "edge", "margin", "leading_edge"]
WOUND_BED_LABELS  = ["wound_bed", "bed", "granulation_tissue", "wound_center",
                     "center", "inflammatory_core", "hypoxic_core"]
INTACT_LABELS     = ["intact", "intact_dermis", "normal_skin",
                     "uninjured", "distal"]

# Pseudotime weighting — tune after inspecting validation plots
SPATIAL_WEIGHT = 0.5

# Validation marker lists — canonical dermal wound healing dynamics
VALIDATION_MARKERS = {
    "early_active":    ["KRT6A", "KRT16", "KRT17", "S100A8", "S100A9",
                         "IL1B", "TNF", "IL6", "CXCL1", "CXCL8"],
    "late_resolved":   ["KRT1", "KRT10", "FLG", "IVL", "LOR",
                         "COL1A1", "COL3A1", "ARG1", "MRC1", "CD163"],
    "proliferation":   ["MKI67", "TOP2A", "PCNA", "CCNB1", "CCND1"],
    "hypoxia":         ["HIF1A", "VEGFA", "HILPDA", "BNIP3", "LDHA",
                         "CA9", "SLC2A1"],
    "senescence_sasp": ["CDKN1A", "CDKN2A", "GLB1", "SERPINE1",
                         "MMP1", "MMP3", "IL8"],
    "angiogenesis":    ["ESM1", "APLN", "PECAM1", "VWF", "CDH5"],
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. LOAD AND INSPECT
# -----------------------------------------------------------------------------
print(f"[LOAD] {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

print(f"  {adata.n_obs} cells, {adata.n_vars} genes")
print(f"  Samples: {adata.obs[SAMPLE_KEY].value_counts().to_dict()}")
print(f"  Conditions: {adata.obs[CONDITION_KEY].value_counts().to_dict()}")
print(f"  Niches: {adata.obs[NICHE_KEY].nunique()}")
print(f"  Cell types: {adata.obs[CELLTYPE_KEY].nunique()}")

# Find wound_region column
wound_region_col = None
for cand in WOUND_REGION_CANDIDATES:
    if cand in adata.obs.columns:
        wound_region_col = cand
        print(f"  Found wound region annotation: obs['{cand}']")
        print(f"    Values: {adata.obs[cand].value_counts().to_dict()}")
        break

if wound_region_col is None:
    print("\n[ERROR] No wound_region column found. Available obs columns:")
    for c in adata.obs.columns:
        print(f"    - {c}")
    raise ValueError(
        "Set wound_region_col manually at the top of this script, or add "
        "the annotation column to the object."
    )

# -----------------------------------------------------------------------------
# 2. SPATIAL PSEUDOTIME — SIGNED DISTANCE FROM WOUND EDGE
# -----------------------------------------------------------------------------
def classify_wound_label(label, edge_labels, bed_labels, intact_labels):
    """Return 'edge', 'bed', 'intact', or None using case-insensitive substring match."""
    if pd.isna(label):
        return None
    ll = str(label).lower().strip().replace(" ", "_")
    if any(e.lower() in ll for e in edge_labels):
        return "edge"
    if any(b.lower() in ll for b in bed_labels):
        return "bed"
    if any(i.lower() in ll for i in intact_labels):
        return "intact"
    return None

adata.obs["_wound_class"] = adata.obs[wound_region_col].apply(
    lambda x: classify_wound_label(x, WOUND_EDGE_LABELS, WOUND_BED_LABELS,
                                     INTACT_LABELS)
)
print("\n  Wound class counts:",
      adata.obs["_wound_class"].value_counts(dropna=False).to_dict())

if adata.obs["_wound_class"].value_counts().get("edge", 0) == 0:
    raise ValueError(
        f"No cells classified as 'edge' from obs['{wound_region_col}']. "
        f"Unique values: {adata.obs[wound_region_col].unique()}. "
        f"Edit WOUND_EDGE_LABELS to match your annotation scheme."
    )

def compute_signed_wound_distance(adata, spatial_key, sample_key):
    """Per-cell signed distance (μm) to nearest edge cell, computed per slide."""
    distances = np.full(adata.n_obs, np.nan)
    for slide in adata.obs[sample_key].unique():
        mask = (adata.obs[sample_key] == slide).values
        idx = np.where(mask)[0]
        coords = adata.obsm[spatial_key][mask]
        wclass = adata.obs["_wound_class"].values[mask]

        edge_mask = wclass == "edge"
        if edge_mask.sum() == 0:
            print(f"  [WARN] {slide}: no edge cells; skipping")
            continue

        tree = cKDTree(coords[edge_mask])
        d, _ = tree.query(coords, k=1)

        # Positive = in wound bed, negative = in intact tissue
        bed_mask = wclass == "bed"
        signed = d.copy()
        signed[~bed_mask] = -signed[~bed_mask]
        distances[idx] = signed
    return distances

adata.obs["wound_distance_um"] = compute_signed_wound_distance(
    adata, SPATIAL_KEY, SAMPLE_KEY
)

# Rank-normalize to [0, 1]: 0 = deepest intact, 1 = deepest wound bed
valid = ~np.isnan(adata.obs["wound_distance_um"])
sp_pt = np.full(adata.n_obs, np.nan)
sp_pt[valid] = pd.Series(
    adata.obs["wound_distance_um"].values[valid]
).rank(pct=True).values
adata.obs["spatial_pseudotime"] = sp_pt

# -----------------------------------------------------------------------------
# 3. TRANSCRIPTOMIC PSEUDOTIME — DPT ON NICHECOMPASS LATENT
# -----------------------------------------------------------------------------
# Important: the NicheCompass script stored the neighbor graph under the key
# 'nichecompass_latent', not the default. We rebuild under the default key
# here so CellRank 2 in Step 4 picks it up without configuration.
print("\n[DPT] Building default-keyed neighbor graph on NicheCompass latent...")
sc.pp.neighbors(adata, use_rep=LATENT_KEY, n_neighbors=30)

sc.tl.diffmap(adata, n_comps=15)

# Root: highest quiescence score among intact-dermis cells
quiescence_markers = ["PI16", "DPT", "CD34", "COL17A1", "TP63", "KRT15"]
quiescence_markers = [g for g in quiescence_markers if g in adata.var_names]

if not quiescence_markers:
    print("  [WARN] No quiescence markers in panel; using first intact cell")
    intact_mask = (adata.obs["_wound_class"] == "intact").values
    adata.uns["iroot"] = int(np.where(intact_mask)[0][0]) if intact_mask.any() else 0
else:
    print(f"  Quiescence markers used: {quiescence_markers}")
    sc.tl.score_genes(adata, quiescence_markers, score_name="quiescence_score")
    intact_mask = (adata.obs["_wound_class"] == "intact").values
    scores = adata.obs["quiescence_score"].values.copy().astype(float)
    if intact_mask.any():
        scores[~intact_mask] = -np.inf
    adata.uns["iroot"] = int(np.argmax(scores))
    print(f"  Root cell idx: {adata.uns['iroot']} "
          f"in sample {adata.obs[SAMPLE_KEY].iloc[adata.uns['iroot']]}")

sc.tl.dpt(adata)
adata.obs["transcriptomic_pseudotime"] = adata.obs["dpt_pseudotime"]

# -----------------------------------------------------------------------------
# 4. HYBRID PSEUDOTIME
# -----------------------------------------------------------------------------
def rank_normalize(x):
    x = np.asarray(x, dtype=float)
    mask = ~np.isnan(x)
    out = np.full_like(x, np.nan, dtype=float)
    out[mask] = pd.Series(x[mask]).rank(pct=True).values
    return out

sp_rank = rank_normalize(adata.obs["spatial_pseudotime"].values)
tx_rank = rank_normalize(adata.obs["transcriptomic_pseudotime"].values)

# Fall back to transcriptomic only where spatial is NaN
hybrid = np.where(
    np.isnan(sp_rank),
    tx_rank,
    SPATIAL_WEIGHT * sp_rank + (1 - SPATIAL_WEIGHT) * tx_rank,
)
adata.obs["hybrid_pseudotime"] = hybrid

# -----------------------------------------------------------------------------
# 5. VALIDATION: MARKER DYNAMICS IN HEALTHY
# -----------------------------------------------------------------------------
print("\n[VALIDATE] Marker dynamics along pseudotime (healthy only)")
healthy_mask = (adata.obs[CONDITION_KEY] == "healthy").values

def get_expr(adata_sub, gene):
    X = adata_sub[:, gene].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X).flatten()

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
for ax, (group, genes) in zip(axes.flat, VALIDATION_MARKERS.items()):
    genes_present = [g for g in genes if g in adata.var_names]
    if not genes_present:
        ax.set_title(f"{group}\n(no markers in panel)")
        continue

    adata_h = adata[healthy_mask]
    pt = adata_h.obs["hybrid_pseudotime"].values
    bins = np.linspace(0, 1, 25)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_idx = np.digitize(pt, bins) - 1

    for gene in genes_present:
        expr = get_expr(adata_h, gene)
        means, centers = [], []
        for b in range(len(bin_centers)):
            m = bin_idx == b
            if m.sum() > 20:
                means.append(expr[m].mean())
                centers.append(bin_centers[b])
        if means:
            ax.plot(centers, means, label=gene, alpha=0.75, linewidth=2)

    ax.set_xlabel("Hybrid pseudotime (0 = intact, 1 = wound bed)")
    ax.set_ylabel("Log-norm expression")
    ax.set_title(group)
    ax.legend(fontsize=7, loc="best")

plt.suptitle("Pseudotime validation (healthy)\n"
             "Early/active should ↑ toward wound bed; resolved should peak away from it",
             fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "pseudotime_validation.pdf"),
            bbox_inches="tight")
plt.close()

# -----------------------------------------------------------------------------
# 6. SPATIAL-TRANSCRIPTOMIC DECOUPLING ANALYSIS
# -----------------------------------------------------------------------------
# This is the core quantitative readout for your trajectory divergence question.
# If healthy shows high spatial-transcriptomic coupling but chronic shows
# decoupling, that is direct evidence chronic cells are spatially positioned
# for healing but transcriptomically stalled.
print("\n[DECOUPLE] Spatial vs transcriptomic pseudotime correlation:")
valid_both = (~np.isnan(adata.obs["spatial_pseudotime"]) &
              ~np.isnan(adata.obs["transcriptomic_pseudotime"]))

decoupling_results = []
for cond in ["healthy", "chronic"]:
    cond_mask = ((adata.obs[CONDITION_KEY] == cond) & valid_both).values
    if cond_mask.sum() < 100:
        continue
    r, p = spearmanr(
        adata.obs.loc[cond_mask, "spatial_pseudotime"],
        adata.obs.loc[cond_mask, "transcriptomic_pseudotime"],
    )
    decoupling_results.append({
        "condition": cond, "n_cells": int(cond_mask.sum()),
        "spearman_r": r, "p_value": p,
    })
    print(f"  {cond}: Spearman r = {r:.3f} (n = {cond_mask.sum()})")

pd.DataFrame(decoupling_results).to_csv(
    os.path.join(OUTPUT_DIR, "pseudotime_decoupling_by_condition.csv"),
    index=False,
)

# Also compute decoupling PER sample, so you can see if the effect is consistent
print("\n[DECOUPLE] Per-sample correlation:")
per_sample = []
for sample in adata.obs[SAMPLE_KEY].unique():
    m = ((adata.obs[SAMPLE_KEY] == sample) & valid_both).values
    if m.sum() < 100:
        continue
    r, p = spearmanr(
        adata.obs.loc[m, "spatial_pseudotime"],
        adata.obs.loc[m, "transcriptomic_pseudotime"],
    )
    cond = adata.obs.loc[m, CONDITION_KEY].iloc[0]
    per_sample.append({"sample": sample, "condition": cond,
                       "n_cells": int(m.sum()), "spearman_r": r, "p": p})
    print(f"  {sample} ({cond}): r = {r:.3f}")

pd.DataFrame(per_sample).to_csv(
    os.path.join(OUTPUT_DIR, "pseudotime_decoupling_per_sample.csv"),
    index=False,
)

# Scatter visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
for ax, cond in zip(axes, ["healthy", "chronic"]):
    m = ((adata.obs[CONDITION_KEY] == cond) & valid_both).values
    if m.sum() == 0:
        continue
    idx = np.where(m)[0]
    if len(idx) > 20000:
        idx = np.random.choice(idx, 20000, replace=False)
    ax.hexbin(
        adata.obs["spatial_pseudotime"].values[idx],
        adata.obs["transcriptomic_pseudotime"].values[idx],
        gridsize=50, cmap="viridis", mincnt=1,
    )
    rvals = [d["spearman_r"] for d in decoupling_results if d["condition"] == cond]
    ax.set_title(f"{cond} (Spearman r = {rvals[0]:.3f})" if rvals else cond)
    ax.set_xlabel("Spatial pseudotime")
    ax.set_ylabel("Transcriptomic pseudotime")
    ax.plot([0, 1], [0, 1], "r--", alpha=0.5, linewidth=1)

plt.suptitle("Spatial vs transcriptomic pseudotime decoupling by condition\n"
             "Red dashed = perfect coupling", fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "pseudotime_decoupling_scatter.pdf"),
            bbox_inches="tight")
plt.close()

# -----------------------------------------------------------------------------
# 7. PER-SLIDE SPATIAL PSEUDOTIME VISUALIZATION
# -----------------------------------------------------------------------------
samples = (adata.obs[SAMPLE_KEY].cat.categories.tolist()
           if hasattr(adata.obs[SAMPLE_KEY], "cat")
           else sorted(adata.obs[SAMPLE_KEY].unique()))

n_samples = len(samples)
n_cols = 3
n_rows = int(np.ceil(n_samples / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
axes = np.atleast_1d(axes).flatten()

for ax, sample in zip(axes, samples):
    m = (adata.obs[SAMPLE_KEY] == sample).values
    coords = adata.obsm[SPATIAL_KEY][m]
    pt = adata.obs["hybrid_pseudotime"].values[m]
    cond = adata.obs[CONDITION_KEY].values[m][0]
    s = ax.scatter(coords[:, 0], coords[:, 1], c=pt, cmap="viridis",
                   s=0.5, alpha=0.7, rasterized=True, vmin=0, vmax=1)
    ax.set_title(f"{sample} ({cond})", fontsize=10)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(s, ax=ax, fraction=0.04, label="pseudotime")

for ax in axes[n_samples:]:
    ax.axis("off")

plt.suptitle("Hybrid pseudotime in space, per slide", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "pseudotime_spatial_per_slide.pdf"),
            bbox_inches="tight", dpi=150)
plt.close()

# -----------------------------------------------------------------------------
# 8. SAVE
# -----------------------------------------------------------------------------
out_path = os.path.join(OUTPUT_DIR, "adata_step3_pseudotime.h5ad")
adata.write_h5ad(out_path)
print(f"\n[SAVE] {out_path}")
print("\n[DONE] Step 3 complete. Inspect:")
print(f"  - {OUTPUT_DIR}/pseudotime_validation.pdf")
print(f"  - {OUTPUT_DIR}/pseudotime_decoupling_scatter.pdf")
print(f"  - {OUTPUT_DIR}/pseudotime_spatial_per_slide.pdf")
print("\nIf validation plots look wrong, adjust:")
print("  - Root cell (uns['iroot'])")
print("  - SPATIAL_WEIGHT (try 0.3 or 0.7)")
print("  - WOUND_*_LABELS to match your annotations")