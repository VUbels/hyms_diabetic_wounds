#!/usr/bin/env python3
"""
reconstruct_skin_reference.py
=============================================================================
Reconstructs the perinatal mouse-skin scRNA reference from Lee et al.
(Exp Mol Med 2026, GSE286328) following the published pipeline in
`#####FigureS1_integration.ipynb`.

This reproduces the integration + Leiden(res=0.9) clustering so that the
per-cluster marker genes can be compared against the paper. The cluster ->
cell-type *labels* are NOT in the repo (the marker-study notebook loads an
already-labelled object), so the final naming step is left to you. The exact
34-label vocabulary and the 64-gene marker panel they annotated against are
embedded at the bottom of this file for that purpose.

SCOPE NOTE (read this):
  The published reference integrates E13.5, E16.5, E18.5, PD0, PD2, PD4.
  The three public datasets you have cover E13.5, E16.5, PD0, PD2, PD4 only.
  E18.5 scRNA lives in the new GEO (GSE286328: GSM8723662/3) and is NOT
  included here. Omitting it means:
    - total cell count will not match the paper's 67,558,
    - the Leiden(0.9) partition and a few perinatal-transition labels may
      shift, and the DP-lineage subcluster index will differ.
  Hooks to add E18.5 later are marked with `ADD_E18` below.

Faithfulness vs. the notebook:
  - QC filtering uses the notebook's HARD thresholds (cells 18-21), NOT the
    MAD approach described in the Methods prose (MAD was the upstream
    per-sample step that produced their *_quality_control.h5ad inputs).
  - SoupX (ambient RNA) is in the Methods but is NOT in the integration
    notebook and cannot be auto-run from GEO *filtered* matrices (no raw
    droplets). It is skipped here; this minimally affects cluster identity
    relative to the scVI integration. Flag: RUN_SOUPX (left False).
  - scDblFinder doublet removal IS replicated (per sample, via rpy2->R).
    Falls back to scrublet if R/scDblFinder is unavailable (a minor deviation).

Dependencies:
  scanpy, anndata, scvi-tools (paper used 0.20.3), scikit-misc (for the
  seurat_v3 HVG flavor), leidenalg. Optional: rpy2 + anndata2ri + R pkg
  scDblFinder (for faithful doublet calls).
=============================================================================
"""

import os
import gzip
import glob
import tarfile
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import io as sio
from scipy import sparse

import scvi

# ----------------------------------------------------------------------------
# 0. CONFIG  -- edit these paths/flags, nothing else should need changing
# ----------------------------------------------------------------------------
DATA_DIR = r"/mnt/d/scRNA_datasets/HYMS_metal_diabetes/Reference"  # WSL2 mount of D:\ (use /mnt/d/..., not D:\...)
WORK_DIR = os.path.join(DATA_DIR, "_extracted")                 # tars get unpacked here
OUT_H5AD = os.path.join(DATA_DIR, "reconstructed_skin_reference.h5ad")

GSE122043_TAR = os.path.join(DATA_DIR, "GSE122043_RAW.tar")          # E13.5
GSE122043_GENES = os.path.join(DATA_DIR, "GSE122043_genes.tsv.gz")  # shared genes file (download
#   separately from GEO if it is not already next to the tar; GSE122043 stores
#   one genes.tsv for all samples rather than one per sample).
GSE131498_TXT = os.path.join(DATA_DIR, "GSE131498_scRNA_seq_skin.txt.gz")  # E16.5 dense matrix
GSE181390_TAR = os.path.join(DATA_DIR, "GSE181390_RAW.tar")              # PD0/PD2/PD4

SCVI_SEED = 0          # paper's model was "...seed0..."
SCVI_MAX_EPOCHS = None # None -> scvi default (matches notebook `vae.train()`)
RUN_SOUPX = False      # see note above; leave False for the public filtered matrices
RUN_SCDBLFINDER = True # faithful doublet removal; falls back to scrublet on failure
# If rpy2 cannot auto-detect R, set this to the output of `R RHOME`. Leave None when
# R lives in the same conda env (recommended) -- do NOT reuse the notebook's hardcoded path.
R_HOME = None
GSE131498_GENES_ARE_ROWS = True  # dense matrix orientation; flip if load looks transposed

# GSM -> (orig.ident, bulk.ident). orig.ident is the scVI per-sample batch;
# bulk.ident is the timepoint batch. Both are scVI categorical covariates.
PD_SAMPLES = {
    "GSM5501542": ("PD0_1", "PD0"),
    "GSM5501543": ("PD0_2", "PD0"),
    "GSM5501545": ("PD2_1", "PD2"),
    "GSM5501546": ("PD2_2", "PD2"),
    "GSM5501548": ("PD4_1", "PD4"),
    "GSM5501550": ("PD4_3", "PD4"),   # NB: paper used PD4 Rep1 + Rep3, not Rep2
}
E13_GSMS = ["GSM3453535", "GSM3453536"]  # both E13.5 controls, merged into one "E13.5" sample
MERGE_E13_REPS = True  # the notebook used a single E13 object; True merges the two control reps

# ADD_E18: to include E18.5 later, drop the two GSE286328 scRNA H5s somewhere and
# load them as samples ("E18.5_1","E18.5") and ("E18.5_2","E18.5"), then append to `samples`.

sc.settings.verbosity = 1
warnings.simplefilter("ignore", category=FutureWarning)


# ----------------------------------------------------------------------------
# 1. Loaders
# ----------------------------------------------------------------------------
def _extract_tar(tar_path, dest):
    os.makedirs(dest, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        try:
            tf.extractall(dest, filter="data")   # py3.12+: silence/secure extraction
        except TypeError:
            tf.extractall(dest)                  # older pythons lack the filter arg
    return dest


def ensure_gse122043_genes():
    """GSE122043 stores one genes.tsv for the whole series (not per sample).
    Download it next to the data if it's missing."""
    if os.path.exists(GSE122043_GENES):
        return GSE122043_GENES
    import urllib.request
    url = ("https://www.ncbi.nlm.nih.gov/geo/download/"
           "?acc=GSE122043&format=file&file=GSE122043%5Fgenes%2Etsv%2Egz")
    print(f">> GSE122043 shared genes file missing; downloading -> {GSE122043_GENES}")
    os.makedirs(os.path.dirname(GSE122043_GENES), exist_ok=True)
    urllib.request.urlretrieve(url, GSE122043_GENES)
    return GSE122043_GENES


def _find(dirpath, gsm, *suffixes):
    """Return the first file under dirpath whose name contains gsm and any suffix."""
    for suf in suffixes:
        hits = glob.glob(os.path.join(dirpath, f"*{gsm}*{suf}"))
        if hits:
            return sorted(hits)[0]
    return None


def load_10x_sample(matrix_path, barcodes_path, features_path):
    """Load a 10x MTX triplet into an AnnData (cells x genes)."""
    mat = sio.mmread(matrix_path).tocsr()           # genes x cells in 10x convention
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].values
    feats = pd.read_csv(features_path, header=None, sep="\t")
    # gene symbol is col 1 (0-based) for v2 genes.tsv and v3 features.tsv; fall back to col 0
    gene_names = feats[1].values if feats.shape[1] > 1 else feats[0].values
    if mat.shape[0] == len(gene_names) and mat.shape[1] == len(barcodes):
        X = mat.T.tocsr()                            # -> cells x genes
    elif mat.shape[0] == len(barcodes) and mat.shape[1] == len(gene_names):
        X = mat.tocsr()
    else:
        raise ValueError(f"Matrix shape {mat.shape} does not match "
                         f"{len(barcodes)} barcodes / {len(gene_names)} genes")
    a = ad.AnnData(X=X.astype(np.float32))
    a.obs_names = barcodes
    a.var_names = gene_names
    a.var_names_make_unique()
    return a


def load_gse122043_e13(workdir):
    """E13.5 from GSE122043 (10x v2; genes.tsv shared at series level)."""
    ensure_gse122043_genes()
    parts = []
    for gsm in E13_GSMS:
        mtx = _find(workdir, gsm, "matrix.mtx.gz", "matrix.mtx")
        bcs = _find(workdir, gsm, "barcodes.tsv.gz", "barcodes.tsv")
        feats = _find(workdir, gsm, "features.tsv.gz", "genes.tsv.gz",
                      "features.tsv", "genes.tsv")
        if feats is None and os.path.exists(GSE122043_GENES):
            feats = GSE122043_GENES   # shared genes file fallback
        if not (mtx and bcs and feats):
            raise FileNotFoundError(f"E13.5 files for {gsm} not found in {workdir} "
                                    f"(matrix={mtx}, barcodes={bcs}, genes={feats})")
        s = load_10x_sample(mtx, bcs, feats)
        s.obs["orig.ident"] = "E13.5"
        s.obs["bulk.ident"] = "E13.5"
        s.obs["_gsm"] = gsm
        parts.append(s)
    if MERGE_E13_REPS:
        e13 = ad.concat(parts, join="outer", index_unique="-")
        e13.obs["orig.ident"] = "E13.5"
        e13.obs["bulk.ident"] = "E13.5"
        return [e13]
    return parts


def load_gse131498_e16(reference_genes=None):
    """E16.5 from GSE131498 (single dense text matrix).

    Orientation is auto-detected by overlapping each axis against known mouse
    gene symbols (passed in from the already-loaded 10x samples), so it does
    not matter whether the file is genes x cells or cells x genes, nor what the
    cell labels look like. Falls back to a size heuristic, then to the manual
    GSE131498_GENES_ARE_ROWS flag, if the overlap is uninformative."""
    # sniff delimiter from the header line
    with gzip.open(GSE131498_TXT, "rt") as fh:
        head = fh.readline()
    if "\t" in head:
        sep = "\t"
    elif head.count(",") > head.count(" "):
        sep = ","
    else:
        sep = r"\s+"
    if sep == r"\s+":
        raw = pd.read_csv(GSE131498_TXT, sep=sep, engine="python",
                          index_col=0, compression="gzip")
    else:
        raw = pd.read_csv(GSE131498_TXT, sep=sep, index_col=0, compression="gzip")
    print(f"   [E16.5] raw matrix {raw.shape} (rows x cols); "
          f"row[:3]={list(raw.index[:3])} col[:3]={list(raw.columns[:3])}")

    genes_are_rows = None
    if reference_genes:
        ov_rows = len(set(map(str, raw.index)) & reference_genes)
        ov_cols = len(set(map(str, raw.columns)) & reference_genes)
        print(f"   [E16.5] gene-symbol overlap: rows={ov_rows}, cols={ov_cols}")
        if max(ov_rows, ov_cols) > 100:
            genes_are_rows = ov_rows >= ov_cols
    if genes_are_rows is None:                       # size-based fallback
        r_is_gene = 8000 <= raw.shape[0] <= 40000
        c_is_gene = 8000 <= raw.shape[1] <= 40000
        if r_is_gene and not c_is_gene:
            genes_are_rows = True
        elif c_is_gene and not r_is_gene:
            genes_are_rows = False
        else:
            genes_are_rows = GSE131498_GENES_ARE_ROWS
            warnings.warn("E16.5 orientation ambiguous; using "
                          f"GSE131498_GENES_ARE_ROWS={GSE131498_GENES_ARE_ROWS}. "
                          "Verify with a manual peek of the file.")

    df = raw.T if genes_are_rows else raw            # -> cells x genes
    a = ad.AnnData(X=sparse.csr_matrix(df.values.astype(np.float32)))
    a.obs_names = df.index.astype(str)
    a.var_names = df.columns.astype(str)
    a.var_names_make_unique()
    a.obs["orig.ident"] = "E16.5"
    a.obs["bulk.ident"] = "E16.5"
    if a.n_obs == 0 or a.n_vars == 0:
        raise ValueError(f"E16.5 loaded empty {a.shape}: orientation/parse wrong. "
                         "Peek the file and set GSE131498_GENES_ARE_ROWS by hand.")
    print(f"   [E16.5] -> {a.n_obs} cells x {a.n_vars} genes "
          f"(genes_are_rows={genes_are_rows})")
    # SANITY: scVI (nb) and scDblFinder need raw integer counts.
    if a.X.size and not np.allclose(a.X.data, np.round(a.X.data)):
        warnings.warn("GSE131498 values are not integers -- scVI/scDblFinder assume raw "
                      "counts. If this is normalized/log data, obtain the raw UMI matrix.")
    return a


def load_gse181390_pd(workdir):
    """PD0/PD2/PD4 WT reps from GSE181390 (10x MTX triplets)."""
    out = []
    for gsm, (orig, bulk) in PD_SAMPLES.items():
        mtx = _find(workdir, gsm, "matrix.mtx.gz", "matrix.mtx")
        bcs = _find(workdir, gsm, "barcodes.tsv.gz", "barcodes.tsv")
        feats = _find(workdir, gsm, "features.tsv.gz", "genes.tsv.gz",
                      "features.tsv", "genes.tsv")
        if not (mtx and bcs and feats):
            raise FileNotFoundError(f"{gsm} ({orig}) files not found in {workdir}")
        s = load_10x_sample(mtx, bcs, feats)
        s.obs["orig.ident"] = orig
        s.obs["bulk.ident"] = bulk
        s.obs["_gsm"] = gsm
        out.append(s)
    return out


# ----------------------------------------------------------------------------
# 2. Per-sample QC metrics + doublet removal
# ----------------------------------------------------------------------------
def annotate_qc(a):
    a.var["mt"] = a.var_names.str.startswith("mt-")
    a.var["ribo"] = a.var_names.str.startswith(("Rps", "Rpl"))
    a.var["hb"] = a.var_names.str.startswith(("Hba", "Hbb"))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo", "hb"],
                               percent_top=[20], log1p=True, inplace=True)
    return a


def run_scdblfinder(a):
    """Per-sample scDblFinder via rpy2->R; returns a boolean 'is_singlet' array.
    Converts a minimal counts-only AnnData to a SingleCellExperiment with
    anndata2ri (the paper's approach), which is more robust than hand-passing a
    sparse matrix. Falls back to scanpy scrublet if R/scDblFinder is unavailable."""
    if a.n_obs == 0:
        raise ValueError("sample has 0 cells before doublet calling -- loader produced an "
                         "empty object; fix orientation/parse upstream.")
    if a.n_obs < 50:
        warnings.warn(f"only {a.n_obs} cells; skipping doublet calling for this sample.")
        return np.ones(a.n_obs, dtype=bool)
    try:
        if R_HOME:                       # only if rpy2 can't auto-detect R
            os.environ["R_HOME"] = R_HOME
        # The project is an renv project: launched from its root, R auto-activates renv
        # and switches .libPaths() to the renv library (which lacks scDblFinder). Disable
        # the autoloader BEFORE R starts so R keeps the conda env's library paths.
        os.environ.setdefault("RENV_CONFIG_AUTOLOADER_ENABLED", "FALSE")
        import anndata2ri
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.conversion import localconverter
        # Belt-and-suspenders: explicitly put the conda R library (R.home()/library,
        # where bioconductor-scdblfinder installs) first, in case renv already ran.
        ro.r('.libPaths(c(file.path(R.home(), "library"), .libPaths()))')
        ro.r('suppressMessages({library(scDblFinder); library(SingleCellExperiment)})')
        counts = a.layers["counts"] if "counts" in a.layers else a.X
        tmp = ad.AnnData(X=counts.copy())          # minimal, counts-only -> clean conversion
        tmp.obs_names = a.obs_names.astype(str)
        tmp.var_names = a.var_names.astype(str)
        with localconverter(anndata2ri.converter):
            globalenv["sce"] = tmp                  # AnnData(cells x genes) -> SCE(genes x cells)
        ro.r('''
          set.seed(0)
          if (!("counts" %in% assayNames(sce))) { assay(sce, "counts") <- assays(sce)[[1]] }
          sce <- scDblFinder(sce)
          cls <- as.character(sce$scDblFinder.class)
        ''')
        cls = np.array(list(globalenv["cls"]))
        return cls == "singlet"
    except Exception as e:                 # noqa: BLE001
        warnings.warn(f"scDblFinder failed ({e}); trying scrublet.")
        try:
            sc.pp.scrublet(a)              # adds a.obs['predicted_doublet']
            return ~a.obs["predicted_doublet"].values
        except Exception as e2:            # noqa: BLE001
            warnings.warn(f"scrublet also failed ({e2}); keeping all cells as singlets.")
            return np.ones(a.n_obs, dtype=bool)


# ----------------------------------------------------------------------------
# 3. Pipeline (mirrors integration notebook cells 14-49)
# ----------------------------------------------------------------------------
def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    print(">> extracting tars (skip if already extracted)")
    if not glob.glob(os.path.join(WORK_DIR, "*GSM34535*")):
        _extract_tar(GSE122043_TAR, WORK_DIR)
    if not glob.glob(os.path.join(WORK_DIR, "*GSM55015*")):
        _extract_tar(GSE181390_TAR, WORK_DIR)

    print(">> loading samples")
    e13 = load_gse122043_e13(WORK_DIR)        # 10x, gene symbols -> reference vocabulary
    pd_samples = load_gse181390_pd(WORK_DIR)  # 10x, gene symbols
    refgenes = set()
    for x in e13 + pd_samples:
        refgenes |= set(map(str, x.var_names))
    e16 = load_gse131498_e16(refgenes)        # orientation auto-detected against refgenes
    samples = e13 + [e16] + pd_samples        # keep notebook order: E13, E16, PD...
    # ADD_E18: samples += load_e18(...)   # see header note

    print(">> per-sample QC + doublet removal")
    cleaned = []
    for s in samples:
        s = annotate_qc(s)
        if RUN_SCDBLFINDER:
            keep = run_scdblfinder(s)
            s = s[keep].copy()
        cleaned.append(s)
        print(f"   {s.obs['orig.ident'][0]:>8}: {s.n_obs} singlets x {s.n_vars} genes")

    print(">> concatenate (notebook: E13.concatenate(E16, ...))")
    adata = ad.concat(cleaned, join="outer", label="batch",
                      index_unique="-", fill_value=0)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # ---- QC filtering: hard thresholds, exact order of notebook cells 18-21 ----
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts > 1000, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    adata = adata[adata.obs.pct_counts_ribo > 5, :]
    adata = adata[adata.obs.pct_counts_hb < 1, :]
    adata = adata[adata.obs.n_genes_by_counts < 7000, :]
    adata = adata[adata.obs.total_counts < 40000, :]
    adata = adata[adata.obs.total_counts > 4000, :].copy()
    print(f">> post-QC: {adata.n_obs} cells x {adata.n_vars} genes "
          f"(paper had 67,558 cells WITH E18.5)")

    # ---- normalize / log1p / freeze raw (notebook cell 24) ----
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # ---- HVG: seurat_v3 on raw counts, 3000, batched by orig.ident (cell 29) ----
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True,
                                layer="counts", flavor="seurat_v3",
                                batch_key="orig.ident")

    # scVI needs plain column names (no dots) for covariate keys (cell 32)
    adata.obs["bulk_ident"] = adata.obs["bulk.ident"].astype("category")
    adata.obs["orig_ident"] = adata.obs["orig.ident"].astype("category")

    # ---- scVI integration (cell 33), identical params ----
    scvi.settings.seed = SCVI_SEED
    scvi.model.SCVI.setup_anndata(
        adata, layer="counts",
        categorical_covariate_keys=["bulk_ident", "orig_ident"],
        continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo"],
    )
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train(max_epochs=SCVI_MAX_EPOCHS, check_val_every_n_epoch=1)
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=10e4)

    # ---- neighbors / UMAP / Leiden 0.9 (cells 38, 47) ----
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    # flavor="leidenalg" reproduces the pre-scanpy-1.10 default the paper used;
    # scanpy >=1.10 switched the default to the igraph flavor (different partitions).
    sc.tl.leiden(adata, key_added="leiden_scVI_0.9", resolution=0.9,
                 flavor="leidenalg")

    # ---- DP-lineage subcluster (cell 48): split DC/DP from DP ----
    # The notebook hard-coded cluster '7' as the DP lineage. That index is
    # run-specific -- in YOUR run, identify the DP-lineage cluster by its markers
    # (Sox2, Sox18, Cd200, Col11a1, Alx4 high) and set DP_CLUSTER accordingly.
    DP_CLUSTER = None  # e.g. "7"
    if DP_CLUSTER is not None:
        sc.tl.leiden(adata, restrict_to=("leiden_scVI_0.9", [DP_CLUSTER]),
                     resolution=0.1, key_added="sub_cluster", flavor="leidenalg")
    else:
        adata.obs["sub_cluster"] = adata.obs["leiden_scVI_0.9"]
        print(">> DP_CLUSTER not set: skipping DC/DP vs DP split. Inspect markers, "
              "set DP_CLUSTER, and re-run the last block to match the paper.")

    # ---- marker check: confirm per-cluster markers match the paper ----
    sc.tl.rank_genes_groups(adata, "leiden_scVI_0.9", method="wilcoxon",
                            key_added="rank_leiden09")
    adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X  # for the dotplot below

    adata.write(OUT_H5AD)
    print(f">> written {OUT_H5AD}")
    print(">> next: open in Seurat (h5ad -> RDS via sceasy/zellkonverter); in Seurat v5 "
          "the 'counts' layer becomes the counts layer and X(lognorm) the data layer of "
          "the RNA assay. Annotate leiden_scVI_0.9 against PAPER_MARKERS using PAPER_LABELS.")


# ----------------------------------------------------------------------------
# Reference constants from the marker-study notebook (for your manual annotation)
# ----------------------------------------------------------------------------
PAPER_LABELS = [  # the 34 final cell types (desired_order in the notebook)
    "APM", "BK-1", "BK-2", "BK-3", "Chond-Fib", "Coch-Fib", "DC/DP", "DP", "DS",
    "FIB-1", "FIB-2", "FIB-3", "FIB-4", "HFKC-1", "HFKC-2", "HP/HG", "LEC",
    "Lymphocyte", "Mac-1", "Mac-2", "Mast", "Mast-low", "Melanocyte", "Merkel",
    "Muscle", "Neutrophil", "Pericyte", "Pre-adipo", "SBK", "Schwann", "VEC",
    "emFIB-1", "emFIB-2", "emK",
]
PAPER_MARKERS = [  # 64-gene dotplot panel used to assign the labels above
    "Crabp1", "Twist2", "Pdgfra", "Dpp4", "Wif1", "Apcdd1", "Dlk1", "Agtr2",
    "Fabp4", "Mfap5", "Ebf2", "Meox2", "Igfbp7", "Col11a1", "Cd200", "Acan",
    "Sox18", "Cxcr4", "Sox2", "Alpl", "Alx4", "Col23a1", "Actg2", "Itga8",
    "Coch", "Matn4", "Dlx5", "Trp63", "Lhx2", "Edar", "Barx2", "Sox9",
    "Il11ra1", "Krt79", "Sostdc1", "Apoe", "Msx2", "Krt71", "Krt5", "Krt14",
    "Krt1", "Krt10", "Cpa3", "Mcpt4", "Cxcr2", "Itgam", "Ptprc", "Cd3g",
    "Cd68", "Cd86", "Mrc1", "Cd163", "Krt8", "Krt18", "Msc", "Ttn", "Pax7",
    "Rgs5", "Sox10", "Dct", "Tyr", "Pecam1", "Cdh5", "Lyve1",
]


if __name__ == "__main__":
    main()