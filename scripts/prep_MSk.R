#!/usr/bin/env Rscript
# =============================================================================
# prep_MSk.R
# Purpose : Load a reference .h5 (AnnData) file and save as Seurat v5 RDS.
#           Optionally build a panel gene list from proseg gene-metadata files
#           and restrict the reference to those genes.
#
# Usage   :
#   # Build panel from proseg dirs + restrict reference
#   Rscript prep_MSk.R -i ref.h5 -o ./reference/reference_prepared.rds \
#                       -d ./proseg_results -m TRUE
#
#   # Use previously generated panel
#   Rscript prep_MSk.R -i ref.h5 -o ./reference/reference_prepared.rds \
#                       -p ./reference/x5_gene_panel.txt
#
#   # Full transcriptome (no panel)
#   Rscript prep_MSk.R -i ref.h5 -o ./reference/reference_prepared.rds
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdf5r)
  library(Matrix)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "Path to input .h5 file [required]"),
  make_option(c("-o", "--output"),
              type = "character", default = "reference_prepared.rds",
              help = "Path for output .rds file [default: reference_prepared.rds]"),
  make_option(c("-d", "--panel-dir"),
              type = "character", default = "./proseg_results",
              help = "Root proseg directory to scan for gene-metadata.csv.gz [default: ./proseg_results]"),
  make_option(c("-m", "--make-panel"),
              action = "store_true", default = FALSE,
              help = "Build panel from proseg gene-metadata files [default: FALSE]"),
  make_option(c("-p", "--panel"),
              type = "character", default = NULL,
              help = "Path to existing panel gene list (plain text or CSV with 'Name' column)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input))      stop("--input is required.", call. = FALSE)
if (!file.exists(opt$input)) stop(sprintf("File not found: %s", opt$input), call. = FALSE)

if (dir.exists(opt$output)) {
  opt$output <- file.path(opt$output, "reference_prepared.rds")
}

message("\n========================================")
message("  prep_MSk.R  |  Seurat v5 reference prep")
message("========================================\n")
message(sprintf("Input      : %s", opt$input))
message(sprintf("Output     : %s", opt$output))

# ===========================================================================
# Panel gene list — build from proseg or load from file
# ===========================================================================
panel_genes <- NULL

if (opt[["make-panel"]]) {
  
  message(sprintf("Panel dir  : %s", opt[["panel-dir"]]))
  
  csv_hits <- list.files(
    opt[["panel-dir"]],
    pattern    = "^gene-metadata\\.csv\\.gz$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(csv_hits) == 0)
    stop("No gene-metadata.csv.gz found under: ", opt$panel_dir, call. = FALSE)
  
  message(sprintf("\n[PANEL] Scanning %d gene-metadata file(s):", length(csv_hits)))
  
  all_genes <- character(0)
  for (f in csv_hits) {
    region <- basename(dirname(f))
    dt     <- fread(f, select = "gene")
    genes  <- unique(dt$gene)
    message(sprintf("  %-20s %d genes", region, length(genes)))
    all_genes <- union(all_genes, genes)
  }
  
  panel_genes <- sort(all_genes)
  message(sprintf("[PANEL] Total unique panel genes: %d", length(panel_genes)))
  
  # Write to fixed location
  panel_out <- "./reference/x5_gene_panel.txt"
  dir.create(dirname(panel_out), showWarnings = FALSE, recursive = TRUE)
  writeLines(panel_genes, panel_out)
  message(sprintf("[PANEL] Saved to: %s\n", panel_out))
  
} else if (!is.null(opt$panel)) {
  
  if (!file.exists(opt$panel))
    stop(sprintf("Panel file not found: %s", opt$panel), call. = FALSE)
  
  first_line <- readLines(opt$panel, n = 1)
  
  if (grepl(",", first_line)) {
    panel_df <- read.csv(opt$panel, stringsAsFactors = FALSE)
    name_col <- grep("^Name$", colnames(panel_df), ignore.case = TRUE, value = TRUE)
    if (length(name_col) == 0) {
      message("  [PANEL] No 'Name' column in CSV; using first column")
      panel_genes <- unique(trimws(panel_df[[1]]))
    } else {
      panel_genes <- unique(trimws(panel_df[[name_col[1]]]))
    }
  } else {
    panel_genes <- unique(trimws(readLines(opt$panel)))
  }
  
  panel_genes <- panel_genes[nzchar(panel_genes)]
  message(sprintf("[PANEL] Loaded %d genes from: %s\n", length(panel_genes), opt$panel))
  
} else {
  message("Panel  : none (full transcriptome)\n")
}

# ===========================================================================
# Read .h5
# ===========================================================================
message("Opening .h5 file...")
h5   <- hdf5r::H5File$new(opt$input, mode = "r")
keys <- h5$ls()$name
message(sprintf("Root groups: %s", paste(keys, collapse = ", ")))

message("Reading var (genes) and obs (cells) identifiers...")

var_grp    <- h5[["var"]]
var_index  <- if (var_grp$attr_exists("_index")) var_grp$attr_open("_index")$read() else "_index"
gene_names <- var_grp[[var_index]][]
ngenes     <- length(gene_names)

obs_grp    <- h5[["obs"]]
obs_index  <- if (obs_grp$attr_exists("_index")) obs_grp$attr_open("_index")$read() else "_index"
cell_names <- obs_grp[[obs_index]][]
ncells     <- length(cell_names)

message(sprintf("  var (genes): %d | obs (cells): %d", ngenes, ncells))

message("Reading count matrix (dense)...")
X <- h5[["layers/raw"]]

if (inherits(X, "H5Group")) {
  indices <- X[["indices"]][] + 1L
  indptr  <- X[["indptr"]][]
  data_v  <- as.numeric(X[["data"]][])
  
  if (X$attr_exists("shape")) {
    shape <- X$attr_open("shape")$read()
  } else if ("shape" %in% X$ls()$name) {
    shape <- X[["shape"]][]
  } else {
    shape <- c(length(indptr) - 1L, max(indices))
  }
  
  encoding <- "csr_matrix"
  if (X$attr_exists("encoding-type")) encoding <- X$attr_open("encoding-type")$read()
  else if (X$attr_exists("encoding_type")) encoding <- X$attr_open("encoding_type")$read()
  
  if (grepl("csc", encoding, ignore.case = TRUE)) {
    counts <- Matrix::sparseMatrix(i = indices, p = indptr, x = data_v,
                                   dims = c(shape[1], shape[2]), repr = "C")
  } else {
    csr <- Matrix::sparseMatrix(j = indices, p = indptr, x = data_v,
                                dims = c(shape[1], shape[2]), repr = "C")
    counts <- Matrix::t(csr)
  }
  if (nrow(counts) == ncells && ncol(counts) == ngenes) {
    counts <- Matrix::t(counts)
  }
} else {
  mat <- X$read()
  message(sprintf("  Raw dense dimensions: %d x %d", nrow(mat), ncol(mat)))
  
  if (nrow(mat) == ngenes && ncol(mat) == ncells) {
    message("  Orientation: genes x cells — no transpose needed")
    counts <- Matrix::Matrix(mat, sparse = TRUE)
  } else if (nrow(mat) == ncells && ncol(mat) == ngenes) {
    message("  Orientation: cells x genes — transposing")
    counts <- Matrix::Matrix(t(mat), sparse = TRUE)
  } else {
    stop(sprintf(
      "Matrix dims [%d, %d] don't match var (%d) x obs (%d) in either orientation.",
      nrow(mat), ncol(mat), ngenes, ncells
    ), call. = FALSE)
  }
}

rownames(counts) <- gene_names
colnames(counts) <- cell_names
message(sprintf("  Final matrix: %d genes x %d cells", nrow(counts), ncol(counts)))

# ===========================================================================
# Panel restriction
# ===========================================================================
use_panel <- FALSE

if (!is.null(panel_genes)) {
  shared  <- intersect(rownames(counts), panel_genes)
  missing <- setdiff(panel_genes, rownames(counts))
  
  message(sprintf("  [PANEL] Matched: %d / %d panel genes",
                  length(shared), length(panel_genes)))
  if (length(missing) > 0) {
    message(sprintf("  [PANEL] %d absent from reference (first 20): %s",
                    length(missing), paste(head(missing, 20), collapse = ", ")))
  }
  if (length(shared) < 100)
    stop(sprintf("[PANEL] Only %d shared genes — check gene name format.",
                 length(shared)), call. = FALSE)
  
  counts <- counts[shared, , drop = FALSE]
  message(sprintf("  [PANEL] Restricted: %d genes x %d cells",
                  nrow(counts), ncol(counts)))
  use_panel <- TRUE
}

# ===========================================================================
# Cell metadata
# ===========================================================================
message("Reading cell metadata...")

col_order_raw <- obs_grp$attr_open("column-order")$read()
if (length(col_order_raw) == 1 && grepl(",", col_order_raw)) {
  obs_cols <- trimws(strsplit(col_order_raw, ",")[[1]])
} else {
  obs_cols <- as.character(col_order_raw)
}
message(sprintf("  Columns to read: %d", length(obs_cols)))

cat_grp <- NULL
cat_names <- character(0)
if ("__categories" %in% obs_grp$ls()$name) {
  cat_grp   <- obs_grp[["__categories"]]
  cat_names <- cat_grp$ls()$name
}
message(sprintf("  Columns with __categories entries: %d", length(cat_names)))

meta_list <- vector("list", length(obs_cols))
names(meta_list) <- obs_cols
failed_cols <- character(0)

for (col in obs_cols) {
  result <- tryCatch({
    obs_child_names <- obs_grp$ls()$name
    if (!(col %in% obs_child_names)) {
      message(sprintf("    WARNING: '%s' not found in obs, skipping", col))
      return(NULL)
    }
    
    ds        <- obs_grp[[col]]
    dtype_str <- ds$get_type()$to_text()
    
    if (col %in% cat_names) {
      codes <- as.integer(ds$read()) + 1L
      levs  <- cat_grp[[col]]$read()
      codes[codes < 1L] <- NA_integer_
      levs_chr <- as.character(levs)
      factor(levs_chr[codes], levels = levs_chr)
    } else if (grepl("H5T_ENUM", dtype_str)) {
      raw_vals <- ds$read()
      as.logical(as.integer(raw_vals))
    } else {
      ds$read()
    }
  }, error = function(e) {
    message(sprintf("    ERROR reading '%s': %s", col, conditionMessage(e)))
    failed_cols <<- c(failed_cols, col)
    rep(NA, ncells)
  })
  
  meta_list[[col]] <- result
}

meta_list <- meta_list[!vapply(meta_list, is.null, logical(1))]
meta_df   <- as.data.frame(meta_list, row.names = cell_names,
                           check.names = FALSE, stringsAsFactors = FALSE)

message(sprintf("  Metadata: %d columns x %d cells", ncol(meta_df), nrow(meta_df)))
if (length(failed_cols) > 0) {
  message(sprintf("  Failed columns (set to NA): %s", paste(failed_cols, collapse = ", ")))
}

# ===========================================================================
# Embeddings (PCA, UMAP) from obsm
# ===========================================================================
message("Reading embeddings from obsm...")

embed_list <- list()
if ("obsm" %in% keys) {
  obsm_grp   <- h5[["obsm"]]
  obsm_names <- obsm_grp$ls()$name
  
  for (emb_name in obsm_names) {
    emb_ds <- obsm_grp[[emb_name]]
    if (!inherits(emb_ds, "H5Group")) {
      emb_mat <- emb_ds$read()
      if (ncol(emb_mat) != ncells && nrow(emb_mat) == ncells) {
        embed_list[[emb_name]] <- emb_mat
      } else {
        embed_list[[emb_name]] <- t(emb_mat)
      }
      message(sprintf("  %s: %d cells x %d dims",
                      emb_name, nrow(embed_list[[emb_name]]), ncol(embed_list[[emb_name]])))
    }
  }
}

h5$close_all()

# ===========================================================================
# Build Seurat object
# ===========================================================================
message("Creating Seurat object...")
obj <- CreateSeuratObject(counts    = counts,
                          meta.data = meta_df,
                          project   = "reference",
                          min.cells = 0, min.features = 0)

if (!use_panel) {
  
  if ("X_pca" %in% names(embed_list)) {
    pca_mat <- embed_list[["X_pca"]]
    rownames(pca_mat) <- cell_names
    colnames(pca_mat) <- paste0("PC_", seq_len(ncol(pca_mat)))
    obj[["pca"]] <- CreateDimReducObject(embeddings = pca_mat, key = "PC_", assay = "RNA")
    message(sprintf("  Added PCA: %d dims", ncol(pca_mat)))
  }
  
  if ("X_umap" %in% names(embed_list)) {
    umap_mat <- embed_list[["X_umap"]]
    rownames(umap_mat) <- cell_names
    colnames(umap_mat) <- paste0("umap_", seq_len(ncol(umap_mat)))
    obj[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "umap_", assay = "RNA")
    message(sprintf("  Added UMAP: %d dims", ncol(umap_mat)))
  }
  
} else {
  message("  [PANEL] Skipping pre-existing PCA/UMAP (computed on full transcriptome)")
}

# ===========================================================================
# Join layers, save
# ===========================================================================
message("Joining layers (Seurat v5)...")
obj <- JoinLayers(obj)

message(sprintf("Saving to: %s", opt$output))
saveRDS(obj, file = opt$output)

message(sprintf("\nDone — %d cells | %d features", ncol(obj), nrow(obj)))
message(sprintf("Metadata columns (%d):", ncol(obj@meta.data)))
print(colnames(obj@meta.data))
message(sprintf("Reductions: %s", paste(Reductions(obj), collapse = ", ")))

if (use_panel) {
  message(sprintf("\n[PANEL] Reference restricted to %d panel genes.", nrow(obj)))
  message("[PANEL] PCA/UMAP will be recomputed by the annotation pipeline.")
}