################################################################################
# helper_functions.R
#
# Proseg + Xenium spatial analysis helpers for Seurat v5.
################################################################################

library(data.table)
library(Seurat)
library(SeuratObject)
library(anndataR)
library(rhdf5)
library(arrow)
library(dplyr)
library(Matrix)
library(scuttle)
library(spacexr)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(SingleR)
library(ggplot2)
library(sf)
library(patchwork)
library(viridis)
library(BiocParallel)
library(ggpmisc)
library(future)
library(BPCells)
library(cowplot)
library(gridExtra)
library(scales)

options(future.globals.maxSize = 200 * 1024^3)
Assays <- SeuratObject::Assays

#######################################################################
#######################################################################
#                                                                     #
#                         SETUP AND UTILITIES                         #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: timestamps, logging, default palette, coordinate column resolution

#TIMESTAMP PREFIX FOR CONSOLE OUTPUT
.ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

#MESSAGE LOGGER
log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), paste0(...), "\n")
}

#FALLBACK PALETTE WHEN NO CONTRAST PALETTE IS AVAILABLE
default_palette <- c(
  "#2166AC", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD", "#17BECF",
  "#E377C2", "#8C564B", "#BCBD22", "#1F77B4", "#AEC7E8", "#FFBB78",
  "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2", "#DBDB8D",
  "#9EDAE5", "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#5254A3", "#6B6ECF", "#9C9EDE", "#E7BA52", "#BD9E39", "#AD494A",
  "#D6616B", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6", "#E6550D",
  "#FD8D3C", "#31A354", "#74C476", "#756BB1", "#FDAE6B", "#A1D99B",
  "#DADAEB", "#636363", "#969696", "#525252", "#FDD0A2", "#C7E9C0"
)

#RESOLVE A COORDINATE COLUMN FROM ALTERNATIVE NAMES
resolve_coord_col <- function(meta_colnames, preferred, alternatives) {
  if (preferred %in% meta_colnames) return(preferred)
  for (alt in alternatives) {
    if (alt %in% meta_colnames) return(alt)
  }
  return(preferred)  # fall through, will be caught later
}

#######################################################################
#######################################################################
#                                                                     #
#                          COHORT DEFINITION                          #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: sample sheet, regions list, region checks, centroid accessor

######################################
# SAMPLE SHEET                       #
######################################

#read the cohort definition: sample_id,condition,proseg_dir
read_sample_sheet <- function(path = "samples.csv", verbose = TRUE) {

  if (!file.exists(path)) stop("sample sheet not found: ", path)
  sheet <- data.table::fread(path, data.table = FALSE,
                             colClasses = "character")

  required <- c("sample_id", "condition", "proseg_dir")
  missing  <- setdiff(required, colnames(sheet))
  if (length(missing))
    stop("sample sheet missing column(s): ", paste(missing, collapse = ", "))
  if (anyDuplicated(sheet$sample_id))
    stop("duplicate sample_id in ", path)
  if (any(grepl("[^A-Za-z0-9._-]", sheet$sample_id)))
    stop("sample_id must contain only [A-Za-z0-9._-]")

  rownames(sheet) <- sheet$sample_id
  if (verbose)
    cat(.ts(), " sample sheet:", nrow(sheet), "samples |",
        paste(names(table(sheet$condition)), table(sheet$condition),
              sep = "=", collapse = " "), "\n")
  sheet
}

#build the regions list from the sample sheet
regions_from_sheet <- function(sheet,
                               rds_name = "{sample_id}_proseg_seurat.rds",
                               require_exists = TRUE) {

  regions <- lapply(seq_len(nrow(sheet)), function(i) {
    sample_id <- sheet$sample_id[i]
    list(sample_id = sample_id,
         condition = sheet$condition[i],
         dir       = sheet$proseg_dir[i],
         rds_path  = file.path(sheet$proseg_dir[i],
                               gsub("{sample_id}", sample_id, rds_name,
                                    fixed = TRUE)))
  })
  names(regions) <- sheet$sample_id

  if (require_exists) {
    missing <- names(regions)[!vapply(regions, function(region)
      file.exists(region$rds_path), logical(1))]
    if (length(missing))
      stop("object not built for: ", paste(missing, collapse = ", "))
  }
  regions
}

#check every region object before any expensive work
check_proseg_regions <- function(regions, verbose = TRUE) {

  report <- lapply(names(regions), function(sample_id) {
    obj    <- readRDS(regions[[sample_id]]$rds_path)
    issues <- character(0)

    if (!"RNA" %in% Assays(obj)) issues <- c(issues, "no RNA assay")
    if (anyDuplicated(colnames(obj))) issues <- c(issues, "duplicate cell names")
    if (!length(Images(obj))) issues <- c(issues, "no FOV")
    if (length(Images(obj))) {
      coords <- proseg_coords(obj)
      if (is.null(coords)) issues <- c(issues, "no centroids boundary")
      else if (nrow(coords) != ncol(obj))
        issues <- c(issues, "centroids do not cover all cells")
    }

    result <- data.frame(sample_id = sample_id, n_cells = ncol(obj),
                         n_genes = nrow(obj), n_fov = length(Images(obj)),
                         issues = paste(issues, collapse = " | "),
                         stringsAsFactors = FALSE)
    rm(obj); gc(verbose = FALSE)
    result
  })
  report <- do.call(rbind, report)

  if (verbose) {
    print(report)
    failed <- report[nzchar(report$issues), ]
    if (nrow(failed)) warning("regions with issues: ",
                              paste(failed$sample_id, collapse = ", "))
  }
  invisible(report)
}

#per-cell centroids, one row per cell, in microns
proseg_coords <- function(obj, fov = NULL, cells = NULL) {

  if (!length(Images(obj))) return(NULL)
  fov <- fov %||% Images(obj)[1]

  coords <- tryCatch(
    GetTissueCoordinates(obj, image = fov, which = "centroids"),
    error = function(e) NULL)
  if (is.null(coords) || !nrow(coords)) return(NULL)

  ids <- if ("cell" %in% colnames(coords)) as.character(coords$cell) else
    rownames(coords)

  #segmentation boundaries give one row per vertex
  if (anyDuplicated(ids)) {
    coords <- data.table::as.data.table(
      list(cell = ids, x = coords$x, y = coords$y))[
        , .(x = mean(x), y = mean(y)), by = cell]
    ids    <- coords$cell
    coords <- data.frame(x = coords$x, y = coords$y)
  } else {
    coords <- data.frame(x = coords$x, y = coords$y)
  }

  rownames(coords) <- ids
  if (!is.null(cells)) coords <- coords[intersect(cells, ids), , drop = FALSE]
  coords
}

#######################################################################
#######################################################################
#                                                                     #
#                      PROSEG OBJECT CONSTRUCTION                     #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: h5ad loading, FOV attachment, per-region object building

######################################
# LOAD PROSEG H5AD                   #
######################################

load_proseg_h5ad <- function(h5ad_path) {
  
  if (!file.exists(h5ad_path)) stop("h5ad not found: ", h5ad_path)
  
  message("Reading h5ad: ", h5ad_path)
  adata <- read_h5ad(h5ad_path)
  seu <- adata$as_Seurat()
  
  # Fix layer naming: anndataR stores as "X", Seurat expects "counts"
  available_layers <- Layers(seu[["RNA"]])
  if ("X" %in% available_layers && !"counts" %in% available_layers) {
    seu[["RNA"]]$counts <- seu[["RNA"]]$X
    seu[["RNA"]]$X <- NULL
    message("Renamed layer 'X' -> 'counts'")
  }
  
  # Store spatial coords in metadata if available
  spatial_coords <- adata$obsm[["spatial"]]
  if (!is.null(spatial_coords)) {
    seu$centroid_x <- spatial_coords[, 1]
    seu$centroid_y <- spatial_coords[, 2]
  }
  
  message("Seurat object: ", ncol(seu), " cells, ", nrow(seu), " genes")
  return(seu)
}

######################################
# ATTACH FOVS                        #
######################################

#PROSEG SEGMENTATION AND CENTROIDS AS A NATIVE FOV
attach_proseg_fov <- function(seu,
                              vertices_csv_path,
                              centroids_csv_path = NULL,
                              fov_name = "proseg") {
  require(SeuratObject)
  require(Seurat)
  require(data.table)
  
  if (!file.exists(vertices_csv_path)) {
    stop("Vertex CSV not found: ", vertices_csv_path)
  }
  
  # ---- Load vertices ----
  message("Reading polygon vertices: ", vertices_csv_path)
  vtx <- fread(vertices_csv_path)
  
  if (!all(c("x", "y", "cell") %in% colnames(vtx))) {
    stop("Vertex CSV must have columns: x, y, cell")
  }
  
  # Filter to cells in Seurat object
  seurat_cells <- colnames(seu)
  vtx <- vtx[vtx$cell %in% seurat_cells, ]
  n_cells_with_seg <- length(unique(vtx$cell))
  message("  ", nrow(vtx), " vertices for ", n_cells_with_seg, " cells")
  
  # Build Segmentation
  seg <- CreateSegmentation(as.data.frame(vtx))
  
  # Build Centroids
  if (!is.null(centroids_csv_path) && file.exists(centroids_csv_path)) {
    message("Reading centroids: ", centroids_csv_path)
    cents_df <- fread(centroids_csv_path)
    cents_df <- as.data.frame(cents_df[cents_df$cell %in% seurat_cells, ])
  } else if (all(c("centroid_x", "centroid_y") %in% colnames(seu@meta.data))) {
    cents_df <- data.frame(
      x    = seu$centroid_x,
      y    = seu$centroid_y,
      cell = seurat_cells
    )
  } else {
    stop("No centroids available. Provide centroids_csv_path or ensure ",
         "centroid_x/centroid_y are in metadata (from load_proseg_h5ad).")
  }
  
  cents <- CreateCentroids(cents_df)
  
  # Build FOV
  fov <- CreateFOV(
    coords = list("segmentation" = seg, "centroids" = cents),
    type   = c("segmentation", "centroids"),
    assay  = DefaultAssay(seu),
    key    = paste0(fov_name, "_")
  )
  
  # Subset to shared cells
  shared <- intersect(Cells(fov[["segmentation"]]), seurat_cells)
  fov <- subset(fov, cells = shared)
  
  seu[[fov_name]] <- fov
  message("Attached FOV '", fov_name, "' with segmentation + centroids (",
          length(shared), " cells).")
  
  return(seu)
}

#XENIUM CELL BOUNDARIES AS A SECOND FOV
attach_xenium_fov <- function(seu,
                              xenium_dir,
                              fov_name = "xenium",
                              boundary_type = "cell") {
  require(arrow)
  require(dplyr)
  require(SeuratObject)
  require(Seurat)
  
  parquet_file <- switch(boundary_type,
                         "cell"    = file.path(xenium_dir, "cell_boundaries.parquet"),
                         "nucleus" = file.path(xenium_dir, "nucleus_boundaries.parquet"),
                         stop("boundary_type must be 'cell' or 'nucleus'")
  )
  if (!file.exists(parquet_file)) stop("Not found: ", parquet_file)
  
  message("Reading Xenium ", boundary_type, " boundaries...")
  bdf <- read_parquet(parquet_file)
  
  required <- c("cell_id", "vertex_x", "vertex_y")
  if (!all(required %in% colnames(bdf))) {
    stop("Parquet must contain: ", paste(required, collapse = ", "))
  }
  
  if ("label_id" %in% colnames(bdf)) {
    bdf <- bdf %>% arrange(cell_id, label_id)
  }
  
  # ---- Map Xenium cell_id -> Seurat colnames ----
  seurat_names <- colnames(seu)
  xenium_ids <- as.character(bdf$cell_id)
  unique_xenium <- unique(xenium_ids)
  
  # Try direct match
  # Try direct match
  direct <- sum(unique_xenium %in% seurat_names)
  
  if (direct > 0) {
    message("  Direct match: ", direct, " / ", length(unique_xenium), " Xenium cells")
    id_map <- setNames(unique_xenium, unique_xenium)
  } else {
    # Map via original_cell_id in metadata
    meta     <- seu@meta.data
    ocid_col <- intersect(c("original_cell_id", "cell_id"), colnames(meta))
    
    if (length(ocid_col) == 0) {
      message("  No 'original_cell_id' / 'cell_id' column in metadata. ",
              "Cannot map Xenium boundaries to proseg cells. ",
              "Skipping Xenium FOV attachment.")
      return(seu)
    }
    
    id_map   <- setNames(seurat_names, as.character(meta[[ocid_col[1]]]))
    n_mapped <- sum(!is.na(id_map[unique_xenium]))
    message("  Mapped via ", ocid_col[1], ": ", n_mapped,
            " / ", length(unique_xenium), " Xenium cells")
    
    if (n_mapped == 0) {
      message("  Zero Xenium IDs match proseg metadata. ",
              "This is expected when proseg ran with --use-xenium-mask: ",
              "proseg's 'original_cell_id' is then an integer mask label, ",
              "not a Xenium alphanumeric cell_id. ",
              "Skipping Xenium FOV attachment.")
      return(seu)
    }
  }
  
  # Apply mapping to all vertex rows
  bdf$cell <- as.character(id_map[xenium_ids])
  bdf <- bdf[!is.na(bdf$cell), ]
  
  # Build vertex data.frame for CreateSegmentation
  vtx_df <- data.frame(
    x    = bdf$vertex_x,
    y    = bdf$vertex_y,
    cell = bdf$cell
  )
  
  n_cells <- length(unique(vtx_df$cell))
  message("  Building segmentation from ", nrow(vtx_df),
          " vertices (", n_cells, " cells)")
  
  seg <- CreateSegmentation(vtx_df)
  
  # Build Centroids
  cells_parquet <- file.path(xenium_dir, "cells.parquet")
  if (file.exists(cells_parquet)) {
    cmeta <- read_parquet(cells_parquet)
    if (all(c("cell_id", "x_centroid", "y_centroid") %in% colnames(cmeta))) {
      cmeta$cell <- as.character(id_map[as.character(cmeta$cell_id)])
      cmeta <- cmeta[!is.na(cmeta$cell) & cmeta$cell %in% seurat_names, ]
      cents_df <- data.frame(
        x    = cmeta$x_centroid,
        y    = cmeta$y_centroid,
        cell = cmeta$cell
      )
    } else {
      cents_df <- NULL
    }
  } else {
    cents_df <- NULL
  }
  
  # Fallback: compute centroids from vertices
  if (is.null(cents_df) || nrow(cents_df) == 0) {
    message("  Computing centroids from vertex means...")
    cents_df <- vtx_df %>%
      group_by(cell) %>%
      summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
      as.data.frame()
  }
  
  cents <- CreateCentroids(cents_df)
  
  # Build FOV
  fov <- CreateFOV(
    coords = list("segmentation" = seg, "centroids" = cents),
    type   = c("segmentation", "centroids"),
    assay  = DefaultAssay(seu),
    key    = paste0(fov_name, "_")
  )
  
  shared <- intersect(Cells(fov[["segmentation"]]), seurat_names)
  fov <- subset(fov, cells = shared)
  
  seu[[fov_name]] <- fov
  message("Attached FOV '", fov_name, "' with segmentation + centroids (",
          length(shared), " cells).")
  
  return(seu)
}

#XENIUM NUCLEUS BOUNDARIES AS A THIRD FOV
attach_xenium_nuclei <- function(seu,
                                 xenium_dir,
                                 fov_name = "xenium",
                                 boundary_name = "nuclei") {
  require(arrow)
  require(dplyr)
  require(SeuratObject)
  
  parquet_file <- file.path(xenium_dir, "nucleus_boundaries.parquet")
  if (!file.exists(parquet_file)) stop("Not found: ", parquet_file)
  if (!fov_name %in% Images(seu)) {
    message("FOV '", fov_name, "' not found (Xenium boundaries were likely ",
            "skipped upstream). Skipping nuclei attachment.")
    return(seu)
  }
  
  message("Reading Xenium nucleus boundaries...")
  bdf <- read_parquet(parquet_file)
  
  # Same cell name mapping as attach_xenium_fov
  seurat_names <- colnames(seu)
  xenium_ids <- as.character(bdf$cell_id)
  unique_xenium <- unique(xenium_ids)
  
  direct <- sum(unique_xenium %in% seurat_names)
  if (direct > 0) {
    id_map <- setNames(unique_xenium, unique_xenium)
  } else {
    meta <- seu@meta.data
    ocid_col <- intersect(c("original_cell_id", "cell_id"), colnames(meta))
    if (length(ocid_col) > 0) {
      id_map <- setNames(seurat_names, as.character(meta[[ocid_col[1]]]))
    } else {
      stop("No cell name mapping found.")
    }
  }
  
  bdf$cell <- as.character(id_map[as.character(bdf$cell_id)])
  bdf <- bdf[!is.na(bdf$cell) & bdf$cell %in% seurat_names, ]
  
  vtx_df <- data.frame(x = bdf$vertex_x, y = bdf$vertex_y, cell = bdf$cell)
  
  nuclei_seg <- CreateSegmentation(vtx_df)
  seu[[fov_name]][[boundary_name]] <- nuclei_seg
  
  n_cells <- length(unique(vtx_df$cell))
  message("Added '", boundary_name, "' boundary to FOV '", fov_name,
          "' (", n_cells, " cells).")
  return(seu)
}

######################################
# BUILD REGION OBJECTS               #
######################################

#COUNTS PLUS ALL AVAILABLE FOVS FOR ONE REGION
load_proseg_full <- function(h5ad_path,
                             vertices_csv_path,
                             centroids_csv_path = NULL,
                             xenium_dir = NULL,
                             proseg_fov_name = "proseg",
                             xenium_fov_name = "xenium") {
  
  # Load counts + metadata
  seu <- load_proseg_h5ad(h5ad_path)
  
  # Attach proseg segmentation FOV
  seu <- attach_proseg_fov(
    seu, vertices_csv_path,
    centroids_csv_path = centroids_csv_path,
    fov_name = proseg_fov_name
  )
  
  # Optionally attach Xenium boundaries
  if (!is.null(xenium_dir)) {
    seu <- attach_xenium_fov(seu, xenium_dir, fov_name = xenium_fov_name)
    seu <- attach_xenium_nuclei(seu, xenium_dir, fov_name = xenium_fov_name)
  }
  
  return(seu)
}

#CONVERT EVERY PROSEG ZARR OUTPUT INTO A SEURAT OBJECT
build_proseg_seurat <- function(proseg_dir,
                                xenium_dir,
                                resolution = 0.7,
                                npcs       = 30,
                                umap_dims  = 1:20,
                                overwrite  = FALSE) {
  
  # Discover all zarr outputs
  zarr_hits <- list.files(
    proseg_dir,
    pattern      = "^proseg-output\\.zarr$",
    recursive    = TRUE,
    full.names   = TRUE,
    include.dirs = TRUE
  )
  
  # zarr sits directly in <proseg_dir>/<sample>/proseg-output.zarr
  sample_dirs <- unique(dirname(zarr_hits))
  
  if (length(sample_dirs) == 0) {
    stop("No proseg-output.zarr found under: ", proseg_dir)
  }
  
  cat("Found", length(sample_dirs), "sample(s)\n")
  
  for (sample_dir in sample_dirs) {
    
    sample_name <- basename(sample_dir)
    rds_path    <- file.path(sample_dir, paste0(sample_name, "_proseg_seurat.rds"))
    
    cat("\n========================================\n")
    cat("Sample:", sample_name, "\n")
    cat("========================================\n")
    
    if (!overwrite && file.exists(rds_path)) {
      cat("  Skipping — RDS already exists:", rds_path, "\n")
      next
    }
    
    zarr_path      <- file.path(sample_dir, "proseg-output.zarr")
    xenium_sample  <- file.path(xenium_dir, sample_name)
    
    # Validate 
    if (!dir.exists(xenium_sample)) {
      warning(sample_name, ": xenium dir not found at ", xenium_sample, " — skipping.")
      next
    }
    
    # Python: zarr -> h5ad + vertex CSV
    cat("  Converting zarr to h5ad...\n")
    zconv$convert(
      zarr_path  = zarr_path,
      output_dir = sample_dir
    )
    
    h5ad_path      <- file.path(sample_dir, "proseg-anndata.h5ad")
    vertices_path  <- file.path(sample_dir, "proseg-seg-vertices.csv.gz")
    centroids_path <- file.path(sample_dir, "proseg-centroids.csv.gz")
    
    missing <- c(h5ad_path, vertices_path, centroids_path)
    missing <- missing[!file.exists(missing)]
    if (length(missing) > 0) {
      warning(sample_name, ": conversion produced missing files — skipping.\n",
              paste("  ", missing, collapse = "\n"))
      next
    }
    
    # Build Seurat with dual FOVs
    cat("  Building Seurat object...\n")
    obj <- load_proseg_full(
      h5ad_path          = h5ad_path,
      vertices_csv_path  = vertices_path,
      centroids_csv_path = centroids_path,
      #xenium_dir         = xenium_sample,
      proseg_fov_name    = "proseg",
      #xenium_fov_name    = "xenium"
    )
    
    obj <- subset(obj, subset = nCount_RNA > 1)  
    
    cat("  Cells after subset:", ncol(obj), "\n")
    
    for (fov_nm in Images(obj)) {
      fov_cells <- Cells(obj[[fov_nm]])
      cat("  FOV ", fov_nm, " cells: ", length(fov_cells), "\n", sep = "")
      
      cent_cells <- Cells(obj[[fov_nm]], boundary = "centroids")
      cat("  ", fov_nm, " FOV cells in Seurat: ",
          sum(cent_cells %in% colnames(obj)), " / ",
          length(cent_cells), "\n", sep = "")
    }
    
    saveRDS(obj, file = rds_path)
    cat("  Saved:", rds_path, "\n")
    
    rm(obj)
    gc()
  }
}

#######################################################################
#######################################################################
#                                                                     #
#                           QUALITY CONTROL                           #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: assay cleaning, count and feature filters, doublet removal, h5 loading

#DROP DERIVED LAYERS AND REDUCTIONS FROM AN ASSAY
ensure_clean_assay <- function(obj, nm) {
  
  # Force default assay to RNA
  if (DefaultAssay(obj) != "RNA") {
    log_msg(nm, " | Switching default assay from '",
            DefaultAssay(obj), "' to 'RNA'")
    DefaultAssay(obj) <- "RNA"
  }
  
  available <- Layers(obj[["RNA"]])
  
  # If counts layer is missing, try to recover
  if (!"counts" %in% available) {
    
    # Some v3/v4 -> v5 conversions store raw counts under "data"
    if ("data" %in% available) {
      log_msg(nm, " | No 'counts' layer found; copying 'data' -> 'counts'")
      obj[["RNA"]]$counts <- LayerData(obj[["RNA"]], layer = "data")
    } else {
      # Check for split layers (counts.X, counts.Y, ...)
      split_counts <- grep("^counts\\.", available, value = TRUE)
      if (length(split_counts) > 0) {
        log_msg(nm, " | Found split count layers (",
                paste(split_counts, collapse = ", "),
                "); joining...")
        obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
      } else {
        stop("[", nm, "] RNA assay has no 'counts' layer and no recoverable ",
             "alternative. Available layers: ",
             paste(available, collapse = ", "))
      }
    }
  }
  
  log_msg(nm, " | RNA layers: ",
          paste(Layers(obj[["RNA"]]), collapse = ", "))
  
  return(obj)
}

#FILTER CELLS ON COUNTS, FEATURES AND MITOCHONDRIAL FRACTION
qc_filter <- function(obj, config, nm) {
  
  n_before <- ncol(obj)
  
  # Compute percent.mt if not already present
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  }
  
  # Apply filters
  keep <- obj$nFeature_RNA >= config$min_genes &
    obj$nFeature_RNA <= config$max_genes &
    obj$percent.mt   <= config$max_mt_pct
  
  n_remove <- sum(!keep)
  
  if (n_remove > 0) {
    obj <- subset(obj, cells = colnames(obj)[keep])
  }
  
  log_msg(nm, " | QC: ", n_before, " -> ", ncol(obj), " cells ",
          "(removed ", n_remove, "; genes [", config$min_genes, "-",
          config$max_genes, "], MT <= ", config$max_mt_pct, "%)")
  
  return(obj)
}

#REMOVE DOUBLETS PER BATCH
remove_doublets <- function(obj, batch_col, nm) {
  
  n_before <- ncol(obj)
  
  # Convert to SCE for scDblFinder
  sce <- as.SingleCellExperiment(obj)
  
  # Run scDblFinder; if batch_col has only 1 level, don't pass samples
  batches <- obj@meta.data[[batch_col]]
  n_batches <- length(unique(batches))
  
  if (n_batches > 1) {
    sce <- scDblFinder(sce, samples = batch_col)
  } else {
    sce <- scDblFinder(sce)
  }
  
  # Transfer calls back to Seurat
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  
  n_doublets <- sum(obj$scDblFinder.class == "doublet")
  
  # Remove doublets
  obj <- subset(obj, scDblFinder.class == "singlet")
  
  log_msg(nm, " | Doublets: ", n_doublets, " / ", n_before,
          " (", round(100 * n_doublets / n_before, 1), "%) removed -> ",
          ncol(obj), " cells")
  
  return(obj)
}

#LOAD A CELLRANGER H5 INTO A SEURAT OBJECT
load_h5_to_seurat <- function(h5_path, sample_name) {
  
  counts <- Read10X_h5(h5_path)
  
  # Read10X_h5 returns a list when multiple modalities exist (e.g. GEX + ADT)
  # Take Gene Expression if so
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  
  obj <- CreateSeuratObject(
    counts       = counts,
    project      = sample_name,
    min.cells    = 0,
    min.features = 0
  )
  
  obj$sample_name <- sample_name
  
  return(obj)
}

#PARSE A SAMPLE NAME FROM AN H5 FILENAME
sample_name_from_h5 <- function(filename) {
  
  name <- sub("\\.h5$", "", filename)
  
  stripped <- sub("_?filtered_feature_bc_matrix$", "", name)
  stripped <- sub("_?raw_feature_bc_matrix$", "", stripped)
  
  if (nchar(stripped) == 0) {
    return(name)
  }
  
  return(stripped)
}

#######################################################################
#######################################################################
#                                                                     #
#                       SKETCH BASED ANNOTATION                       #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: merge regions, sketch, cluster, markers, manual labels, projection

######################################
# MERGE PROSEG REGIONS               #
######################################

#merge all regions into one object, keeping segmentation FOVs
merge_proseg_regions <- function(regions,
                                 sample_sheet = NULL,
                                 min_counts   = 10,
                                 on_disk      = FALSE,
                                 on_disk_dir  = "./bpcells_counts",
                                 sample_col   = "sample_id",
                                 verbose      = TRUE) {

  object_list <- list()

  for (sample_id in names(regions)) {

    if (verbose) cat(.ts(), " reading", sample_id, "\n")
    obj <- readRDS(regions[[sample_id]]$rds_path)
    DefaultAssay(obj) <- "RNA"

    #drop low-count cells before merging
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
    keep   <- colnames(obj)[Matrix::colSums(counts) >= min_counts]
    if (length(keep) < ncol(obj)) obj <- subset(obj, cells = keep)
    rm(counts)

    #prefix cell names so they are unique across the cohort
    obj$proseg_cell_id <- colnames(obj)
    new_names <- paste0(sample_id, "_", colnames(obj))
    
    #detach FOVs, rename object and boundaries separately, reattach as one
    #FOV named after the region. RenameCells on the object does not
    #propagate into boundary cell names.
    fov_list <- Images(obj)
    boundaries <- list()
    for (fov in fov_list) {
      boundaries[[fov]] <- obj[[fov]]
      obj[[fov]] <- NULL
    }
    
    obj <- RenameCells(obj, new.names = new_names)
    
    if (length(boundaries)) {
      fov <- boundaries[[1]]
      fov <- RenameCells(fov, new.names = paste0(sample_id, "_", Cells(fov)))
      DefaultAssay(fov) <- "RNA"
      obj[[sample_id]] <- fov
    }

    obj[[sample_col]] <- sample_id
    condition <- regions[[sample_id]]$condition
    if (is.null(condition) && !is.null(sample_sheet))
      condition <- sample_sheet[sample_id, "condition"]
    obj$condition <- if (is.null(condition)) NA_character_ else condition

    if (verbose) cat(.ts(), "  ", ncol(obj), "cells |", length(Images(obj)),
                     "FOV\n")
    object_list[[sample_id]] <- obj
  }

  if (verbose) cat(.ts(), " merging\n")
  obj <- if (length(object_list) == 1) object_list[[1]] else
    merge(x = object_list[[1]], y = object_list[-1])
  rm(object_list); gc(verbose = FALSE)

  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

  if (on_disk) {
    if (verbose) cat(.ts(), " writing counts to", on_disk_dir, "\n")
    dir.create(dirname(on_disk_dir), recursive = TRUE, showWarnings = FALSE)
    BPCells::write_matrix_dir(
      mat = LayerData(obj, assay = "RNA", layer = "counts"),
      dir = on_disk_dir, overwrite = TRUE)
    LayerData(obj, assay = "RNA", layer = "counts") <-
      BPCells::open_matrix_dir(dir = on_disk_dir)
    gc(verbose = FALSE)
  }

  #split by region so sketching draws a fixed number of cells per region
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[sample_col]])

  if (verbose) {
    cat(.ts(), " merged:", ncol(obj), "cells,", nrow(obj), "genes\n")
    cat(.ts(), " FOVs:", paste(Images(obj), collapse = ", "), "\n")
  }
  obj
}


######################################
# SKETCH                             #
######################################

#hvg selection with a dispersion fallback for fractional proseg counts
find_variable_features_safe <- function(obj, assay = NULL, n_features = 2000) {

  assay   <- assay %||% DefaultAssay(obj)
  warnings_seen <- NULL

  obj <- withCallingHandlers(
    FindVariableFeatures(obj, assay = assay, nfeatures = n_features,
                         verbose = FALSE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  degenerate <- FALSE
  tryCatch({
    hvf <- HVFInfo(obj, assay = assay, method = "vst")
    if ("variance.standardized" %in% colnames(hvf))
      degenerate <- sum(is.nan(hvf$variance.standardized)) > 0.5 * nrow(hvf)
  }, error = function(e)
    degenerate <<- any(grepl("NaN", warnings_seen, fixed = TRUE)))

  if (degenerate) {
    cat(.ts(), " vst degenerate, using dispersion\n")
    obj <- FindVariableFeatures(obj, assay = assay, nfeatures = n_features,
                                selection.method = "dispersion",
                                verbose = FALSE)
  }
  obj
}

#add a sketch assay holding a fixed number of cells per region
sketch_proseg_object <- function(obj,
                                 sketch_cells   = 20000,
                                 method         = c("leverage", "uniform"),
                                 n_features     = 2000,
                                 score_features = 4000,
                                 sketch_assay   = "sketch",
                                 seed           = 123,
                                 verbose        = TRUE) {

  method <- match.arg(method)
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- find_variable_features_safe(obj, "RNA", n_features)

  #LeverageScore aborts above 5000 features, so cap the scoring features
  features <- head(VariableFeatures(obj), score_features)
  if (verbose) cat(.ts(), " sketching", sketch_cells, "cells per region on",
                   length(features), "features\n")

  obj <- SketchData(
    object         = obj,
    ncells         = sketch_cells,
    method         = if (method == "leverage") "LeverageScore" else "Uniform",
    sketched.assay = sketch_assay,
    features       = features,
    seed           = seed,
    verbose        = FALSE)

  DefaultAssay(obj) <- sketch_assay
  obj[[sketch_assay]] <- JoinLayers(obj[[sketch_assay]])

  if (verbose) cat(.ts(), " sketch assay:", ncol(obj[[sketch_assay]]),
                   "cells\n")
  obj
}


######################################
# CLUSTER SKETCH                     #
######################################

#cluster the sketch assay and record cluster labels
cluster_sketch <- function(obj,
                           dims         = 1:30,
                           resolutions  = c(0.3, 0.5, 0.8, 1.2),
                           resolution   = 0.5,
                           n_features   = 2000,
                           integrate    = FALSE,
                           sample_col   = "sample_id",
                           sketch_assay = "sketch",
                           cluster_col  = "sketch_cluster",
                           verbose      = TRUE) {

  DefaultAssay(obj) <- sketch_assay

  obj <- NormalizeData(obj, assay = sketch_assay, verbose = FALSE)
  obj <- find_variable_features_safe(obj, sketch_assay, n_features)

  counts   <- LayerData(obj, assay = sketch_assay, layer = "counts")
  detected <- rownames(counts)[Matrix::rowSums(counts) > 0]
  features <- intersect(VariableFeatures(obj, assay = sketch_assay), detected)
  rm(counts)

  if (verbose) cat(.ts(), " pca on", length(features), "features\n")
  obj <- ScaleData(obj, assay = sketch_assay, features = features,
                   verbose = FALSE)
  obj <- RunPCA(obj, assay = sketch_assay, features = features,
                npcs = max(50, max(dims)), verbose = FALSE)

  reduction <- "pca"
  if (integrate) {
    if (!requireNamespace("harmony", quietly = TRUE))
      stop("integrate = TRUE requires harmony")
    if (verbose) cat(.ts(), " harmony on", sample_col, "\n")
    obj <- harmony::RunHarmony(obj, group.by.vars = sample_col,
                               reduction.use = "pca",
                               reduction.save = "harmony",
                               assay.use = sketch_assay, verbose = FALSE)
    reduction <- "harmony"
  }

  obj <- FindNeighbors(obj, reduction = reduction, dims = dims,
                       verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolutions, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction, dims = dims,
                 reduction.name = "umap", return.model = TRUE,
                 verbose = FALSE)

  obj <- set_sketch_resolution(obj, resolution, sketch_assay, cluster_col,
                               verbose = verbose)

  #report whether clusters are region specific
  composition <- prop.table(
    table(obj@meta.data[[cluster_col]], obj@meta.data[[sample_col]]), 1)
  region_specific <- sum(apply(composition, 1, max) > 0.8)
  if (verbose)
    cat(.ts(), " clusters >80% from one region:", region_specific, "of",
        nrow(composition), "\n")

  obj
}

#switch to a precomputed clustering resolution
set_sketch_resolution <- function(obj, resolution, sketch_assay = "sketch",
                                  cluster_col = "sketch_cluster",
                                  verbose = TRUE) {

  column <- paste0(sketch_assay, "_snn_res.", resolution)
  if (!column %in% colnames(obj@meta.data)) {
    obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
    column <- grep(paste0("_snn_res\\.", resolution, "$"),
                   colnames(obj@meta.data), value = TRUE)[1]
  }
  obj@meta.data[[cluster_col]] <- as.character(obj@meta.data[[column]])
  Idents(obj) <- cluster_col
  if (verbose)
    cat(.ts(), " resolution", resolution, "->",
        length(unique(obj@meta.data[[cluster_col]])), "clusters\n")
  obj
}


######################################
# CLUSTER MARKERS                    #
######################################

#markers per sketch cluster
sketch_cluster_markers <- function(obj,
                                   cluster_col  = "sketch_cluster",
                                   sketch_assay = "sketch",
                                   output_dir   = "./annotated_data",
                                   only_pos     = TRUE,
                                   min_pct      = 0.1,
                                   logfc        = 0.25,
                                   format       = "csv",
                                   cache        = TRUE) {

  dir  <- file.path(output_dir, "_sketch")
  file <- paste0("markers_", cluster_col)
  path <- file.path(dir, paste0(file, ".csv"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  if (cache && file.exists(path)) {
    cat(.ts(), " cached markers:", path, "\n")
    return(data.table::fread(path, data.table = FALSE))
  }

  DefaultAssay(obj) <- sketch_assay
  Idents(obj) <- cluster_col
  markers <- FindAllMarkers(obj, assay = sketch_assay, group.by = cluster_col,
                            only.pos = only_pos,
                            min.pct = min_pct, logfc.threshold = logfc,
                            verbose = FALSE)

  written <- .write_markers(markers, dir, file, format = format)
  cat(.ts(), " markers:", paste(written, collapse = ", "), "\n")
  markers
}

#top n markers per cluster as a named list
top_cluster_markers <- function(markers, n = 15, order_by = "avg_log2FC") {
  markers <- markers[order(markers$cluster, -markers[[order_by]]), ]
  lapply(split(markers$gene, markers$cluster), head, n)
}


######################################
# MANUAL ANNOTATION                  #
######################################

#write a csv with one row per cluster for manual labelling
write_cluster_template <- function(obj,
                                   markers     = NULL,
                                   cluster_col = "sketch_cluster",
                                   n_genes     = 15,
                                   output_dir  = "./annotated_data",
                                   path        = NULL) {

  path <- path %||% file.path(output_dir, "_sketch", "cluster_annotation.csv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  clusters <- sort(unique(as.character(obj@meta.data[[cluster_col]])))
  counts   <- table(as.character(obj@meta.data[[cluster_col]]))
  top      <- if (!is.null(markers)) top_cluster_markers(markers, n_genes)

  template <- data.frame(
    cluster   = clusters,
    n_cells   = as.integer(counts[clusters]),
    top_genes = vapply(clusters, function(cluster)
      if (is.null(top[[cluster]])) "" else
        paste(top[[cluster]], collapse = ", "), character(1)),
    cell_type = "",
    stringsAsFactors = FALSE)

  data.table::fwrite(template, path)
  cat(.ts(), " fill the cell_type column of:", path, "\n")
  invisible(template)
}

#read a completed template as a named vector
read_cluster_template <- function(output_dir = "./annotated_data",
                                  path       = NULL,
                                  label_col  = "cell_type") {

  path <- path %||% file.path(output_dir, "_sketch", "cluster_annotation.csv")
  template <- data.table::fread(path, data.table = FALSE)
  if (any(is.na(template[[label_col]]) | template[[label_col]] == ""))
    stop("unfilled rows in ", path)
  setNames(as.character(template[[label_col]]), as.character(template$cluster))
}

#apply cluster labels to the sketch cells
label_sketch_clusters <- function(obj,
                                  labels      = NULL,
                                  cluster_col = "sketch_cluster",
                                  label_col   = "cell_type",
                                  output_dir  = "./annotated_data",
                                  template    = NULL,
                                  verbose     = TRUE) {

  if (is.null(labels))
    labels <- read_cluster_template(output_dir, template)
  if (is.list(labels)) labels <- unlist(labels)

  clusters <- as.character(obj@meta.data[[cluster_col]])
  missing  <- setdiff(unique(na.omit(clusters)), names(labels))
  if (length(missing))
    stop("no label for cluster(s): ", paste(missing, collapse = ", "))

  obj@meta.data[[label_col]] <- unname(labels[clusters])
  Idents(obj) <- label_col

  if (verbose) print(sort(table(obj@meta.data[[label_col]]), decreasing = TRUE))
  obj
}

zoom_region <- function(obj, fov, x, y, features = NULL,
                        group_col = "cell_type", colours = NULL,
                        size = 1.5, flip_xy = TRUE) {
  
  obj[["zoom_tmp"]] <- Crop(obj[[fov]], x = x, y = y, coords = "tissue")
  n <- length(Cells(obj[["zoom_tmp"]]))
  cat("cells in window:", n, "\n")
  if (n == 0) stop("empty window, check the bbox")
  
  plot <- if (is.null(features)) {
    ImageDimPlot(obj, fov = "zoom_tmp", group.by = group_col,
                 cols = colours, size = size, border.color = NA,
                 border.size = 0, flip_xy = flip_xy)
  } else {
    ImageFeaturePlot(obj, fov = "zoom_tmp", features = features,
                     size = size, border.color = NA, border.size = 0)
  }
  plot
}

######################################
# PROJECT LABELS TO ALL CELLS        #
######################################

#extend sketch labels and embeddings to every cell in the full assay
project_sketch_labels <- function(obj,
                                  label_col    = "cell_type",
                                  assay        = "RNA",
                                  sketch_assay = "sketch",
                                  dims         = 1:30,
                                  k_weight     = 50,
                                  verbose      = TRUE) {

  if (!label_col %in% colnames(obj@meta.data))
    stop(label_col, " not found, run label_sketch_clusters() first")

  #Seurat's CreateCategoryMatrix() rewrites underscores to dashes in
  #identity names, and TransferLablesNN() reads predicted labels straight
  #off those column names, so "Mural_Muscle" comes back "Mural-Muscle" and
  #splits into two compartments. keep the originals to map them back.
  originals <- unique(na.omit(as.character(obj@meta.data[[label_col]])))
  mangled   <- gsub("_", "-", originals)
  if (anyDuplicated(mangled))
    warning("labels differing only in _ vs - cannot be told apart after ",
            "projection: ",
            paste(originals[duplicated(mangled) | duplicated(mangled,
                  fromLast = TRUE)], collapse = ", "))

  DefaultAssay(obj) <- sketch_assay
  refdata <- setNames(list(label_col), paste0(label_col, "_projected"))

  if (verbose) cat(.ts(), " projecting to", ncol(obj[[assay]]), "cells\n")
  obj <- ProjectData(
    object             = obj,
    assay              = assay,
    full.reduction     = "pca.full",
    sketched.assay     = sketch_assay,
    sketched.reduction = "pca",
    umap.model         = "umap",
    dims               = dims,
    refdata            = refdata,
    k.weight           = k_weight,
    verbose            = FALSE)

  #single label column covering every cell
  projected <- paste0(label_col, "_projected")

  #undo the underscore to dash rewrite described above
  lookup  <- setNames(originals, mangled)
  guessed <- as.character(obj@meta.data[[projected]])
  hit     <- !is.na(guessed) & guessed %in% names(lookup)
  guessed[hit] <- unname(lookup[guessed[hit]])
  obj@meta.data[[projected]] <- guessed

  full_labels <- obj@meta.data[[projected]]
  sketch_labels <- as.character(obj@meta.data[[label_col]])
  full_labels[!is.na(sketch_labels)] <- sketch_labels[!is.na(sketch_labels)]
  obj@meta.data[[label_col]] <- full_labels

  DefaultAssay(obj) <- assay
  Idents(obj) <- label_col

  if (verbose) {
    cat(.ts(), " labelled:", sum(!is.na(obj@meta.data[[label_col]])), "of",
        ncol(obj), "cells\n")
    print(sort(table(obj@meta.data[[label_col]]), decreasing = TRUE))
  }
  obj
}

#######################################################################
#######################################################################
#                                                                     #
#                        TWO LEVEL ANNOTATION                         #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: subsetting a broad compartment, reclustering it, and folding the
#fine labels back into the parent object. Labels can be supplied either as
#named vectors in the calling script or as filled csv templates on disk.

######################################
# PATHS                              #
######################################

subcluster_dir <- function(output_dir, broad_type) {
  file.path(output_dir, "_subclusters",
            gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", broad_type)))
}


######################################
# SUBSET AND STRIP                   #
######################################

#clean subset with everything derived from the parent clustering removed.
#the parent pca.full in particular has to go, because ProjectData() skips
#computation whenever full.reduction already exists in the object
prepare_subcluster <- function(obj,
                               cells,
                               sample_col     = "sample_id",
                               keep_images    = FALSE,
                               min_per_region = 500,
                               verbose        = TRUE) {

  sub <- subset(obj, cells = cells)

  if (!keep_images) for (nm in Images(sub))     sub[[nm]] <- NULL
  for (nm in setdiff(Assays(sub), "RNA"))       sub[[nm]] <- NULL
  for (nm in Reductions(sub))                   sub[[nm]] <- NULL
  for (nm in Graphs(sub))                       sub[[nm]] <- NULL
  for (nm in Neighbors(sub))                    sub[[nm]] <- NULL
  slot(sub, "tools")$TransferSketchLabels <- NULL

  drop <- grep(paste0("^sketch_snn_res\\.|^RNA_snn_res\\.|^seurat_clusters$|",
                      "^sketch_cluster$|^subcluster$|^sub_type$|",
                      "_projected$|_projected\\.score$"),
               colnames(sub@meta.data), value = TRUE)
  for (nm in drop) sub@meta.data[[nm]] <- NULL

  DefaultAssay(sub) <- "RNA"
  sub[["RNA"]] <- JoinLayers(sub[["RNA"]])

  #only split by region if every region contributes enough cells to
  #normalise and pick hvgs on
  per_region <- table(sub@meta.data[[sample_col]])
  per_region <- per_region[per_region > 0]

  if (length(per_region) > 1 && min(per_region) >= min_per_region) {
    sub[["RNA"]] <- split(
      sub[["RNA"]],
      f = droplevels(factor(sub@meta.data[[sample_col]])))
    if (verbose) cat(.ts(), "   split into", length(per_region), "layers\n")
  } else if (verbose) {
    cat(.ts(), "   single RNA layer\n")
  }

  if (verbose) cat(.ts(), "   subset:", ncol(sub), "cells\n")
  sub
}


######################################
# CLUSTER A SUBSET                   #
######################################

#small subsets cluster directly on RNA, large ones go through the same
#sketch -> cluster -> project route as the parent. either way cluster_col
#ends up covering every cell in the subset
subcluster_broad <- function(sub,
                             dims         = 1:20,
                             resolutions  = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
                             resolution   = 0.6,
                             n_features   = 2000,
                             sketch_cells = 20000,
                             direct_max   = 75000,
                             cluster_col  = "subcluster",
                             verbose      = TRUE) {

  DefaultAssay(sub) <- "RNA"

  if (ncol(sub) <= direct_max) {

    if (verbose) cat(.ts(), "   clustering", ncol(sub), "cells directly\n")

    sub <- NormalizeData(sub, verbose = FALSE)
    sub <- find_variable_features_safe(sub, "RNA", n_features)

    layer_sums <- sapply(
      Layers(sub[["RNA"]], search = "counts"),
      function(l) Matrix::rowSums(LayerData(sub, assay = "RNA", layer = l)))
    detected <- rownames(sub)[rowSums(as.matrix(layer_sums)) > 0]
    features <- intersect(VariableFeatures(sub, assay = "RNA"), detected)

    npcs <- min(50, length(features) - 1, ncol(sub) - 1)
    dims <- 1:min(max(dims), npcs)

    sub <- ScaleData(sub, features = features, verbose = FALSE)
    sub <- RunPCA(sub, features = features, npcs = npcs, verbose = FALSE)
    sub <- FindNeighbors(sub, reduction = "pca", dims = dims, verbose = FALSE)
    sub <- FindClusters(sub, resolution = resolutions, verbose = FALSE)
    sub <- RunUMAP(sub, reduction = "pca", dims = dims,
                   reduction.name = "umap", verbose = FALSE)

    column <- grep(paste0("_snn_res\\.", resolution, "$"),
                   colnames(sub@meta.data), value = TRUE)[1]
    sub@meta.data[[cluster_col]] <- as.character(sub@meta.data[[column]])

  } else {

    if (verbose) cat(.ts(), "   sketching", ncol(sub), "cells\n")

    sub <- sketch_proseg_object(sub,
                                sketch_cells = sketch_cells,
                                method       = "leverage",
                                n_features   = n_features,
                                verbose      = verbose)

    sub <- cluster_sketch(sub,
                          dims        = dims,
                          resolutions = resolutions,
                          resolution  = resolution,
                          n_features  = n_features,
                          integrate   = FALSE,
                          cluster_col = cluster_col,
                          verbose     = verbose)

    k_weight <- max(5, min(50, floor(ncol(sub[["sketch"]]) / 20)))
    sub <- project_sketch_labels(sub,
                                 label_col = cluster_col,
                                 dims      = dims,
                                 k_weight  = k_weight,
                                 verbose   = verbose)
  }

  DefaultAssay(sub) <- "RNA"
  Idents(sub) <- cluster_col

  if (verbose)
    cat(.ts(), "  ", length(unique(sub@meta.data[[cluster_col]])),
        "subclusters\n")
  sub
}


######################################
# MARKERS FOR A SUBSET               #
######################################

#writes a marker table as csv, xlsx, or both. xlsx gets one sheet per
#cluster plus an "all" sheet, which is the format worth opening when you
#are eyeballing markers cluster by cluster
.write_markers <- function(markers,
                           dir,
                           file,
                           format    = "csv",
                           extra     = NULL,
                           split_col = "cluster") {

  format <- match.arg(format, c("csv", "xlsx"), several.ok = TRUE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character()

  if ("csv" %in% format) {
    path <- file.path(dir, paste0(file, ".csv"))
    data.table::fwrite(markers, path)
    paths <- c(paths, path)
    if (!is.null(extra))
      for (nm in names(extra)) {
        side <- file.path(dir, paste0(nm, ".csv"))
        data.table::fwrite(extra[[nm]], side)
        paths <- c(paths, side)
      }
  }

  if ("xlsx" %in% format) {

    if (!requireNamespace("openxlsx", quietly = TRUE))
      stop("format \"xlsx\" needs the openxlsx package")

    sheets <- list()
    if (!is.null(extra)) sheets <- extra
    sheets[["all"]] <- markers

    if (nrow(markers) && split_col %in% colnames(markers)) {
      keys <- as.character(markers[[split_col]])
      by   <- split(markers, factor(keys, levels = unique(keys)))
      names(by) <- paste0("c", names(by))
      sheets <- c(sheets, by)
    }

    #excel caps sheet names at 31 characters and bans a few symbols
    names(sheets) <- make.unique(substr(
      gsub("[\\[\\]:*?/\\\\]", "_", names(sheets)), 1, 28))

    path <- file.path(dir, paste0(file, ".xlsx"))
    openxlsx::write.xlsx(sheets, path, overwrite = TRUE)
    paths <- c(paths, path)
  }

  invisible(paths)
}


#dir = NULL computes markers without touching disk, which is how you test a
#candidate resolution before committing to it. format may be "csv", "xlsx",
#or c("csv", "xlsx")
subcluster_markers <- function(sub,
                               cluster_col = "subcluster",
                               dir         = NULL,
                               only_pos    = TRUE,
                               logfc       = 0.25,
                               min_pct     = 0.2,
                               format      = "csv",
                               cache       = FALSE) {

  if (!cluster_col %in% colnames(sub@meta.data))
    stop(cluster_col, " not found. available groupings: ",
         paste(subcluster_resolutions(sub)$column, collapse = ", "))

  path <- if (!is.null(dir)) file.path(dir, "markers.csv") else NULL

  if (cache && !is.null(path) && file.exists(path))
    return(data.table::fread(path, data.table = FALSE))

  assay <- if ("sketch" %in% Assays(sub)) "sketch" else "RNA"

  if (length(Layers(sub[[assay]], search = "data")) > 1)
    sub[[assay]] <- JoinLayers(sub[[assay]])

  DefaultAssay(sub) <- assay
  Idents(sub) <- cluster_col

  markers <- FindAllMarkers(sub,
                            assay           = assay,
                            group.by        = cluster_col,
                            only.pos        = only_pos,
                            logfc.threshold = logfc,
                            min.pct         = min_pct,
                            verbose         = FALSE)

  if (!is.null(dir))
    .write_markers(markers, dir, "markers", format = format)

  markers
}


######################################
# COMPARE TWO CLUSTERS DIRECTLY      #
######################################

#subcluster_markers() is one vs rest, so two similar clusters will both
#return the same genes against the pooled remainder. this tests them
#against each other only. a near empty result means the two clusters are
#not separable and should probably be merged
compare_clusters <- function(sub,
                             ident_1,
                             ident_2,
                             cluster_col = "subcluster",
                             assay       = NULL,
                             logfc       = 0.2,
                             min_pct     = 0.1,
                             n           = 30,
                             dir         = NULL,
                             format      = "csv",
                             verbose     = TRUE) {

  if (!cluster_col %in% colnames(sub@meta.data))
    stop(cluster_col, " not found. available groupings: ",
         paste(subcluster_resolutions(sub)$column, collapse = ", "))

  groups <- as.character(sub@meta.data[[cluster_col]])
  ident_1 <- as.character(ident_1)
  ident_2 <- as.character(ident_2)

  for (id in c(ident_1, ident_2))
    if (!id %in% groups)
      stop("cluster ", id, " not present in ", cluster_col)

  assay <- assay %||% if ("sketch" %in% Assays(sub)) "sketch" else "RNA"
  if (length(Layers(sub[[assay]], search = "data")) > 1)
    sub[[assay]] <- JoinLayers(sub[[assay]])
  DefaultAssay(sub) <- assay

  res <- FindMarkers(sub,
                     assay           = assay,
                     group.by        = cluster_col,
                     ident.1         = ident_1,
                     ident.2         = ident_2,
                     only.pos        = FALSE,
                     logfc.threshold = logfc,
                     min.pct         = min_pct,
                     verbose         = FALSE)

  if (nrow(res)) {
    res$gene      <- rownames(res)
    res$higher_in <- ifelse(res$avg_log2FC > 0, ident_1, ident_2)
    res <- res[order(-abs(res$avg_log2FC)),
               c("gene", "avg_log2FC", "pct.1", "pct.2",
                 "p_val_adj", "higher_in")]
    rownames(res) <- NULL
  }

  if (!is.null(dir) && nrow(res))
    .write_markers(res, dir,
                   paste0("compare_", ident_1, "_vs_", ident_2),
                   format = format, split_col = "higher_in")

  if (verbose) {
    cat(.ts(), " ", ident_1, " (n=", sum(groups == ident_1), ") vs ",
        ident_2, " (n=", sum(groups == ident_2), "): ",
        nrow(res), " genes at logfc>", logfc, "\n", sep = "")
    if (nrow(res)) print(head(res, n))
    else cat(.ts(), " nothing separates them, consider merging\n")
  }
  invisible(res)
}


######################################
# NEAREST NEIGHBOUR MARKERS          #
######################################

#one vs rest drowns out the differences between clusters that sit next to
#each other, because the pooled remainder is dominated by distant cell
#types. this locates each cluster's k closest neighbours by centroid
#distance in the reduction and tests against those only, so the genes
#returned are the ones that actually separate look-alikes
neighbour_markers <- function(sub,
                              cluster_col = "subcluster",
                              reduction   = "pca",
                              dims        = NULL,
                              k           = 1,
                              assay       = NULL,
                              logfc       = 0.2,
                              min_pct     = 0.1,
                              only_pos    = TRUE,
                              dir         = NULL,
                              format      = "csv",
                              verbose     = TRUE) {

  if (!cluster_col %in% colnames(sub@meta.data))
    stop(cluster_col, " not found. available groupings: ",
         paste(subcluster_resolutions(sub)$column, collapse = ", "))
  if (!reduction %in% Reductions(sub))
    stop("no reduction called '", reduction, "'. available: ",
         paste(Reductions(sub), collapse = ", "))

  assay <- assay %||% if ("sketch" %in% Assays(sub)) "sketch" else "RNA"
  if (length(Layers(sub[[assay]], search = "data")) > 1)
    sub[[assay]] <- JoinLayers(sub[[assay]])
  DefaultAssay(sub) <- assay

  #the sketch reduction only covers sketch cells, which is also all the
  #de can use, so the intersection is the right cell set either way
  emb   <- Embeddings(sub[[reduction]])
  cells <- intersect(rownames(emb), colnames(sub[[assay]]))
  emb   <- emb[cells, , drop = FALSE]

  if (is.null(dims)) dims <- 1:min(ncol(emb), 20)
  emb <- emb[, dims, drop = FALSE]

  groups <- as.character(sub@meta.data[cells, cluster_col])
  ok     <- !is.na(groups)
  emb    <- emb[ok, , drop = FALSE]
  groups <- groups[ok]

  centroids <- t(sapply(split(seq_along(groups), groups),
                        function(i) colMeans(emb[i, , drop = FALSE])))
  ids <- rownames(centroids)
  if (length(ids) < 2) stop("need at least two clusters")

  distances <- as.matrix(stats::dist(centroids))
  diag(distances) <- Inf
  k <- min(k, length(ids) - 1)

  tables <- list()
  map    <- list()

  for (id in ids) {

    nb  <- names(sort(distances[id, ]))[1:k]
    tag <- paste(nb, collapse = "+")

    res <- tryCatch(
      FindMarkers(sub,
                  assay           = assay,
                  group.by        = cluster_col,
                  ident.1         = id,
                  ident.2         = nb,
                  only.pos        = only_pos,
                  logfc.threshold = logfc,
                  min.pct         = min_pct,
                  verbose         = FALSE),
      error = function(e) NULL)

    n_genes <- if (is.null(res)) NA_integer_ else nrow(res)

    if (!is.null(res) && nrow(res)) {
      res$gene    <- rownames(res)
      res$cluster <- id
      res$versus  <- tag
      res <- res[order(-res$avg_log2FC), ]
      rownames(res) <- NULL
      tables[[id]] <- res
    }

    map[[id]] <- data.frame(
      cluster  = id,
      n_cells  = sum(groups == id),
      nearest  = tag,
      distance = round(min(distances[id, ]), 2),
      n_genes  = n_genes,
      stringsAsFactors = FALSE)
  }

  pairs <- do.call(rbind, map)
  rownames(pairs) <- NULL

  markers <- if (length(tables)) do.call(rbind, tables) else data.frame()
  if (nrow(markers))
    markers <- markers[, c("gene", "cluster", "versus", "avg_log2FC",
                           "pct.1", "pct.2", "p_val", "p_val_adj")]
  rownames(markers) <- NULL
  attr(markers, "pairs") <- pairs

  if (!is.null(dir))
    .write_markers(markers, dir, "markers_neighbour",
                   format = format,
                   extra  = list(neighbour_pairs = pairs))

  if (verbose) {
    cat(.ts(), " each cluster vs its", k, "nearest in", reduction, "\n")
    print(pairs)
    weak <- pairs$cluster[!is.na(pairs$n_genes) & pairs$n_genes == 0]
    if (length(weak))
      cat(.ts(), " nothing separates", paste(weak, collapse = ", "),
          "from their neighbour\n")
  }
  markers
}


######################################
# PASS 1: RECLUSTER EVERY COMPARTMENT#
######################################

#loops the broad labels, reclusters each one, writes markers, an empty
#annotation template and the subcluster object to
#<output_dir>/_subclusters/<compartment>/
run_subclustering <- function(obj,
                              broad_col    = "broad_type",
                              output_dir   = "./annotated_data",
                              skip_types   = NULL,
                              min_cells    = 500,
                              dims         = 1:20,
                              resolutions  = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
                              resolution   = 0.6,
                              n_features   = 2000,
                              sketch_cells = 20000,
                              direct_max   = 75000,
                              cluster_col  = "subcluster",
                              n_genes      = 15,
                              markers      = FALSE,
                              overwrite    = FALSE,
                              verbose      = TRUE) {

  if (!broad_col %in% colnames(obj@meta.data))
    stop(broad_col, " not found, run label_sketch_clusters() first")

  types <- setdiff(
    sort(unique(na.omit(as.character(obj@meta.data[[broad_col]])))),
    skip_types)

  rows <- list()

  for (type in types) {

    cells <- colnames(obj)[which(obj@meta.data[[broad_col]] == type)]
    dir   <- subcluster_dir(output_dir, type)
    path  <- file.path(dir, "subcluster.rds")

    if (!overwrite && file.exists(path)) {
      if (verbose) cat(.ts(), " ", type, ": cached\n", sep = "")
      sub <- readRDS(path)
      rows[[type]] <- data.frame(
        broad_type    = type,
        n_cells       = ncol(sub),
        n_subclusters = length(unique(sub@meta.data[[cluster_col]])),
        status        = "cached",
        dir           = dir,
        stringsAsFactors = FALSE)
      rm(sub); gc(verbose = FALSE)
      next
    }

    if (length(cells) < min_cells) {
      if (verbose) cat(.ts(), " ", type, ": only ", length(cells),
                       " cells, skipped\n", sep = "")
      rows[[type]] <- data.frame(
        broad_type = type, n_cells = length(cells), n_subclusters = NA_integer_,
        status = "too small", dir = NA_character_, stringsAsFactors = FALSE)
      next
    }

    if (verbose) cat("\n", .ts(), " ==== ", type, " ====\n", sep = "")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    sub <- prepare_subcluster(obj, cells = cells, verbose = verbose)
    sub <- subcluster_broad(sub,
                            dims         = dims,
                            resolutions  = resolutions,
                            resolution   = resolution,
                            n_features   = n_features,
                            sketch_cells = sketch_cells,
                            direct_max   = direct_max,
                            cluster_col  = cluster_col,
                            verbose      = verbose)

    #markers are deferred by default, they depend on the resolution you
    #settle on after looking at the umap. see finalise_subcluster()
    if (isTRUE(markers)) {
      marker_table <- subcluster_markers(sub, cluster_col = cluster_col,
                                         dir = dir)
      write_cluster_template(sub, marker_table, cluster_col = cluster_col,
                             n_genes = n_genes,
                             path = file.path(dir, "annotation.csv"))
      rm(marker_table)
    }

    saveRDS(sub, path)

    rows[[type]] <- data.frame(
      broad_type    = type,
      n_cells       = ncol(sub),
      n_subclusters = length(unique(sub@meta.data[[cluster_col]])),
      status        = "done",
      dir           = dir,
      stringsAsFactors = FALSE)

    rm(sub); gc(verbose = FALSE)
  }

  summary_table <- do.call(rbind, rows)
  rownames(summary_table) <- NULL
  if (verbose) print(summary_table)
  invisible(summary_table)
}


######################################
# INSPECT ONE COMPARTMENT            #
######################################

#every resolution stored on an object, with the cluster count for each
subcluster_resolutions <- function(object) {

  cols <- grep("_snn_res\\.", colnames(object@meta.data), value = TRUE)
  if (!length(cols)) return(data.frame())

  out <- data.frame(
    column     = cols,
    resolution = as.numeric(gsub("^.*_snn_res\\.", "", cols)),
    n_clusters = as.integer(sapply(cols, function(x)
      length(unique(na.omit(object@meta.data[[x]]))))),
    stringsAsFactors = FALSE)

  out <- out[order(out$resolution), ]
  rownames(out) <- NULL
  out
}


#pull a subcluster back into the session for labelling in RStudio.
#returns list(obj = , markers = , top = , dir = , broad_type = ,
#             resolutions = )
load_subcluster <- function(broad_type,
                            output_dir  = "./annotated_data",
                            cluster_col = "subcluster",
                            n_genes     = 15,
                            plot        = FALSE) {

  dir  <- subcluster_dir(output_dir, broad_type)
  path <- file.path(dir, "subcluster.rds")
  if (!file.exists(path)) stop("no subcluster object at ", path)

  sub <- readRDS(path)
  res <- subcluster_resolutions(sub)

  marker_path <- file.path(dir, "markers.csv")
  markers <- if (file.exists(marker_path))
    data.table::fread(marker_path, data.table = FALSE) else NULL
  top <- if (!is.null(markers)) top_cluster_markers(markers, n_genes) else NULL

  cat(.ts(), " ", broad_type, ": ", ncol(sub), " cells\n", sep = "")
  print(res)

  if (is.null(markers))
    cat(.ts(), " no markers yet. pick a resolution from the umap, then",
        "finalise_subcluster()\n")

  if (plot && nrow(res)) print(DimPlot(sub, group.by = res$column, label = TRUE))

  invisible(list(obj         = sub,
                 markers     = markers,
                 top         = top,
                 dir         = dir,
                 broad_type  = broad_type,
                 resolutions = res))
}


######################################
# SET THE FINAL RESOLUTION           #
######################################

#points cluster_col at one of the stored resolutions. on a sketched subset
#the labels are re-projected to every cell, which is cheap because pca.full
#and the transfer neighbours are already cached on the object
set_subcluster_resolution <- function(object,
                                      resolution,
                                      cluster_col = "subcluster",
                                      dims        = 1:20,
                                      project     = TRUE,
                                      verbose     = TRUE) {

  sketched <- "sketch" %in% Assays(object)
  column   <- paste0(if (sketched) "sketch" else "RNA", "_snn_res.", resolution)

  if (!column %in% colnames(object@meta.data)) {
    hits <- grep(paste0("_snn_res\\.", resolution, "$"),
                 colnames(object@meta.data), value = TRUE)
    if (!length(hits))
      stop("no clustering at resolution ", resolution, ". available: ",
           paste(subcluster_resolutions(object)$resolution, collapse = ", "))
    column <- hits[1]
  }

  object@meta.data[[cluster_col]] <- as.character(object@meta.data[[column]])

  if (sketched && project) {
    k_weight <- max(5, min(50, floor(ncol(object[["sketch"]]) / 20)))
    object <- project_sketch_labels(object,
                                    label_col = cluster_col,
                                    dims      = dims,
                                    k_weight  = k_weight,
                                    verbose   = verbose)
  }

  DefaultAssay(object) <- "RNA"
  Idents(object) <- cluster_col

  if (verbose)
    cat(.ts(), " ", column, " -> ", cluster_col, ": ",
        length(unique(na.omit(object@meta.data[[cluster_col]]))),
        " clusters\n", sep = "")
  object
}


######################################
# LOCK IN A RESOLUTION AND GET MARKERS#
######################################

#set the resolution, compute markers for it, write the annotation template
#and save the object back. x is the list from load_subcluster(), or a Seurat
#object plus broad_type
finalise_subcluster <- function(x,
                                resolution,
                                broad_type  = NULL,
                                output_dir  = "./annotated_data",
                                cluster_col = "subcluster",
                                dims        = 1:20,
                                n_genes     = 15,
                                only_pos    = TRUE,
                                logfc       = 0.25,
                                min_pct     = 0.2,
                                format      = "csv",
                                save        = TRUE,
                                verbose     = TRUE) {

  if (inherits(x, "Seurat")) {
    if (is.null(broad_type))
      stop("supply broad_type when passing a Seurat object")
    x <- list(obj = x, dir = subcluster_dir(output_dir, broad_type),
              broad_type = broad_type)
  }

  sub <- set_subcluster_resolution(x$obj, resolution,
                                   cluster_col = cluster_col,
                                   dims        = dims,
                                   verbose     = verbose)

  if (verbose) cat(.ts(), " markers at resolution", resolution, "\n")

  markers <- subcluster_markers(sub,
                                cluster_col = cluster_col,
                                dir         = x$dir,
                                only_pos    = only_pos,
                                logfc       = logfc,
                                min_pct     = min_pct,
                                format      = format,
                                cache       = FALSE)

  write_cluster_template(sub, markers,
                         cluster_col = cluster_col,
                         n_genes     = n_genes,
                         path        = file.path(x$dir, "annotation.csv"))

  if (save) saveRDS(sub, file.path(x$dir, "subcluster.rds"))

  if (verbose)
    print(sort(table(sub@meta.data[[cluster_col]]), decreasing = TRUE))

  invisible(list(obj         = sub,
                 markers     = markers,
                 top         = top_cluster_markers(markers, n_genes),
                 dir         = x$dir,
                 broad_type  = x$broad_type,
                 resolutions = subcluster_resolutions(sub)))
}


######################################
# WRITE LABELS BACK TO THE PARENT    #
######################################

write_back_labels <- function(obj,
                              sub,
                              from   = "sub_type",
                              to     = "cell_type",
                              prefix = NULL,
                              cells  = NULL) {

  if (!to %in% colnames(obj@meta.data)) obj@meta.data[[to]] <- NA_character_

  cells  <- cells %||% colnames(sub)
  cells  <- intersect(cells, colnames(sub))
  if (!length(cells)) return(obj)

  labels <- as.character(sub@meta.data[cells, from])
  if (!is.null(prefix)) labels <- paste(prefix, labels, sep = " - ")

  obj@meta.data[cells, to] <- labels
  obj
}


######################################
# STAGE MISPLACED CELLS FOR RESCUE   #
######################################

#subclustering always throws off clusters that belong somewhere else:
#fibroblasts sitting inside Immune, macrophages inside Endothelial, and so
#on. label those Rem.<something> in sub_labels, then call this to move them
#into their own compartment so they can be reclustered together, free of
#the dominant cell type that was drowning them out.
#
#labels are matched on the fine part only, so "Immune - Rem.Fibro" and
#"Rem.Fibro" both count.
#
#include_unlabelled also picks up cells that never received a fine label:
#either NA, or still carrying the bare compartment name because their
#broad_type changed after their compartment's subcluster.rds was cached.
stage_rescue_compartment <- function(obj,
                                     label_col          = "cell_type",
                                     broad_col          = "broad_type",
                                     prefix             = "Rem.",
                                     rescue_as          = "Rescue",
                                     include_unlabelled = TRUE,
                                     verbose            = TRUE) {

  if (!label_col %in% colnames(obj@meta.data))
    stop(label_col, " not found, run apply_subcluster_labels() first")

  full  <- as.character(obj@meta.data[[label_col]])
  broad <- as.character(obj@meta.data[[broad_col]])
  fine  <- sub("^.*? - ", "", full, perl = TRUE)

  is_rem <- !is.na(fine) & startsWith(fine, prefix)

  is_bare <- rep(FALSE, length(full))
  if (include_unlabelled) {
    already <- !is.na(broad) & broad == rescue_as
    is_bare <- (is.na(full) |
                (!is.na(full) & !is.na(broad) & full == broad)) & !already
  }

  take <- is_rem | is_bare

  if (!any(take)) {
    if (verbose) cat(.ts(), " nothing to stage\n")
    return(obj)
  }

  if (verbose) {
    cat(.ts(), " staging ", sum(take), " cells as '", rescue_as, "': ",
        sum(is_rem), " matching '", prefix, "', ",
        sum(is_bare), " unlabelled\n", sep = "")
    source <- ifelse(is.na(broad[take]), "<NA>", broad[take])
    print(sort(table(paste0(source, " / ",
                            ifelse(is_rem[take], fine[take], "<unlabelled>"))),
               decreasing = TRUE))
  }

  obj@meta.data[[broad_col]][take] <- rescue_as

  #cell_type is rebuilt from broad_type by apply_subcluster_labels(), so
  #clear the stale labels now
  obj@meta.data[[label_col]][take] <- rescue_as

  if (verbose) print(sort(table(obj@meta.data[[broad_col]]), decreasing = TRUE))
  obj
}


######################################
# PASS 2: APPLY THE FINE LABELS      #
######################################

#labels may be
#  NULL                     read every compartment from its annotation.csv
#  a named list of named
#  character vectors        e.g. list("Fibroblast" = c("0" = "Papillary"))
#compartments missing from the list fall back to their csv
apply_subcluster_labels <- function(obj,
                                    labels        = NULL,
                                    broad_col     = "broad_type",
                                    label_col     = "cell_type",
                                    cluster_col   = "subcluster",
                                    sub_label_col = "sub_type",
                                    output_dir    = "./annotated_data",
                                    prefix        = TRUE,
                                    skip_types    = NULL,
                                    save          = TRUE,
                                    verbose       = TRUE) {

  #compartments with no subclustering keep their broad label
  obj@meta.data[[label_col]] <- as.character(obj@meta.data[[broad_col]])

  types <- setdiff(
    sort(unique(na.omit(as.character(obj@meta.data[[broad_col]])))),
    skip_types)

  for (type in types) {

    dir  <- subcluster_dir(output_dir, type)
    path <- file.path(dir, "subcluster.rds")
    if (!file.exists(path)) next

    sub  <- readRDS(path)
    here <- if (!is.null(labels)) labels[[type]] else NULL

    sub <- tryCatch(
      label_sketch_clusters(sub,
                            labels      = here,
                            cluster_col = cluster_col,
                            label_col   = sub_label_col,
                            output_dir  = output_dir,
                            template    = file.path(dir, "annotation.csv"),
                            verbose     = FALSE),
      error = function(e)
        stop("[", type, "] ", conditionMessage(e), call. = FALSE))

    #a cell that has since been moved to another compartment, e.g. by
    #stage_rescue_compartment(), must not be relabelled by its old one
    still_here <- colnames(sub)[
      which(as.character(obj@meta.data[colnames(sub), broad_col]) == type)]

    use_prefix <- if (isTRUE(prefix)) TRUE
                  else if (isFALSE(prefix)) FALSE
                  else type %in% prefix

    obj <- write_back_labels(obj, sub,
                             from   = sub_label_col,
                             to     = label_col,
                             prefix = if (use_prefix) type else NULL,
                             cells  = still_here)

    if (verbose)
      cat(.ts(), " ", type, ": ",
          length(unique(sub@meta.data[[sub_label_col]])), " labels applied to ",
          length(still_here), " of ", ncol(sub), " cells\n", sep = "")

    if (save) saveRDS(sub, path)
    rm(sub); gc(verbose = FALSE)
  }

  Idents(obj) <- label_col
  if (verbose) print(sort(table(obj@meta.data[[label_col]]), decreasing = TRUE))
  obj
}


#######################################################################
#######################################################################
#                                                                     #
#                              COLOUR MAP                             #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: cohort wide colour assignment and persistence

######################################
# COLOUR MAP                         #
######################################

#one colour per label for the whole cohort, persisted to disk
build_colour_map <- function(obj,
                             label_col        = "cell_type",
                             output_dir       = "./annotated_data",
                             spatial_contrast = TRUE,
                             seed             = 1,
                             overwrite        = TRUE) {

  path   <- file.path(output_dir, "cell_type_colours.csv")
  labels <- sort(unique(na.omit(as.character(obj@meta.data[[label_col]]))))

  if (file.exists(path) && !overwrite) {
    colours <- load_colour_map(output_dir)
    for (label in setdiff(labels, names(colours))) {
      index <- length(colours) + 1
      colours[label] <- if (index <= length(default_palette))
        default_palette[index] else hcl.colors(index, "Dark 3")[index]
    }
  } else if (spatial_contrast &&
             exists("assign_contrast_palette", mode = "function")) {
    adjacency <- NULL
    if (exists("cluster_spatial_adjacency", mode = "function") &&
        length(Images(obj))) {
      adjacency <- tryCatch(
        Reduce(merge_adjacency,
               lapply(Images(obj), function(fov)
                 cluster_spatial_adjacency(obj, group_col = label_col,
                                           image = fov, k = 12))),
        error = function(e) NULL)
    }
    colours <- assign_contrast_palette(labels, adjacency = adjacency,
                                       seed = seed, verbose = FALSE)
  } else {
    ranked  <- names(sort(table(as.character(obj@meta.data[[label_col]])),
                          decreasing = TRUE))
    colours <- setNames(vapply(seq_along(ranked), function(index)
      if (index <= length(default_palette)) default_palette[index] else
        hcl.colors(index, "Dark 3")[index], character(1)), ranked)
    colours <- colours[labels]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.frame(cell_type = names(colours), hex = unname(colours)), path)
  cat(.ts(), " colours:", path, "\n")
  colours
}

#load the persisted colour map
load_colour_map <- function(output_dir = "./annotated_data") {
  path <- file.path(output_dir, "cell_type_colours.csv")
  if (!file.exists(path)) stop("no colour map at ", path)
  table <- data.table::fread(path, data.table = FALSE)
  setNames(table$hex, table$cell_type)
}

#combine two adjacency matrices over the union of their labels
merge_adjacency <- function(a, b) {
  labels <- union(rownames(a), rownames(b))
  merged <- matrix(0, length(labels), length(labels),
                   dimnames = list(labels, labels))
  merged[rownames(a), colnames(a)] <- merged[rownames(a), colnames(a)] + a
  merged[rownames(b), colnames(b)] <- merged[rownames(b), colnames(b)] + b
  merged
}

#######################################################################
#######################################################################
#                                                                     #
#                       SPATIAL PLOTS AND EXPORT                      #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: region and cohort spatial plots, per region rds/h5ad/csv export

######################################
# SPATIAL PLOTS                      #
######################################

#cell types on the proseg segmentation of one region
plot_proseg_spatial <- function(obj,
                                group_col  = "cell_type",
                                fov        = NULL,
                                colours    = NULL,
                                output_dir = "./annotated_data",
                                size       = 0.6,
                                border     = NA,
                                dark       = FALSE,
                                legend     = TRUE,
                                title      = NULL) {

  colours <- colours %||% load_colour_map(output_dir)
  fov     <- fov %||% Images(obj)[1]
  labels  <- setdiff(unique(as.character(obj@meta.data[[group_col]])), NA)
  missing <- setdiff(labels, names(colours))
  if (length(missing)) stop("no colour for: ", paste(missing, collapse = ", "))

  ImageDimPlot(obj, fov = fov, group.by = group_col, cols = colours[labels],
               size = size, border.color = border, border.size = 0,
               dark.background = dark) +
    ggtitle(title %||% fov) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = if (legend) "right" else "none",
          legend.text = element_text(size = 8))
}

#all regions side by side from centroid coordinates
plot_cohort_spatial <- function(obj,
                                group_col  = "cell_type",
                                sample_col = "sample_id",
                                colours    = NULL,
                                size       = 0.3,
                                ncol       = 3) {

  coords <- do.call(rbind, lapply(Images(obj), function(fov) {
    region <- proseg_coords(obj, fov)
    region$fov <- fov
    region
  }))

  meta  <- obj@meta.data[rownames(coords), c(group_col, sample_col),
                         drop = FALSE]
  frame <- cbind(coords, meta)
  frame <- frame[!is.na(frame[[group_col]]), ]

  ggplot(frame, aes(x, y, colour = .data[[group_col]])) +
    geom_point(size = size, stroke = 0) +
    (if (!is.null(colours)) scale_colour_manual(values = colours) else NULL) +
    coord_fixed() +
    facet_wrap(as.formula(paste("~", sample_col)), ncol = ncol,
               scales = "free") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_void(base_size = 10)
}

######################################
# EXPORT ANNOTATED REGIONS           #
######################################

#write rds, h5ad, metadata and plots for each region
#image names are not guaranteed to equal sample_id: Seurat sanitises them
#on assignment, so "XE791-D" can land as "XE791.D". resolve the match
#rather than assuming it
.match_fov <- function(region, sample_id) {

  images <- Images(region)
  if (!length(images)) return(NULL)
  if (sample_id %in% images) return(sample_id)

  flatten <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", x))
  hit <- images[flatten(images) == flatten(sample_id)]
  if (length(hit)) return(hit[1])

  #last resort: the only image still holding cells after subsetting
  n <- vapply(images, function(i)
    length(intersect(Cells(region[[i]]), colnames(region))), integer(1))
  if (sum(n > 0) == 1) return(images[which(n > 0)])
  NULL
}

export_annotated_regions <- function(obj,
                                     output_dir = "./annotated_data",
                                     label_col  = "cell_type",
                                     sample_col = "sample_id",
                                     colours    = NULL,
                                     write_rds  = TRUE,
                                     write_h5ad = TRUE,
                                     plots      = TRUE,
                                     size       = 0.6,
                                     verbose    = TRUE) {

  colours <- colours %||% build_colour_map(obj, label_col, output_dir)
  DefaultAssay(obj) <- "RNA"
  samples <- unique(obj@meta.data[[sample_col]])
  summary_rows <- list()

  for (sample_id in samples) {

    if (verbose) cat(.ts(), " exporting", sample_id, "\n")
    region_dir <- file.path(output_dir, sample_id)
    dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)

    cells  <- colnames(obj)[obj@meta.data[[sample_col]] == sample_id]
    region <- subset(obj, cells = cells)

    fov_name <- .match_fov(region, sample_id)
    for (img in setdiff(Images(region), fov_name)) region[[img]] <- NULL

    region[["RNA"]] <- JoinLayers(region[["RNA"]])

    if (write_rds)
      saveRDS(region, file.path(region_dir,
                                paste0(sample_id, "_annotated.rds")))

    data.table::fwrite(region@meta.data,
                       file.path(region_dir, "cell_annotations.csv.gz"))

    if (write_h5ad)
      export_h5ad(region, file.path(region_dir,
                                    paste0(sample_id, "_annotated.h5ad")))

    if (plots) {
      if (is.null(fov_name)) {
        cat(.ts(), "  no image matching", sample_id,
            "- available:", paste(Images(obj), collapse = ", "), "\n")
      } else {
        pdf(file.path(region_dir, "spatial_cell_types.pdf"),
            width = 14, height = 12)
        print(plot_proseg_spatial(region, label_col, fov = fov_name,
                                  colours = colours, size = size,
                                  title = sample_id))
        dev.off()
      }
    }

    summary_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      n_cells   = ncol(region),
      labelled  = sum(!is.na(region@meta.data[[label_col]])),
      n_types   = length(setdiff(unique(region@meta.data[[label_col]]), NA)),
      stringsAsFactors = FALSE)

    rm(region); gc(verbose = FALSE)
  }

  summary_table <- do.call(rbind, summary_rows)
  data.table::fwrite(summary_table,
                     file.path(output_dir, "annotation_summary.csv"))
  if (verbose) print(summary_table)
  invisible(summary_table)
}

#write a Seurat object to h5ad with spatial coordinates
export_h5ad <- function(obj, path, assay = "RNA", layer = "counts") {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cat(.ts(), " reticulate unavailable, skipping h5ad\n")
    return(invisible(FALSE))
  }

  written <- tryCatch({
    anndata <- reticulate::import("anndata", convert = FALSE)
    sparse  <- reticulate::import("scipy.sparse", convert = FALSE)

    counts <- as(LayerData(obj, assay = assay, layer = layer), "dgCMatrix")
    matrix <- sparse$csr_matrix(reticulate::r_to_py(Matrix::t(counts)))

    meta <- obj@meta.data
    meta[] <- lapply(meta, function(column)
      if (is.factor(column)) as.character(column) else column)

    adata <- anndata$AnnData(
      X   = matrix,
      obs = reticulate::r_to_py(meta),
      var = reticulate::r_to_py(data.frame(row.names = rownames(counts))))

    coords <- proseg_coords(obj, cells = colnames(counts))
    if (!is.null(coords) && nrow(coords) == ncol(counts))
      adata$obsm$update(reticulate::dict(
        spatial = reticulate::r_to_py(as.matrix(coords[colnames(counts), ]))))

    adata$write_h5ad(path)
    TRUE
  }, error = function(e) {
    cat(.ts(), " h5ad failed:", conditionMessage(e), "\n")
    FALSE
  })

  if (isTRUE(written)) cat(.ts(), " h5ad:", path, "\n")
  invisible(written)
}

#######################################################################
#######################################################################
#                                                                     #
#                      REFERENCE BASED ANNOTATION                     #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: label transfer from an external reference, RCTD and SingleR workflows

######################################
# SEURAT LABEL TRANSFER              #
######################################

annotate_proseg_seurat_10x <- function(
    proseg_dir = "./proseg_results",
    reference_dir = "./reference",
    ref_label_col = "cell_type",
    output_dir = "./annotated_data",
    ref_method = "pca",
    dims = 1:30,
    k.weight = 50,
    sketch_ncells = 15000,
    plots = TRUE,
    overwrite = TRUE,
    replot_only = FALSE
) {
  
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  library(ggpmisc)
  library(scales)
  
  ref_method <- match.arg(ref_method, c("pca", "cca"))
  
  ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  
  # --- Stable colour map --------------------------------------------------
  base_palette <- c(
    "#2166AC", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#17BECF", "#E377C2", "#8C564B", "#BCBD22", "#1F77B4",
    "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
    "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5", "#393B79",
    "#637939", "#8C6D31", "#843C39", "#7B4173", "#5254A3",
    "#6B6ECF", "#9C9EDE", "#E7BA52", "#BD9E39", "#AD494A",
    "#D6616B", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
    "#E6550D", "#FD8D3C", "#31A354", "#74C476", "#756BB1",
    "#FDAE6B", "#A1D99B", "#DADAEB", "#636363", "#969696",
    "#525252", "#FDD0A2", "#C7E9C0", "#B3CDE3", "#FDDBC7"
  )
  
  global_cmap <- list()
  
  assign_colours <- function(cell_types) {
    for (ct in cell_types) {
      if (is.null(global_cmap[[ct]])) {
        idx <- length(global_cmap) + 1
        if (idx <= length(base_palette)) {
          global_cmap[[ct]] <<- base_palette[idx]
        } else {
          global_cmap[[ct]] <<- hcl.colors(idx, palette = "Dark 3")[idx]
        }
      }
    }
    return(setNames(
      vapply(cell_types, function(ct) global_cmap[[ct]], character(1)),
      cell_types
    ))
  }
  
  cat(ts(), "============================================================\n")
  cat(ts(), " 10X Xenium Annotation + Reference QC (Memory-Optimised)\n")
  cat(ts(), " Method:", ref_method, "\n")
  cat(ts(), "============================================================\n\n")
  
  ##########################################################################
  # Discover regions
  ##########################################################################
  rds_hits <- list.files(
    proseg_dir,
    pattern = "_proseg_seurat_mask\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(rds_hits) == 0) stop("[ERROR] No *_proseg_seurat.rds found.")
  
  regions <- list()
  for (p in rds_hits) {
    nm <- basename(dirname(p))
    regions[[nm]] <- list(rds_path = p)
    cat(ts(), " Found:", nm, "\n")
  }
  
  ##########################################################################
  # Load reference — BPCells on-disk counts, then preprocess
  ##########################################################################
  ref_path <- file.path(reference_dir, "consensus_reference.rds")
  cat(ts(), "[REF] Loading:", ref_path, "\n")
  reference_obj <- readRDS(ref_path)
  DefaultAssay(reference_obj) <- "RNA"
  
  # Keep an in-memory sparse copy of the original counts BEFORE BPCells
  # replaces them.  The BPCells cache can silently return corrupt data
  # when row-subsetted (stale cache, version mismatch, etc.).  This copy
  # is used inside run_10x_annotation to build the per-region ref_sub.
  ref_counts_inmem <- as(reference_obj[["RNA"]]$counts, "dgCMatrix")
  cat(ts(), "[REF] In-memory counts snapshot:",
      nrow(ref_counts_inmem), "genes x",
      ncol(ref_counts_inmem), "cells\n")
  
  bpcells_dir <- file.path(reference_dir, "ref_counts_bpcells")
  if (!dir.exists(bpcells_dir)) {
    cat(ts(), "[REF] Writing counts to disk (BPCells):", bpcells_dir, "\n")
    write_matrix_dir(
      mat = reference_obj[["RNA"]]$counts,
      dir = bpcells_dir
    )
  } else {
    cat(ts(), "[REF] BPCells cache found:", bpcells_dir, "\n")
  }
  
  counts_on_disk <- open_matrix_dir(dir = bpcells_dir)
  reference_obj[["RNA"]]$counts <- counts_on_disk
  rm(counts_on_disk)
  gc(verbose = FALSE)
  
  if (!replot_only) {
    cat(ts(), "[REF] Counts now on-disk. Preprocessing...\n")
    
    reference_obj <- NormalizeData(reference_obj, verbose = FALSE)
    reference_obj <- FindVariableFeatures(reference_obj, verbose = FALSE)
    reference_obj <- ScaleData(reference_obj, verbose = FALSE)
    reference_obj <- RunPCA(reference_obj, verbose = FALSE)
    
    reference_obj[["RNA"]]$scale.data <- NULL
    gc(verbose = FALSE)
    
    cat(ts(), "[REF] Reference preprocessed. PCA stored; scale.data dropped.\n")
  } else {
    cat(ts(), "[REF] replot_only = TRUE. Skipping preprocessing.\n")
  }
  
  ref_labels   <- reference_obj[[ref_label_col]][, 1]
  ref_genes    <- rownames(reference_obj)
  ref_metadata <- reference_obj@meta.data
  
  ref_type_counts <- sort(table(ref_labels), decreasing = TRUE)
  ranked_types    <- names(ref_type_counts)
  cat(ts(), "[CMAP] Pre-assigning colours for", length(ranked_types),
      "reference types (ranked by abundance)\n")
  for (i in seq_along(ranked_types)) {
    ct <- ranked_types[i]
    if (i <= length(base_palette)) {
      global_cmap[[ct]] <- base_palette[i]
    } else {
      global_cmap[[ct]] <- hcl.colors(i, palette = "Dark 3")[i]
    }
  }
  cat(ts(), "[CMAP] Top 5:",
      paste(sprintf("%s (%s, n=%d)",
                    ranked_types[1:min(5, length(ranked_types))],
                    vapply(ranked_types[1:min(5, length(ranked_types))],
                           function(ct) global_cmap[[ct]], character(1)),
                    ref_type_counts[1:min(5, length(ranked_types))]),
            collapse = ", "), "\n")
  
  ##########################################################################
  # Robust HVG selection with VST-NaN fallback
  #
  # Proseg outputs fractional expected counts. The VST method in
  # FindVariableFeatures fits a loess curve to log10(mean) vs log10(var)
  # assuming Poisson/NB count data.  Fractional counts violate this
  # assumption and produce NaN standardised variances, which corrupts
  # the entire HVG ranking and can cause downstream ScaleData/RunPCA
  # failures.
  #
  # This wrapper tries VST first (for compatibility with integer-count
  # references), checks whether the result contains NaN, and falls back
  # to selection.method = "dispersion" if so.
  ##########################################################################
  .safe_find_hvgs <- function(obj, assay, ts_fn = NULL) {
    
    # Attempt default VST, capturing warnings
    hvg_warnings <- NULL
    obj <- withCallingHandlers(
      FindVariableFeatures(obj, assay = assay, verbose = FALSE),
      warning = function(w) {
        hvg_warnings <<- c(hvg_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    # Check for NaN in the standardised variance column.
    # HVFInfo can differ across Seurat versions, so wrap defensively.
    vst_failed <- FALSE
    tryCatch({
      hvf_info <- HVFInfo(obj, assay = assay, method = "vst")
      if ("variance.standardized" %in% colnames(hvf_info)) {
        n_nan   <- sum(is.nan(hvf_info$variance.standardized))
        n_total <- nrow(hvf_info)
        vst_failed <- (n_nan > n_total * 0.5)
      }
    }, error = function(e) {
      # HVFInfo not available or method mismatch — check via warnings
      vst_failed <<- any(grepl("NaN", hvg_warnings, fixed = TRUE))
    })
    
    if (vst_failed) {
      if (!is.null(ts_fn)) {
        cat(ts_fn(), "  [HVG] VST produced degenerate results.",
            "Retrying with method = 'dispersion'\n")
      }
      obj <- FindVariableFeatures(
        obj, assay = assay,
        selection.method = "dispersion",
        verbose = FALSE
      )
    } else {
      # Replay any non-NaN warnings
      for (w in hvg_warnings) {
        if (!grepl("NaN", w, fixed = TRUE) &&
            !grepl("not a multiple of replacement length", w, fixed = TRUE))
          warning(w, call. = FALSE)
      }
    }
    
    return(obj)
  }
  
  ##########################################################################
  # Annotation with sketch workflow
  ##########################################################################
  run_10x_annotation <- function(obj, region_name) {
    
    DefaultAssay(obj) <- "RNA"
    
    # ------------------------------------------------------------------
    # 1. Establish the shared gene set (non-zero in BOTH datasets)
    # ------------------------------------------------------------------
    query_counts_raw <- LayerData(obj, assay = "RNA", layer = "counts")
    query_genes_nz   <- rownames(query_counts_raw)[rowSums(query_counts_raw) > 0]
    
    shared_genes <- intersect(query_genes_nz, ref_genes)
    if (length(shared_genes) < 50)
      stop("[ERROR] Too few shared genes: ", length(shared_genes))
    
    # Also require non-zero in reference (some panel genes may be absent)
    ref_counts_shared <- ref_counts_inmem[shared_genes, ]
    ref_nz <- rowSums(ref_counts_shared) > 0
    shared_genes <- shared_genes[ref_nz]
    rm(ref_counts_shared, ref_nz); gc(verbose = FALSE)
    
    n_cells <- ncol(obj)
    cat(ts(), "  Shared non-zero genes:", length(shared_genes),
        "| Cells:", n_cells, "\n")
    
    # ------------------------------------------------------------------
    # 2. Build a CLEAN query object for annotation.
    #    Strip FOV/spatial data — the "Not validating" warnings from
    #    Seurat's spatial slots cause cascading layer-state corruption
    #    that makes ScaleData/RunPCA silently produce degenerate output.
    #    Predictions are transferred back to the spatial object later.
    # ------------------------------------------------------------------
    query_counts <- as(query_counts_raw[shared_genes, ], "dgCMatrix")
    rm(query_counts_raw)
    
    use_sketch <- n_cells > sketch_ncells * 1.5
    query_bp_dir <- file.path(output_dir, ".bpcells_tmp", region_name)
    
    if (use_sketch) {
      # ---- SKETCH PATH ----
      dir.create(query_bp_dir, recursive = TRUE, showWarnings = FALSE)
      write_matrix_dir(mat = query_counts, dir = query_bp_dir, overwrite = TRUE)
      query_counts_disk <- open_matrix_dir(dir = query_bp_dir)
      
      query_obj <- CreateSeuratObject(counts = query_counts_disk)
      rm(query_counts, query_counts_disk)
      gc(verbose = FALSE)
      cat(ts(), "  Query on disk (BPCells) for sketching\n")
      
      query_obj <- NormalizeData(query_obj, verbose = FALSE)
      # LeverageScore needs VariableFeatures; use all shared genes.
      VariableFeatures(query_obj) <- shared_genes
      query_obj <- FindVariableFeatures(query_obj, verbose = FALSE)
      
      cat(ts(), "  Sketching", sketch_ncells, "cells...\n")
      query_obj <- SketchData(
        object = query_obj, ncells = sketch_ncells,
        method = "LeverageScore", sketched.assay = "sketch",
        verbose = FALSE
      )
      DefaultAssay(query_obj) <- "sketch"
      query_assay <- "sketch"
      
      # Sketch may zero-out some genes; drop them.
      sk_keep <- rowSums(LayerData(query_obj, assay = "sketch",
                                   layer = "counts")) > 0
      sketch_genes <- names(sk_keep)[sk_keep]
      shared_use <- intersect(sketch_genes, shared_genes)
      
      query_obj <- NormalizeData(query_obj, assay = "sketch", verbose = FALSE)
      
    } else {
      # ---- IN-MEMORY PATH ----
      query_obj <- CreateSeuratObject(counts = query_counts)
      rm(query_counts)
      gc(verbose = FALSE)
      cat(ts(), "  Clean query object built (no FOV, in-memory)\n")
      
      query_assay <- "RNA"
      query_obj <- NormalizeData(query_obj, verbose = FALSE)
      shared_use <- shared_genes
    }
    
    # ------------------------------------------------------------------
    # 3. Prepare the query for PCA on ALL shared genes
    #    (panel data: every gene is a curated marker — skip HVG selection)
    # ------------------------------------------------------------------
    VariableFeatures(query_obj, assay = query_assay) <- shared_use
    
    query_obj <- ScaleData(query_obj, assay = query_assay,
                           features = shared_use, verbose = FALSE)
    cat(ts(), "  Query ScaleData done on", length(shared_use), "genes\n")
    
    query_obj <- RunPCA(query_obj, assay = query_assay,
                        features = shared_use,
                        npcs = min(max(dims), 50), verbose = FALSE)
    cat(ts(), "  Query PCA done\n")
    if (use_sketch) cat(ts(), "  (sketch PCA)\n")
    
    # ------------------------------------------------------------------
    # 4. Build the reference sub on the SAME shared gene set
    #    Uses the in-memory counts snapshot (ref_counts_inmem) to bypass
    #    the BPCells cache, which can return corrupt data on row subsets.
    # ------------------------------------------------------------------
    ref_counts_sub <- ref_counts_inmem[shared_use, ]
    ref_sub <- CreateSeuratObject(
      counts    = ref_counts_sub,
      meta.data = ref_metadata
    )
    rm(ref_counts_sub)
    
    ref_sub <- NormalizeData(ref_sub, verbose = FALSE)
    VariableFeatures(ref_sub) <- shared_use
    ref_sub <- ScaleData(ref_sub, features = shared_use, verbose = FALSE)
    ref_sub <- RunPCA(ref_sub, features = shared_use, verbose = FALSE)
    cat(ts(), "  Reference sub PCA done on", length(shared_use), "genes\n")
    
    if (ref_method != "cca") ref_sub[["RNA"]]$scale.data <- NULL
    
    # Clamp dims
    max_pc   <- ncol(Embeddings(ref_sub, "pca"))
    use_dims <- dims[dims <= max_pc]
    if (length(use_dims) < length(dims))
      cat(ts(), "  Clamped dims to 1:", max(use_dims),
          "(reference has", max_pc, "PCs)\n")
    
    # ------------------------------------------------------------------
    # 5. Find transfer anchors
    # ------------------------------------------------------------------
    # Materialise any lazy layers for anchor-finding.
    query_obj[[query_assay]]$counts <- as(
      query_obj[[query_assay]]$counts, "dgCMatrix"
    )
    if (ref_method == "cca") {
      query_obj[[query_assay]]$data <- as(
        query_obj[[query_assay]]$data, "dgCMatrix"
      )
      if ("scale.data" %in% Layers(query_obj[[query_assay]]))
        query_obj[[query_assay]]$scale.data <- NULL
    }
    
    cat(ts(), "  Finding transfer anchors (", ref_method, ")...\n")
    if (ref_method == "cca") {
      anchors <- FindTransferAnchors(
        reference            = ref_sub,
        query                = query_obj,
        query.assay          = query_assay,
        normalization.method = "LogNormalize",
        reduction            = "cca",
        features             = shared_use,
        dims                 = use_dims,
        verbose              = FALSE
      )
    } else {
      anchors <- FindTransferAnchors(
        reference            = ref_sub,
        query                = query_obj,
        query.assay          = query_assay,
        normalization.method = "LogNormalize",
        reduction            = "pcaproject",
        reference.reduction  = "pca",
        features             = shared_use,
        dims                 = use_dims,
        verbose              = FALSE
      )
    }
    
    # ------------------------------------------------------------------
    # 6. Transfer labels
    # ------------------------------------------------------------------
    cat(ts(), "  Transferring labels...\n")
    wr <- if (ref_method == "cca") "cca" else "pcaproject"
    
    n_anchors <- tryCatch(nrow(slot(anchors, "anchors")),
                          error = function(e) NA_integer_)
    kw <- k.weight
    if (!is.na(n_anchors)) kw <- min(k.weight, max(5L, n_anchors - 1L))
    if (kw < k.weight)
      cat(ts(), "  Clamped k.weight", k.weight, "->", kw,
          "(", n_anchors, "anchors)\n")
    
    label_transfer <- TransferData(
      anchorset        = anchors,
      refdata          = ref_sub[[ref_label_col]][, 1],
      weight.reduction = wr,
      dims             = use_dims,
      k.weight         = kw,
      verbose          = FALSE
    )
    
    rm(anchors, ref_sub); gc(verbose = FALSE)
    
    # ------------------------------------------------------------------
    # 7. Sketch projection (if applicable) or direct extraction
    # ------------------------------------------------------------------
    if (use_sketch) {
      query_obj <- AddMetaData(query_obj, label_transfer)
      
      cat(ts(), "  Projecting sketch labels to full dataset...\n")
      query_obj[["RNA"]]$counts <- as(
        query_obj[["RNA"]]$counts, "dgCMatrix"
      )
      
      DefaultAssay(query_obj) <- "sketch"
      query_obj <- ProjectData(
        object             = query_obj,
        assay              = "RNA",
        full.reduction     = "pca.full",
        sketched.assay     = "sketch",
        sketched.reduction = "pca",
        dims               = use_dims,
        verbose            = FALSE
      )
      
      DefaultAssay(query_obj) <- "RNA"
      query_obj <- TransferSketchLabels(
        query_obj,
        sketched.assay      = "sketch",
        reduction           = "pca.full",
        dims                = use_dims,
        refdata             = list(predicted.id_full = "predicted.id"),
        k                   = 50,
        recompute.neighbors = FALSE,
        recompute.weights   = FALSE,
        verbose             = TRUE
      )
      
      all_cells  <- colnames(obj)
      pred_id    <- query_obj$predicted.id_full[all_cells]
      pred_score <- query_obj$predicted.id_full.score[all_cells]
      
    } else {
      all_cells <- colnames(obj)
      # label_transfer is a data.frame whose ROW NAMES are the query barcodes.
      # Use [rows, col] indexing — `$col[keys]` returns an UNNAMED vector and
      # silently yields all-NA when subset by barcode.
      pred_id    <- as.character(label_transfer[all_cells, "predicted.id"])
      score_cols <- grep("^prediction\\.score\\.", colnames(label_transfer),
                         value = TRUE)
      pred_score <- apply(label_transfer[all_cells, score_cols, drop = FALSE],
                          1, max)
    }
    
    pred_meta <- data.frame(
      predicted_cell_type  = pred_id,
      prediction_score_max = pred_score,
      row.names = all_cells
    )
    
    unlink(query_bp_dir, recursive = TRUE)
    rm(query_obj); gc(verbose = FALSE)
    
    return(pred_meta)
  }
  
  ##########################################################################
  # QC plotting — v5-safe layer access
  ##########################################################################
  getCellMeans <- function(celltype, obj) {
    
    xen_counts <- tryCatch(
      LayerData(obj, assay = "RNA", layer = "counts"),
      error = function(e) { return(NULL) }
    )
    if (is.null(xen_counts)) return(NULL)
    
    ref_counts <- tryCatch(
      LayerData(reference_obj, assay = "RNA", layer = "counts"),
      error = function(e) { return(NULL) }
    )
    if (is.null(ref_counts)) return(NULL)
    
    common_genes <- intersect(rownames(xen_counts), rownames(ref_counts))
    
    xen_cells <- colnames(obj)[obj$predicted_cell_type == celltype]
    ref_cells <- colnames(reference_obj)[ref_labels == celltype]
    
    if (length(xen_cells) < 10 || length(ref_cells) < 10) return(NULL)
    
    xen_cells <- intersect(xen_cells, colnames(xen_counts))
    ref_cells <- intersect(ref_cells, colnames(ref_counts))
    
    if (length(xen_cells) < 10 || length(ref_cells) < 10) return(NULL)
    
    xen_mat <- xen_counts[common_genes, xen_cells, drop = FALSE]
    ref_mat <- ref_counts[common_genes, ref_cells, drop = FALSE]
    
    df <- data.frame(
      Xenium    = rowMeans(xen_mat),
      Reference = rowMeans(ref_mat)
    )
    
    df <- df[df$Xenium > 0 & df$Reference > 0, ]
    return(df)
  }
  
  plotCor <- function(celltype, obj) {
    
    df <- getCellMeans(celltype, obj)
    if (is.null(df) || nrow(df) < 20) return(NULL)
    
    fit <- lm(log10(Xenium) ~ log10(Reference), data = df)
    r2  <- summary(fit)$r.squared
    
    p <- ggplot(df, aes(x = Reference, y = Xenium)) +
      geom_point(size = 0.8, alpha = 0.5, color = "steelblue") +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      stat_poly_eq() +
      scale_x_log10(labels = label_log(digits = 1)) +
      scale_y_log10(labels = label_log(digits = 1)) +
      xlab("Reference") +
      ylab("Xenium") +
      theme_bw(base_size = 8)
    
    title <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste0(celltype, " (R\u00B2=", round(r2, 2), ")"),
        x = 0, hjust = 0, size = 10
      )
    
    cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
  }
  
  ##########################################################################
  # Per-region processing
  ##########################################################################
  summaries <- list()
  
  for (nm in names(regions)) {
    cat(ts(), "\n----------------------------------------\n")
    cat(ts(), " REGION:", nm, "\n")
    
    info <- regions[[nm]]
    if (!replot_only) {
      obj <- readRDS(info$rds_path)
    }
    
    region_out <- file.path(output_dir, nm)
    dir.create(region_out, recursive = TRUE, showWarnings = FALSE)
    output_rds   <- file.path(region_out, paste0(nm, "_annotated.rds"))
    metadata_csv <- file.path(region_out, "cell_annotations.csv.gz")
    plot_pdf     <- file.path(region_out, "reference_correlation_plots.pdf")
    
    if (!overwrite && file.exists(output_rds) && !replot_only) {
      cat(ts(), " Skipping (already annotated)\n")
      next
    }
    
    if (replot_only) {
      if (!file.exists(output_rds)) {
        cat(ts(), " [WARN] No annotated RDS found at", output_rds,
            "-- skipping region\n")
        next
      }
      cat(ts(), " Loading annotated object for replot...\n")
      obj <- readRDS(output_rds)
      Idents(obj) <- "predicted_cell_type"
    } else {
      pred_meta <- run_10x_annotation(obj, region_name = nm)
      
      obj <- AddMetaData(obj, pred_meta)
      
      cat(ts(), "  FOVs in annotated object:",
          paste(Images(obj), collapse = ", "), "\n")
      
      Idents(obj) <- "predicted_cell_type"
      
      saveRDS(obj, output_rds)
      fwrite(obj@meta.data, metadata_csv)
    }
    
    region_types <- unique(obj$predicted_cell_type)
    region_types <- region_types[!is.na(region_types)]
    cmap <- assign_colours(region_types)
    
    if (plots) {
      
      cat(ts(), "[QC] Generating correlation plots...\n")
      plist <- lapply(region_types, function(ct) plotCor(ct, obj))
      plist <- Filter(Negate(is.null), plist)
      if (length(plist) > 0) {
        pdf(plot_pdf, width = 12, height = 10)
        grid.arrange(grobs = plist, nrow = 4)
        dev.off()
      }
      
      cat(ts(), "[QC] Generating spatial cell type plot...\n")
      spatial_pdf <- file.path(region_out, "spatial_cell_types.pdf")
      
      avail_fovs <- Images(obj)
      if (length(avail_fovs) > 0) {
        fov_name <- avail_fovs[1]
        
        p_spatial <- ImageDimPlot(
          obj,
          fov          = fov_name,
          group.by     = "predicted_cell_type",
          cols         = cmap,
          boundaries   = "segmentation",
          border.color = NA,
          border.size  = 0,
          dark.background = FALSE
        ) +
          ggtitle(paste0(nm, " \u2014 Predicted Cell Types")) +
          theme(
            plot.title  = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 8)
          )
        
        pdf(spatial_pdf, width = 14, height = 12)
        print(p_spatial)
        dev.off()
        
        cat(ts(), "  Spatial plot saved:", spatial_pdf, "\n")
      } else {
        cat(ts(), "  [WARN] No FOV found \u2014 skipping spatial plot.\n")
      }
    }
    
    summaries[[nm]] <- data.frame(
      region      = nm,
      total_cells = ncol(obj),
      labeled     = sum(!is.na(obj$predicted_cell_type)),
      n_types     = length(unique(obj$predicted_cell_type))
    )
    
    rm(obj)
    if (exists("pred_meta", inherits = FALSE)) rm(pred_meta)
    gc(verbose = FALSE)
    cat(ts(), " Region", nm, "complete. Memory released.\n")
  }
  
  unlink(file.path(output_dir, ".bpcells_tmp"), recursive = TRUE)
  
  if (length(summaries) > 0) {
    summary_df <- do.call(rbind, summaries)
    fwrite(summary_df, file.path(output_dir, "annotation_10x_summary.csv"))
  }
  
  if (length(global_cmap) > 0) {
    cmap_df <- data.frame(
      cell_type = names(global_cmap),
      hex       = unlist(global_cmap),
      stringsAsFactors = FALSE
    )
    cmap_path <- file.path(output_dir, "cell_type_colour_map.csv")
    fwrite(cmap_df, cmap_path)
    cat(ts(), "[CMAP] Colour map saved:", cmap_path, "\n")
  }
  
  cat(ts(), "\n Annotation complete")
}

######################################
# RCTD AND SINGLER                   #
######################################

annotate_proseg_seurat <- function(proseg_dir,
                                   reference_dir,
                                   # Reference label column — v1 hardcoded "consensus_label"
                                   ref_label_col      = "consensus_label",
                                   # Method selection: "both", "singler", "rctd"
                                   annotation_method  = "both",
                                   # QC
                                   xen_min_features = 5,
                                   xen_min_counts   = 10,
                                   coord_x          = "x",
                                   coord_y          = "y",
                                   # Pre-clustering for scalable SingleR
                                   cluster_resolution = 3.0,
                                   cluster_npcs       = 30,
                                   cluster_dims       = 1:20,
                                   # SingleR
                                   singler_quantile  = 0.8,
                                   singler_fine_tune = FALSE,
                                   singler_de_method = "wilcox",
                                   # RCTD
                                   rctd_mode              = "full",
                                   rctd_max_cores         = 32,
                                   rctd_gene_cutoff       = 0.000125,
                                   rctd_fc_cutoff         = 0.5,
                                   rctd_fc_cutoff_reg     = 0.75,
                                   rctd_UMI_min           = 0,
                                   rctd_counts_MIN        = 0,
                                   rctd_CELL_MIN_INSTANCE = 10,
                                   # Consensus
                                   singler_delta_floor            = 0.05,
                                   singler_delta_pctile_threshold = 0.95,
                                   rctd_weight_threshold          = 0.70,
                                   rctd_confident_classes         = c("singlet", "doublet_certain"),
                                   # Caching
                                   use_cache = TRUE,
                                   overwrite = FALSE) {
  
  # Validate annotation_method
  annotation_method <- match.arg(annotation_method, c("both", "singler", "rctd"))
  run_singler <- annotation_method %in% c("both", "singler")
  run_rctd    <- annotation_method %in% c("both", "rctd")
  
  config <- as.list(environment())
  config <- config[!names(config) %in% c("proseg_dir", "reference_dir", "overwrite")]
  
  ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  
  cat(ts(), "Annotation method:", annotation_method, "\n")
  cat(ts(), "  Run SingleR:", run_singler, "| Run RCTD:", run_rctd, "\n")
  cat(ts(), "  Reference label column:", ref_label_col, "\n")
  cat(ts(), "  RCTD mode:", rctd_mode, "| cores:", rctd_max_cores, "\n")
  cat(ts(), "  Cluster resolution:", cluster_resolution, "\n\n")
  
  # ============================================================================
  # DISCOVER BUILT SEURAT OBJECTS
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Discovering built Seurat objects in:", proseg_dir, "\n")
  cat(ts(), "================================================================\n\n")
  
  rds_hits <- list.files(
    proseg_dir,
    pattern    = "_proseg_seurat\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(rds_hits) == 0) {
    stop("[ERROR] No *_proseg_seurat.rds files found under: ", proseg_dir,
         "\n  Run build_proseg_seurat() first.")
  }
  
  regions <- list()
  for (rds_path in rds_hits) {
    sample_dir  <- dirname(rds_path)
    region_name <- basename(sample_dir)
    regions[[region_name]] <- list(
      dir      = sample_dir,
      rds_path = rds_path
    )
    cat(ts(), "  Found:", region_name, "->", rds_path, "\n")
  }
  
  cat(ts(), "\nDiscovered", length(regions), "regions:",
      paste(names(regions), collapse = ", "), "\n\n")
  
  # ============================================================================
  # PANEL GENES (union across regions)
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Discovering Xenium panel genes\n")
  cat(ts(), "================================================================\n\n")
  
  panel_genes <- character(0)
  for (nm in names(regions)) {
    gm_path <- file.path(regions[[nm]]$dir, "gene-metadata.csv.gz")
    if (file.exists(gm_path)) {
      gm <- fread(gm_path)
      region_genes <- gm$gene
    } else {
      cat(ts(), " ", nm, ": gene-metadata.csv.gz not found, reading RDS header...\n")
      tmp <- readRDS(regions[[nm]]$rds_path)
      region_genes <- rownames(tmp)
      rm(tmp); gc(verbose = FALSE)
    }
    cat(ts(), " ", nm, ":", length(region_genes), "genes\n")
    panel_genes <- union(panel_genes, region_genes)
  }
  
  panel_genes <- sort(unique(panel_genes))
  cat(ts(), "\nPanel gene union:", length(panel_genes), "genes\n\n")
  
  # ============================================================================
  # REFERENCE + PANEL-SPECIFIC TRAINING
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Building panel-specific reference models\n")
  cat(ts(), "================================================================\n\n")
  
  cache_singler <- file.path(proseg_dir, "singler_trained_panel.rds")
  cache_rctd    <- file.path(proseg_dir, "rctd_reference_panel.rds")
  cache_genes   <- file.path(proseg_dir, "panel_genes.txt")
  
  # Determine what's cached
  genes_match <- FALSE
  if (config$use_cache && file.exists(cache_genes)) {
    cached_genes <- readLines(cache_genes)
    genes_match  <- identical(sort(cached_genes), sort(panel_genes))
    if (!genes_match) {
      cat(ts(), "[CACHE] Panel genes changed. Retraining needed models.\n")
    }
  }
  
  singler_cached <- genes_match && file.exists(cache_singler)
  rctd_cached    <- genes_match && file.exists(cache_rctd)
  
  # Load or train SingleR
  singler_trained <- NULL
  if (run_singler) {
    if (singler_cached) {
      singler_trained <- readRDS(cache_singler)
      cat(ts(), "[CACHE] SingleR panel model loaded.\n")
    }
  }
  
  # Load or prepare RCTD ref
  rctd_ref <- NULL
  if (run_rctd) {
    if (rctd_cached) {
      rctd_ref <- readRDS(cache_rctd)
      cat(ts(), "[CACHE] RCTD panel components loaded.\n")
    }
  }
  
  # If anything still needs training, load the reference
  need_singler_train <- run_singler && is.null(singler_trained)
  need_rctd_train    <- run_rctd && is.null(rctd_ref)
  
  if (need_singler_train || need_rctd_train) {
    ref_path <- file.path(reference_dir, "consensus_reference.rds")
    cat(ts(), "[REF] Loading full reference:", ref_path, "\n")
    ref_obj <- readRDS(ref_path)
    
    # --- Validate ref_label_col exists in reference ---
    if (!ref_label_col %in% colnames(ref_obj@meta.data)) {
      avail_cols <- colnames(ref_obj@meta.data)
      stop("[ERROR] ref_label_col '", ref_label_col,
           "' not found in reference metadata.\n",
           "  Available columns: ", paste(avail_cols, collapse = ", "))
    }
    
    cat(ts(), "[REF] Full reference:", ncol(ref_obj), "cells,",
        nrow(ref_obj), "genes,",
        length(unique(ref_obj@meta.data[[ref_label_col]])), "types (",
        ref_label_col, ")\n")
    
    shared_genes <- intersect(rownames(ref_obj), panel_genes)
    panel_only   <- setdiff(panel_genes, rownames(ref_obj))
    
    cat(ts(), "[REF] Panel genes in reference:", length(shared_genes), "/",
        length(panel_genes), "\n")
    if (length(panel_only) > 0) {
      cat(ts(), "[REF] Panel genes NOT in reference (", length(panel_only),
          "):", paste(head(panel_only, 20), collapse = ", "),
          if (length(panel_only) > 20) " ..." else "", "\n")
    }
    if (length(shared_genes) < 50) {
      stop("[ERROR] Only ", length(shared_genes),
           " panel genes found in reference.")
    }
    
    cat(ts(), "[REF] Subsetting reference to", length(shared_genes), "panel genes...\n")
    ref_panel <- subset(ref_obj, features = shared_genes)
    ref_panel <- ensure_clean_assay(ref_panel, "ref_panel")
    ref_panel <- NormalizeData(ref_panel, verbose = FALSE)
    rm(ref_obj); gc(verbose = FALSE)
    
    if (need_singler_train) {
      cat(ts(), "[SINGLER] Training on", length(shared_genes), "panel genes...\n")
      ref_counts_mat <- GetAssayData(ref_panel, assay = "RNA", layer = "counts")
      ref_data_mat   <- GetAssayData(ref_panel, assay = "RNA", layer = "data")
      ref_sce <- SingleCellExperiment(
        assays  = list(counts = ref_counts_mat, logcounts = ref_data_mat),
        colData = ref_panel@meta.data
      )
      rm(ref_counts_mat, ref_data_mat)
      
      singler_trained <- trainSingleR(
        ref       = ref_sce,
        labels    = ref_sce[[ref_label_col]],
        de.method = config$singler_de_method,
        BPPARAM   = MulticoreParam(config$rctd_max_cores)
      )
      cat(ts(), "[SINGLER] Training complete.\n")
      saveRDS(singler_trained, cache_singler)
      rm(ref_sce); gc(verbose = FALSE)
    }
    
    if (need_rctd_train) {
      cat(ts(), "[RCTD] Building panel-gene reference components...\n")
      ref_counts <- GetAssayData(ref_panel, assay = "RNA", layer = "counts")
      ref_types  <- setNames(
        factor(ref_panel@meta.data[[ref_label_col]]),
        colnames(ref_panel)
      )
      ref_numi <- setNames(colSums(ref_counts), colnames(ref_panel))
      
      rctd_ref <- list(
        counts     = ref_counts,
        cell_types = ref_types,
        nUMI       = ref_numi
      )
      
      cat(ts(), "[RCTD] Panel reference:", ncol(ref_counts), "cells,",
          nrow(ref_counts), "genes,",
          length(levels(ref_types)), "types\n")
      saveRDS(rctd_ref, cache_rctd)
      rm(ref_counts); gc(verbose = FALSE)
    }
    
    writeLines(panel_genes, cache_genes)
    rm(ref_panel); gc(verbose = FALSE)
  }
  
  # RCTD drop warning
  if (run_rctd) {
    ref_type_counts <- table(rctd_ref$cell_types)
    rctd_would_drop <- names(ref_type_counts[ref_type_counts < config$rctd_CELL_MIN_INSTANCE])
    if (length(rctd_would_drop) > 0) {
      cat(ts(), "[WARN] RCTD CELL_MIN_INSTANCE =", config$rctd_CELL_MIN_INSTANCE,
          "will drop", length(rctd_would_drop), "types:",
          paste(rctd_would_drop, collapse = ", "), "\n")
    }
  }
  
  # ============================================================================
  # HELPERS
  # ============================================================================
  
  resolve_coord_col <- function(meta_colnames, preferred, alternatives) {
    if (preferred %in% meta_colnames) return(preferred)
    for (alt in alternatives) {
      if (alt %in% meta_colnames) return(alt)
    }
    return(preferred)
  }
  
  # ============================================================================
  # PER-REGION ANNOTATION
  # ============================================================================
  
  annotate_region <- function(region_name, region_info) {
    
    cat(ts(), "\n========================================\n")
    cat(ts(), " REGION:", region_name, "\n")
    cat(ts(), "========================================\n\n")
    
    region_outdir <- file.path(region_info$dir, "annotation")
    dir.create(region_outdir, showWarnings = FALSE, recursive = TRUE)
    
    annotated_rds <- region_info$rds_path
    
    # --- Load ---
    if (!overwrite) {
      tmp <- readRDS(annotated_rds)
      
      if ("SCT" %in% names(tmp@assays)) {
        tmp@assays[["SCT"]] <- NULL
        tmp@active.assay <- "RNA"
      }
      
      if ("consensus_label" %in% colnames(tmp@meta.data)) {
        cat(ts(), "  Skipping -- already annotated:", annotated_rds, "\n")
        
        region_summary <- data.frame(
          region        = region_name,
          total_cells   = ncol(tmp),
          consensus     = sum(!is.na(tmp$consensus_label)),
          unlabeled     = sum(is.na(tmp$consensus_label)),
          pct_consensus = round(100 * sum(!is.na(tmp$consensus_label)) / ncol(tmp), 1),
          n_types       = length(unique(na.omit(tmp$consensus_label))),
          n_disagree_both_confident = if (annotation_method == "both") {
            sum(tmp$unlabeled_reason == "disagree_both_confident", na.rm = TRUE)
          } else { NA_integer_ },
          stringsAsFactors = FALSE
        )
        
        rm(tmp); gc(verbose = FALSE)
        return(region_summary)
      }
      xen <- tmp
      rm(tmp)
    } else {
      xen <- readRDS(annotated_rds)
      if ("SCT" %in% names(xen@assays)) {
        xen@assays[["SCT"]] <- NULL
        xen@active.assay <- "RNA"
      }
    }
    
    cat(ts(), "[LOAD] Loaded:", ncol(xen), "cells,", nrow(xen), "genes\n")
    cat(ts(), "[LOAD] Assays:", paste(names(xen@assays), collapse = ", "), "\n")
    cat(ts(), "[LOAD] Default assay:", xen@active.assay, "\n")
    
    # --- Normalize RNA ---
    
    available_layers <- Layers(xen[["RNA"]])
    if (!"data" %in% available_layers) {
      cat(ts(), "[NORM] No data layer found. Running NormalizeData...\n")
      xen <- NormalizeData(xen, assay = "RNA", verbose = FALSE)
    } else {
      rna_counts <- GetAssayData(xen, assay = "RNA", layer = "counts")
      rna_data   <- GetAssayData(xen, assay = "RNA", layer = "data")
      if (identical(rna_counts[1:min(100, nrow(rna_counts)), 1:min(10, ncol(rna_counts))],
                    rna_data[1:min(100, nrow(rna_counts)), 1:min(10, ncol(rna_counts))])) {
        cat(ts(), "[NORM] Data layer matches counts. Running NormalizeData...\n")
        xen <- NormalizeData(xen, assay = "RNA", verbose = FALSE)
      } else {
        cat(ts(), "[NORM] RNA data layer already log-normalized.\n")
      }
      rm(rna_counts, rna_data); gc(verbose = FALSE)
    }
    
    # ========================================================================
    # PRE-CLUSTERING for scalable SingleR
    #
    # High-resolution Louvain produces ~1000-3000 internally homogeneous
    # clusters. SingleR classifies cluster pseudobulk profiles instead of
    # individual cells, reducing runtime from days to minutes.
    # ========================================================================
    
    if (run_singler) {
      cat(ts(), "[CLUSTER] Pre-clustering for scalable annotation...\n")
      
      xen <- FindVariableFeatures(xen, assay = "RNA", verbose = FALSE)
      xen <- ScaleData(xen, assay = "RNA", verbose = FALSE)
      xen <- RunPCA(xen, assay = "RNA",
                    npcs = config$cluster_npcs, verbose = FALSE)
      
      # Drop scale.data immediately — it's the largest dense matrix
      xen[["RNA"]]$scale.data <- NULL
      gc(verbose = FALSE)
      
      xen <- FindNeighbors(xen, dims = config$cluster_dims, verbose = FALSE)
      xen <- FindClusters(xen, resolution = config$cluster_resolution,
                          verbose = FALSE)
      
      n_clusters <- length(unique(xen$seurat_clusters))
      cat(ts(), "[CLUSTER]", n_clusters, "clusters at resolution",
          config$cluster_resolution, "\n")
      
      # Warn if clustering is too coarse for accurate pseudobulk
      if (n_clusters < 100) {
        cat(ts(), "[WARN] Only", n_clusters,
            "clusters — consider increasing cluster_resolution for better",
            "SingleR accuracy.\n")
      }
    }
    
    # ========================================================================
    # SingleR — CLUSTER-LEVEL CLASSIFICATION
    # ========================================================================
    
    if (run_singler) {
      cat(ts(), "[SINGLER] Building SCE for classification...\n")
      
      counts_mat <- GetAssayData(xen, assay = "RNA", layer = "counts")
      data_mat   <- GetAssayData(xen, assay = "RNA", layer = "data")
      
      xen_sce <- SingleCellExperiment(
        assays = list(counts = counts_mat, logcounts = data_mat),
        colData = xen@meta.data
      )
      
      test_genes  <- rownames(xen_sce)
      train_genes <- rownames(singler_trained$ref)
      shared <- intersect(train_genes, test_genes)
      missing_in_test <- setdiff(train_genes, test_genes)
      
      cat(ts(), "[SINGLER] Trained on:", length(train_genes),
          "| This region:", length(test_genes),
          "| Shared:", length(shared), "\n")
      
      # ---- Memory-efficient gene alignment ----
      # v1 used rbind() + reindex which creates full intermediate copies.
      # This constructs the aligned matrix directly.
      
      lc_mat <- logcounts(xen_sce)
      ct_mat <- if ("counts" %in% assayNames(xen_sce)) counts(xen_sce) else NULL
      xen_coldata <- colData(xen_sce)
      n_test_cells <- ncol(xen_sce)
      
      if (length(missing_in_test) > 0) {
        cat(ts(), "[SINGLER] Padding", length(missing_in_test),
            "absent genes with zeros\n")
        
        # Build aligned logcounts directly
        lc_aligned <- Matrix(0, nrow = length(train_genes), ncol = n_test_cells,
                             sparse = TRUE,
                             dimnames = list(train_genes, colnames(xen_sce)))
        present <- train_genes[train_genes %in% test_genes]
        lc_aligned[present, ] <- lc_mat[present, ]
        lc_mat <- lc_aligned
        rm(lc_aligned)
        
        # Build aligned counts directly
        if (!is.null(ct_mat)) {
          ct_aligned <- Matrix(0, nrow = length(train_genes), ncol = n_test_cells,
                               sparse = TRUE,
                               dimnames = list(train_genes, colnames(xen_sce)))
          ct_aligned[present, ] <- ct_mat[present, ]
          ct_mat <- ct_aligned
          rm(ct_aligned)
        }
        rm(present)
      } else {
        lc_mat <- lc_mat[train_genes, ]
        if (!is.null(ct_mat)) {
          ct_mat <- ct_mat[train_genes, ]
        }
      }
      
      gc(verbose = FALSE)
      
      assay_list <- list(logcounts = lc_mat)
      if (!is.null(ct_mat)) assay_list$counts <- ct_mat
      
      xen_sce <- SingleCellExperiment(
        assays  = assay_list,
        colData = xen_coldata
      )
      rm(lc_mat, ct_mat, xen_coldata, assay_list)
      gc(verbose = FALSE)
      # ---- End gene alignment ----
      
      # ---- Cluster-level classification ----
      # Instead of classifying 220k individual cells, we pass the
      # `clusters` argument which tells SingleR to aggregate each
      # cluster into a pseudobulk profile and classify those.
      # Result: one row per cluster, mapped back to cells below.
      
      cluster_ids <- as.character(xen$seurat_clusters)
      names(cluster_ids) <- colnames(xen_sce)
      colData(xen_sce)$cluster <- cluster_ids
      
      n_to_classify <- length(unique(cluster_ids))
      cat(ts(), "[SINGLER] Classifying", n_to_classify,
          "cluster pseudobulk profiles (not", ncol(xen_sce), "individual cells)...\n")
      
      agg <- scuttle::aggregateAcrossCells(xen_sce, ids = colData(xen_sce)$cluster)
      logcounts(agg) <- log1p(counts(agg))
      
      singler_results <- classifySingleR(
        test      = agg,
        trained   = singler_trained,
        quantile  = config$singler_quantile,
        fine.tune = config$singler_fine_tune,
        BPPARAM   = MulticoreParam(config$rctd_max_cores)
      )
      
      rm(agg)
      
      cat(ts(), "[SINGLER] Classification complete.\n")
      
      # ---- Map cluster-level results back to individual cells ----
      # singler_results has one row per cluster (rownames = cluster IDs)
      cluster_names <- rownames(singler_results)
      
      cluster_labels    <- singler_results$labels
      names(cluster_labels) <- cluster_names
      
      cluster_pruned    <- singler_results$pruned.labels
      names(cluster_pruned) <- cluster_names
      
      cluster_delta     <- singler_results$delta.next
      names(cluster_delta) <- cluster_names
      
      cluster_max_score <- apply(singler_results$scores, 1, max)
      names(cluster_max_score) <- cluster_names
      
      # Map back to each cell via its cluster assignment
      xen$singler_label      <- as.character(cluster_labels[cluster_ids])
      xen$singler_pruned     <- as.character(cluster_pruned[cluster_ids])
      xen$singler_delta_next <- as.numeric(cluster_delta[cluster_ids])
      xen$singler_max_score  <- as.numeric(cluster_max_score[cluster_ids])
      
      # Per-cell delta percentile within each label (unchanged from v1)
      xen$singler_delta_pctile <- ave(
        xen$singler_delta_next, xen$singler_label,
        FUN = function(x) rank(x) / length(x)
      )
      
      n_pruned <- sum(is.na(xen$singler_pruned))
      n_types_found <- length(unique(na.omit(xen$singler_label)))
      cat(ts(), "[SINGLER]", n_pruned, "cells in pruned clusters.",
          n_types_found, "types assigned.\n")
      
      saveRDS(singler_results, file.path(region_outdir, "singler_raw_results.rds"))
      rm(xen_sce, singler_results, counts_mat, data_mat,
         cluster_ids, cluster_labels, cluster_pruned,
         cluster_delta, cluster_max_score)
      gc(verbose = FALSE)
    }
    
    # ========================================================================
    # RCTD
    # ========================================================================
    
    if (run_rctd) {
      cat(ts(), "[RCTD] Running", config$rctd_mode, "mode deconvolution",
          "with", config$rctd_max_cores, "cores...\n")
      
      common_genes <- intersect(rownames(xen), rownames(rctd_ref$counts))
      cat(ts(), "[RCTD] Gene intersection:", length(common_genes), "genes\n")
      
      spatial_counts <- GetAssayData(xen, assay = "RNA", layer = "counts")[common_genes, ]
      
      cx <- resolve_coord_col(colnames(xen@meta.data), config$coord_x,
                              c("x_centroid", "X", "centroid_x"))
      cy <- resolve_coord_col(colnames(xen@meta.data), config$coord_y,
                              c("y_centroid", "Y", "centroid_y"))
      
      if (!(cx %in% colnames(xen@meta.data))) stop("[ERROR] Cannot find x-coordinate column.")
      if (!(cy %in% colnames(xen@meta.data))) stop("[ERROR] Cannot find y-coordinate column.")
      
      coords <- data.frame(
        x = xen@meta.data[[cx]],
        y = xen@meta.data[[cy]],
        row.names = colnames(xen)
      )
      
      ref_rctd_obj <- Reference(
        counts     = rctd_ref$counts[common_genes, ],
        cell_types = rctd_ref$cell_types,
        nUMI       = rctd_ref$nUMI
      )
      
      query_rctd <- SpatialRNA(
        coords = coords,
        counts = spatial_counts,
        nUMI   = setNames(colSums(spatial_counts), colnames(spatial_counts))
      )
      
      rctd_obj <- create.RCTD(
        spatialRNA        = query_rctd,
        reference         = ref_rctd_obj,
        max_cores         = config$rctd_max_cores,
        gene_cutoff       = config$rctd_gene_cutoff,
        fc_cutoff         = config$rctd_fc_cutoff,
        fc_cutoff_reg     = config$rctd_fc_cutoff_reg,
        UMI_min           = config$rctd_UMI_min,
        counts_MIN        = config$rctd_counts_MIN,
        CELL_MIN_INSTANCE = config$rctd_CELL_MIN_INSTANCE
      )
      
      rctd_obj <- run.RCTD(rctd_obj, doublet_mode = config$rctd_mode)
      
      xen$rctd_spot_class   <- NA_character_
      xen$rctd_first_type   <- NA_character_
      xen$rctd_second_type  <- NA_character_
      xen$rctd_first_weight <- NA_real_
      
      if (config$rctd_mode == "full") {
        rctd_weights <- rctd_obj@results$weights
        rctd_cells   <- rownames(rctd_weights)
        matching     <- intersect(colnames(xen), rctd_cells)
        match_idx    <- match(matching, colnames(xen))
        
        row_sums <- rowSums(rctd_weights[matching, , drop = FALSE])
        norm_weights <- rctd_weights[matching, , drop = FALSE] / row_sums
        
        top_type   <- colnames(norm_weights)[apply(norm_weights, 1, which.max)]
        top_weight <- apply(norm_weights, 1, max)
        
        xen$rctd_first_type[match_idx]   <- top_type
        xen$rctd_first_weight[match_idx] <- top_weight
        xen$rctd_spot_class[match_idx]   <- "full"
        
        rm(rctd_weights, norm_weights)
        
      } else if (config$rctd_mode == "doublet") {
        rctd_df      <- rctd_obj@results$results_df
        rctd_weights <- rctd_obj@results$weights
        rctd_cells   <- rownames(rctd_df)
        matching     <- intersect(colnames(xen), rctd_cells)
        match_idx    <- match(matching, colnames(xen))
        
        xen$rctd_first_type[match_idx]   <- as.character(rctd_df[matching, "first_type"])
        xen$rctd_spot_class[match_idx]   <- as.character(rctd_df[matching, "spot_class"])
        xen$rctd_second_type[match_idx]  <- as.character(rctd_df[matching, "second_type"])
        
        if (!is.null(rctd_weights) && nrow(rctd_weights) > 0) {
          weight_cells <- intersect(matching, rownames(rctd_weights))
          if (length(weight_cells) > 0) {
            first_types <- as.character(rctd_df[weight_cells, "first_type"])
            valid <- first_types %in% colnames(rctd_weights)
            if (any(valid)) {
              vc <- weight_cells[valid]
              ft <- first_types[valid]
              xen$rctd_first_weight[match(vc, colnames(xen))] <- rctd_weights[cbind(vc, ft)]
            }
          }
        }
        rm(rctd_df, rctd_weights)
        
      } else {
        # Multi mode
        rctd_weights <- rctd_obj@results$weights
        rctd_cells   <- rownames(rctd_weights)
        matching     <- intersect(colnames(xen), rctd_cells)
        match_idx    <- match(matching, colnames(xen))
        
        row_sums <- rowSums(rctd_weights[matching, , drop = FALSE])
        norm_weights <- rctd_weights[matching, , drop = FALSE] / row_sums
        
        top_type   <- colnames(norm_weights)[apply(norm_weights, 1, which.max)]
        top_weight <- apply(norm_weights, 1, max)
        
        xen$rctd_first_type[match_idx]   <- top_type
        xen$rctd_first_weight[match_idx] <- top_weight
        xen$rctd_spot_class[match_idx]   <- "multi"
        
        rm(rctd_weights, norm_weights)
      }
      
      n_rctd_dropped <- ncol(xen) - length(matching)
      if (n_rctd_dropped > 0) {
        cat(ts(), "[RCTD] Dropped", n_rctd_dropped, "cells (below UMI/gene thresholds)\n")
      }
      
      cat(ts(), "[RCTD] Complete.", length(matching), "cells assigned.\n")
      cat(ts(), "[RCTD] Mode:", config$rctd_mode, "\n")
      
      rm(rctd_obj, ref_rctd_obj, query_rctd, spatial_counts)
      gc(verbose = FALSE)
    }
    
    # ========================================================================
    # CONSENSUS / LABELING
    # ========================================================================
    
    if (annotation_method == "both") {
      # --- Dual-method consensus ---
      cat(ts(), "[CONSENSUS] Applying dual-method confidence filter...\n")
      
      singler_pass <- (
        !is.na(xen$singler_pruned) &
          !is.na(xen$singler_delta_next) &
          xen$singler_delta_next >= config$singler_delta_floor &
          xen$singler_delta_pctile >= (1 - config$singler_delta_pctile_threshold)
      )
      singler_pass[is.na(singler_pass)] <- FALSE
      
      # For full mode, rctd_spot_class is always "full" so skip class check
      if (config$rctd_mode == "doublet") {
        rctd_pass <- (
          xen$rctd_spot_class %in% config$rctd_confident_classes &
            !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      } else {
        rctd_pass <- (
          !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      }
      rctd_pass[is.na(rctd_pass)] <- FALSE
      
      labels_agree <- (!is.na(xen$singler_label) &
                         !is.na(xen$rctd_first_type) &
                         xen$singler_label == xen$rctd_first_type)
      labels_agree[is.na(labels_agree)] <- FALSE
      
      consensus_pass <- labels_agree & singler_pass & rctd_pass
      xen$consensus_label <- ifelse(consensus_pass, xen$singler_label, NA_character_)
      
      idx <- !consensus_pass
      xen$unlabeled_reason <- NA_character_
      if (any(idx)) {
        xen$unlabeled_reason[idx] <- case_when(
          is.na(xen$rctd_first_type[idx])
          ~ "rctd_dropped",
          labels_agree[idx] & !singler_pass[idx] & !rctd_pass[idx]
          ~ "agree_both_low_confidence",
          labels_agree[idx] & !singler_pass[idx] & rctd_pass[idx]
          ~ "agree_singler_low_confidence",
          labels_agree[idx] & singler_pass[idx] & !rctd_pass[idx]
          ~ "agree_rctd_low_confidence",
          !labels_agree[idx] & singler_pass[idx] & rctd_pass[idx]
          ~ "disagree_both_confident",
          !labels_agree[idx] & !singler_pass[idx] & !rctd_pass[idx]
          ~ "disagree_both_low_confidence",
          !labels_agree[idx] & singler_pass[idx] & !rctd_pass[idx]
          ~ "disagree_only_singler_confident",
          !labels_agree[idx] & !singler_pass[idx] & rctd_pass[idx]
          ~ "disagree_only_rctd_confident",
          TRUE ~ "other"
        )
      }
      
    } else if (annotation_method == "singler") {
      # --- SingleR only ---
      cat(ts(), "[LABEL] Applying SingleR confidence filter...\n")
      
      singler_pass <- (
        !is.na(xen$singler_pruned) &
          !is.na(xen$singler_delta_next) &
          xen$singler_delta_next >= config$singler_delta_floor &
          xen$singler_delta_pctile >= (1 - config$singler_delta_pctile_threshold)
      )
      singler_pass[is.na(singler_pass)] <- FALSE
      
      xen$consensus_label <- ifelse(singler_pass, xen$singler_label, NA_character_)
      
      xen$unlabeled_reason <- NA_character_
      idx <- !singler_pass
      if (any(idx)) {
        xen$unlabeled_reason[idx] <- case_when(
          is.na(xen$singler_pruned[idx])
          ~ "singler_pruned",
          !is.na(xen$singler_delta_next[idx]) &
            xen$singler_delta_next[idx] < config$singler_delta_floor
          ~ "singler_low_delta",
          !is.na(xen$singler_delta_pctile[idx]) &
            xen$singler_delta_pctile[idx] < (1 - config$singler_delta_pctile_threshold)
          ~ "singler_low_delta_pctile",
          TRUE ~ "other"
        )
      }
      
    } else {
      # --- RCTD only ---
      cat(ts(), "[LABEL] Applying RCTD confidence filter...\n")
      
      if (config$rctd_mode == "doublet") {
        rctd_pass <- (
          xen$rctd_spot_class %in% config$rctd_confident_classes &
            !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      } else {
        rctd_pass <- (
          !is.na(xen$rctd_first_weight) &
            xen$rctd_first_weight >= config$rctd_weight_threshold
        )
      }
      rctd_pass[is.na(rctd_pass)] <- FALSE
      
      xen$consensus_label <- ifelse(rctd_pass, xen$rctd_first_type, NA_character_)
      
      xen$unlabeled_reason <- NA_character_
      idx <- !rctd_pass
      if (any(idx)) {
        xen$unlabeled_reason[idx] <- case_when(
          is.na(xen$rctd_first_type[idx])
          ~ "rctd_dropped",
          !is.na(xen$rctd_first_weight[idx]) &
            xen$rctd_first_weight[idx] < config$rctd_weight_threshold
          ~ "rctd_low_weight",
          config$rctd_mode == "doublet" &
            !(xen$rctd_spot_class[idx] %in% config$rctd_confident_classes)
          ~ "rctd_uncertain_spot_class",
          TRUE ~ "other"
        )
      }
    }
    
    n_total     <- ncol(xen)
    n_consensus <- sum(!is.na(xen$consensus_label))
    pct         <- round(100 * n_consensus / n_total, 1)
    cat(ts(), "[LABEL]", region_name, ":", n_consensus, "/", n_total,
        "(", pct, "%) labeled\n")
    
    reason_tbl <- table(xen$unlabeled_reason, useNA = "ifany")
    cat(ts(), "  Reason breakdown:\n")
    for (r in names(reason_tbl)) {
      cat(ts(), "    ", r, ":", reason_tbl[r], "\n")
    }
    
    consensus_tbl <- sort(table(xen$consensus_label), decreasing = TRUE)
    cat(ts(), "  Label composition:\n")
    for (ct in names(consensus_tbl)) {
      cat(ts(), "    ", ct, ":", consensus_tbl[ct], "\n")
    }
    
    # --- Exports ---
    cat(ts(), "[EXPORT] Writing tables...\n")
    
    cx <- resolve_coord_col(colnames(xen@meta.data), config$coord_x,
                            c("x_centroid", "X", "centroid_x"))
    cy <- resolve_coord_col(colnames(xen@meta.data), config$coord_y,
                            c("y_centroid", "Y", "centroid_y"))
    
    export_cols <- intersect(
      c("consensus_label", "unlabeled_reason",
        "singler_label", "singler_pruned", "singler_delta_next",
        "singler_delta_pctile", "singler_max_score",
        "rctd_spot_class", "rctd_first_type", "rctd_second_type",
        "rctd_first_weight",
        "nCount_RNA", "nFeature_RNA", cx, cy),
      colnames(xen@meta.data)
    )
    
    fwrite(xen@meta.data[, export_cols, drop = FALSE],
           file.path(region_outdir, "all_cells_annotation.csv.gz"))
    
    unlabeled_mask <- is.na(xen$consensus_label)
    if (any(unlabeled_mask)) {
      fwrite(xen@meta.data[unlabeled_mask, export_cols, drop = FALSE],
             file.path(region_outdir, "unlabeled_cells.csv.gz"))
    }
    
    # Confusion matrix only if both methods ran
    if (annotation_method == "both") {
      labels_agree <- (!is.na(xen$singler_label) &
                         !is.na(xen$rctd_first_type) &
                         xen$singler_label == xen$rctd_first_type)
      labels_agree[is.na(labels_agree)] <- FALSE
      
      confusion <- table(SingleR = xen$singler_label,
                         RCTD = xen$rctd_first_type, useNA = "ifany")
      fwrite(as.data.frame.matrix(confusion),
             file.path(region_outdir, "singler_vs_rctd_confusion.csv"),
             row.names = TRUE)
      
      disagree_mask <- !labels_agree & !is.na(xen$singler_label) &
        !is.na(xen$rctd_first_type)
      if (any(disagree_mask)) {
        dtab <- table(SingleR = xen$singler_label[disagree_mask],
                      RCTD = xen$rctd_first_type[disagree_mask])
        fwrite(as.data.frame.matrix(dtab),
               file.path(region_outdir, "disagreement_matrix.csv"),
               row.names = TRUE)
      }
    }
    
    # --- Diagnostic plots ---
    cat(ts(), "[DIAG] Generating plots...\n")
    
    has_coords <- cx %in% colnames(xen@meta.data) && cy %in% colnames(xen@meta.data)
    
    if (has_coords) {
      plot_df <- data.frame(
        x = xen@meta.data[[cx]], y = xen@meta.data[[cy]],
        label = ifelse(is.na(xen$consensus_label), "Unlabeled", xen$consensus_label),
        reason = xen$unlabeled_reason, stringsAsFactors = FALSE
      )
      
      p1 <- ggplot(plot_df, aes(x, y, color = label)) +
        geom_point(size = 0.1, alpha = 0.5) +
        coord_fixed() + theme_minimal(base_size = 10) +
        ggtitle(paste(region_name, "--", annotation_method, "labels")) +
        theme(legend.key.size = unit(0.3, "cm"))
      ggsave(file.path(region_outdir, "spatial_consensus.pdf"),
             p1, width = 14, height = 10)
      
      unlabeled_df <- plot_df[plot_df$label == "Unlabeled", ]
      if (nrow(unlabeled_df) > 0) {
        p2 <- ggplot(unlabeled_df, aes(x, y, color = reason)) +
          geom_point(size = 0.2, alpha = 0.6) +
          coord_fixed() + theme_minimal(base_size = 10) +
          ggtitle(paste(region_name, "-- Unlabeled: reason codes"))
        ggsave(file.path(region_outdir, "spatial_unlabeled_reasons.pdf"),
               p2, width = 14, height = 10)
      }
    }
    
    # SingleR delta plot (only if SingleR ran)
    if (run_singler) {
      p3 <- ggplot(xen@meta.data, aes(x = singler_label, y = singler_delta_next)) +
        geom_violin(fill = "lightblue", alpha = 0.4, scale = "width") +
        geom_jitter(aes(color = ifelse(is.na(consensus_label), "Unlabeled", "Labeled")),
                    size = 0.05, alpha = 0.2, width = 0.2) +
        scale_color_manual(values = c(Labeled = "grey30", Unlabeled = "red")) +
        geom_hline(yintercept = config$singler_delta_floor, lty = 2, color = "darkred") +
        coord_flip() + theme_minimal(base_size = 9) +
        ggtitle(paste(region_name, "-- SingleR delta.next by label"))
      ggsave(file.path(region_outdir, "singler_delta_distribution.pdf"),
             p3, width = 10, height = max(6, 0.4 * length(unique(xen$singler_label))))
    }
    
    # Confidence scatter (only if both ran)
    if (annotation_method == "both") {
      p4 <- ggplot(xen@meta.data, aes(x = singler_delta_next, y = rctd_first_weight)) +
        geom_point(aes(color = ifelse(is.na(consensus_label), "Unlabeled", "Labeled")),
                   size = 0.15, alpha = 0.2) +
        scale_color_manual(values = c(Labeled = "steelblue", Unlabeled = "red")) +
        geom_hline(yintercept = config$rctd_weight_threshold, lty = 2, color = "grey40") +
        geom_vline(xintercept = config$singler_delta_floor, lty = 2, color = "grey40") +
        annotate("rect",
                 xmin = config$singler_delta_floor, xmax = Inf,
                 ymin = config$rctd_weight_threshold, ymax = Inf,
                 fill = "steelblue", alpha = 0.05) +
        theme_minimal(base_size = 10) +
        ggtitle(paste(region_name, "-- Confidence space"))
      ggsave(file.path(region_outdir, "confidence_scatter.pdf"),
             p4, width = 8, height = 7)
    }
    
    # --- Save annotated object ---
    saveRDS(xen, annotated_rds)
    cat(ts(), "[DONE]", region_name, "->", annotated_rds, "\n")
    
    # --- Build summary, then free object ---
    region_summary <- data.frame(
      region        = region_name,
      total_cells   = ncol(xen),
      consensus     = sum(!is.na(xen$consensus_label)),
      unlabeled     = sum(is.na(xen$consensus_label)),
      pct_consensus = pct,
      n_types       = length(unique(na.omit(xen$consensus_label))),
      n_disagree_both_confident = if (annotation_method == "both") {
        sum(xen$unlabeled_reason == "disagree_both_confident", na.rm = TRUE)
      } else { NA_integer_ },
      stringsAsFactors = FALSE
    )
    
    rm(xen); gc(verbose = FALSE)
    cat("\n")
    
    return(region_summary)
  }
  
  # ============================================================================
  # RUN ALL REGIONS
  # ============================================================================
  
  cat(ts(), "================================================================\n")
  cat(ts(), " Annotating", length(regions), "Xenium regions\n")
  cat(ts(), "================================================================\n")
  
  summaries <- mapply(
    annotate_region,
    region_name = names(regions),
    region_info = regions,
    SIMPLIFY = FALSE
  )
  
  summaries <- Filter(Negate(is.null), summaries)
  
  # ============================================================================
  # GLOBAL SUMMARY
  # ============================================================================
  
  if (length(summaries) > 0) {
    cat(ts(), "================================================================\n")
    cat(ts(), " GLOBAL SUMMARY\n")
    cat(ts(), "================================================================\n\n")
    
    summary_df <- do.call(rbind, summaries)
    fwrite(summary_df, file.path(proseg_dir, "annotation_summary.csv"))
    
    for (i in seq_len(nrow(summary_df))) {
      r <- summary_df[i, ]
      cat(ts(), " ", r$region, ":", r$consensus, "/", r$total_cells,
          "(", r$pct_consensus, "%) labeled,",
          r$n_types, "types")
      if (annotation_method == "both") {
        cat(",", r$n_disagree_both_confident, "high-confidence disagreements")
      }
      cat("\n")
    }
  }
  
  cat(ts(), "\nAnnotation complete. Method:", annotation_method, "\n")
  invisible(summary_df)
}

#######################################################################
#######################################################################
#                                                                     #
#                       SEGMENTATION DIAGNOSTICS                      #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: polygon plots, segmentation comparisons, nuclei per cell

#PLOT SEGMENTATION POLYGONS
plot_polygons <- function(seu, polygon_slot, fill_by = NULL, line_col = "grey30",
                          line_width = 0.1) {
  require(ggplot2)
  require(sf)
  
  poly_sf <- seu@misc[[polygon_slot]]
  if (is.null(poly_sf)) stop("No polygons in @misc$", polygon_slot)
  
  if (!is.null(fill_by)) {
    if (!fill_by %in% colnames(seu@meta.data)) {
      stop("Metadata column '", fill_by, "' not found.")
    }
    meta_df <- data.frame(
      cell = rownames(seu@meta.data),
      fill_var = seu@meta.data[[fill_by]]
    )
    poly_sf <- merge(poly_sf, meta_df, by = "cell", all.x = TRUE, sort = FALSE)
    poly_sf$fill_var <- as.factor(poly_sf$fill_var)
    
    p <- ggplot(poly_sf) +
      geom_sf(aes(fill = fill_var), color = line_col, linewidth = line_width) +
      theme_void() + coord_sf() +
      labs(fill = fill_by)
  } else {
    p <- ggplot(poly_sf) +
      geom_sf(fill = NA, color = line_col, linewidth = line_width) +
      theme_void() + coord_sf()
  }
  
  return(p)
}

#PLOT A FEATURE ON SEGMENTATION POLYGONS
plot_spatial_feature <- function(seu, polygon_slot, feature,
                                 palette = viridis::viridis(100),
                                 line_col = "grey50", line_width = 0.05) {
  require(ggplot2)
  require(sf)
  
  if (!feature %in% rownames(seu)) stop("Feature '", feature, "' not found.")
  
  poly_sf <- seu@misc[[polygon_slot]]
  if (is.null(poly_sf)) stop("No polygons in @misc$", polygon_slot)
  
  expr <- FetchData(seu, vars = feature)
  expr$cell <- rownames(expr)
  poly_sf <- merge(poly_sf, expr, by = "cell", all.x = TRUE, sort = FALSE)
  
  ggplot(poly_sf) +
    geom_sf(aes(fill = .data[[feature]]), color = line_col, linewidth = line_width) +
    scale_fill_gradientn(colors = palette, na.value = "#eeeeee", name = feature) +
    theme_void() + coord_sf()
}

#SIDE BY SIDE SEGMENTATION COMPARISON
compare_segmentations <- function(seu,
                                  slot1 = "proseg_polygons",
                                  slot2 = "xenium_polygons",
                                  fill_by = NULL,
                                  labels = c("ProSeg", "Xenium")) {
  require(patchwork)
  
  p1 <- plot_polygons(seu, slot1, fill_by = fill_by) + ggtitle(labels[1])
  p2 <- plot_polygons(seu, slot2, fill_by = fill_by) + ggtitle(labels[2])
  
  p1 + p2
}

#ZOOMED SEGMENTATION COMPARISON
compare_segmentations_zoomed <- function(seu,
                                         slot1 = "proseg_polygons",
                                         slot2 = "xenium_polygons",
                                         x_range, y_range,
                                         fill_by = NULL,
                                         labels = c("ProSeg", "Xenium")) {
  require(patchwork)
  require(sf)
  
  crop_and_plot <- function(slot, label) {
    poly_sf <- seu@misc[[slot]]
    poly_sf <- st_crop(poly_sf,
                       xmin = x_range[1], xmax = x_range[2],
                       ymin = y_range[1], ymax = y_range[2])
    
    if (!is.null(fill_by)) {
      meta_df <- data.frame(
        cell = rownames(seu@meta.data),
        fill_var = as.factor(seu@meta.data[[fill_by]])
      )
      poly_sf <- merge(poly_sf, meta_df, by = "cell", all.x = TRUE, sort = FALSE)
      ggplot(poly_sf) +
        geom_sf(aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
        theme_void() + coord_sf() + ggtitle(label) + labs(fill = fill_by)
    } else {
      ggplot(poly_sf) +
        geom_sf(fill = NA, color = "grey30", linewidth = 0.2) +
        theme_void() + coord_sf() + ggtitle(label)
    }
  }
  
  p1 <- crop_and_plot(slot1, labels[1])
  p2 <- crop_and_plot(slot2, labels[2])
  p1 + p2
}

#TWO PANEL SEGMENTATION COMPARISON
plot_two_segmentation_panels <- function(seu,
                                         cell_seg_slot1,
                                         cell_seg_slot2,
                                         nuclei_slot,
                                         labels = c("ProSeg", "Xenium"),
                                         nuclei_fill = "#4a90e2",
                                         nuclei_border = "#d9d9d9",
                                         cell_border_color = "black",
                                         cell_fill1 = "#fdbf6f",
                                         cell_fill2 = "#b2df8a",
                                         zoom_bbox = NULL) {
  require(ggplot2)
  require(sf)
  require(patchwork)
  
  cells1 <- seu@misc[[cell_seg_slot1]]
  cells2 <- seu@misc[[cell_seg_slot2]]
  nuclei <- seu@misc[[nuclei_slot]]
  
  if (!is.null(zoom_bbox)) {
    cells1 <- st_crop(cells1, xmin = zoom_bbox[1], xmax = zoom_bbox[2],
                      ymin = zoom_bbox[3], ymax = zoom_bbox[4])
    cells2 <- st_crop(cells2, xmin = zoom_bbox[1], xmax = zoom_bbox[2],
                      ymin = zoom_bbox[3], ymax = zoom_bbox[4])
    nuclei <- st_crop(nuclei, xmin = zoom_bbox[1], xmax = zoom_bbox[2],
                      ymin = zoom_bbox[3], ymax = zoom_bbox[4])
  }
  
  plot_panel <- function(cells, cell_fill, label) {
    ggplot() +
      geom_sf(data = cells, fill = cell_fill,
              color = cell_border_color, alpha = 0.4) +
      geom_sf(data = nuclei, fill = nuclei_fill,
              color = nuclei_border, alpha = 0.6) +
      theme_minimal() + ggtitle(label) + coord_sf(expand = FALSE)
  }
  
  p1 <- plot_panel(cells1, cell_fill1, labels[1])
  p2 <- plot_panel(cells2, cell_fill2, labels[2])
  p1 + p2
}

#COUNT NUCLEI PER SEGMENTED CELL
compute_nuclei_per_cell <- function(seu,
                                    cell_slot,
                                    nuclei_slot) {
  require(sf)
  
  cells_sf  <- seu@misc[[cell_slot]]
  nuclei_sf <- seu@misc[[nuclei_slot]]
  
  if (is.null(cells_sf))  stop("No polygons in @misc$", cell_slot)
  if (is.null(nuclei_sf)) stop("No polygons in @misc$", nuclei_slot)
  
  cells_sf  <- st_make_valid(cells_sf)
  nuclei_sf <- st_make_valid(nuclei_sf)
  
  message("Computing nuclei per cell...")
  overlaps <- st_within(nuclei_sf, cells_sf, sparse = TRUE)
  
  cell_hits <- sapply(overlaps, function(hit) {
    if (length(hit) == 0) return(NA_integer_)
    hit[1]
  })
  
  counts <- table(na.omit(cell_hits))
  data.frame(
    cell     = cells_sf$cell[as.integer(names(counts))],
    n_nuclei = as.integer(counts)
  )
}

#COMPARE NUCLEI ASSIGNMENT BETWEEN SEGMENTATIONS
compare_nuclei_segmentation <- function(seu,
                                        cell_slot1, cell_slot2,
                                        nuclei_slot,
                                        labels = c("Segmentation 1", "Segmentation 2")) {
  require(patchwork)
  require(ggplot2)
  
  counts1 <- compute_nuclei_per_cell(seu, cell_slot1, nuclei_slot)
  counts2 <- compute_nuclei_per_cell(seu, cell_slot2, nuclei_slot)
  
  plot_dist <- function(counts_df, title) {
    ggplot(counts_df, aes(x = n_nuclei)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
      theme_minimal() + ggtitle(title) +
      labs(x = "Nuclei per cell", y = "Count")
  }
  
  p1 <- plot_dist(counts1, labels[1])
  p2 <- plot_dist(counts2, labels[2])
  p1 + p2
}

#######################################################################
#######################################################################
#                                                                     #
#                          MARKER COMPARISON                          #
#                                                                     #
#######################################################################
#######################################################################

#CONTENTS: top markers per label, jaccard overlap, heatmaps

#TOP MARKERS PER LABEL
get_top_markers <- function(obj, label_col, n_top = 200) {
  
  Idents(obj) <- label_col
  
  if (requireNamespace("presto", quietly = TRUE)) {
    cat("  Using presto for marker detection...\n")
    
    expr_mat <- GetAssayData(obj, layer = "data")
    labels   <- obj@meta.data[[label_col]]
    
    res <- presto::wilcoxauc(expr_mat, labels)
    res <- as.data.table(res)
    
    res <- res[logFC > 0 & padj < 0.05]
    res <- res[order(group, -auc)]
    markers <- res[, head(.SD, n_top), by = group]
    
    marker_list <- split(markers$feature, markers$group)
    
  } else {
    cat("  presto not available. Using FindAllMarkers (this may be slow)...\n")
    
    all_markers <- FindAllMarkers(obj,
                                  only.pos        = TRUE,
                                  min.pct         = 0.1,
                                  logfc.threshold = 0.25,
                                  test.use        = "wilcox",
                                  verbose         = FALSE)
    all_markers <- as.data.table(all_markers)
    all_markers <- all_markers[order(cluster, -avg_log2FC)]
    top <- all_markers[, head(.SD, n_top), by = cluster]
    
    marker_list <- split(top$gene, top$cluster)
  }
  
  return(marker_list)
}

#JACCARD OVERLAP BETWEEN TWO MARKER SETS
compute_jaccard <- function(markers_row, markers_col) {
  
  row_names <- names(markers_row)
  col_names <- names(markers_col)
  
  mat <- matrix(0,
                nrow = length(row_names),
                ncol = length(col_names),
                dimnames = list(row_names, col_names))
  
  for (i in seq_along(row_names)) {
    set_a <- markers_row[[i]]
    for (j in seq_along(col_names)) {
      set_b <- markers_col[[j]]
      
      n_intersect <- length(intersect(set_a, set_b))
      n_union     <- length(union(set_a, set_b))
      
      mat[i, j] <- if (n_union == 0) 0 else n_intersect / n_union
    }
  }
  
  return(mat)
}

#JACCARD HEATMAP
plot_jaccard_heatmap <- function(jaccard_mat, nm_row, nm_col,
                                 pair_tag, output_dir, n_top_markers,
                                 filtered = TRUE) {
  
  if (filtered) {
    max_per_col <- apply(jaccard_mat, 2, max)
    relevant_cols <- names(max_per_col[max_per_col >= 0.02])
    if (length(relevant_cols) < 3) {
      relevant_cols <- names(sort(max_per_col, decreasing = TRUE))[
        1:min(30, length(max_per_col))]
    }
    plot_mat <- jaccard_mat[, relevant_cols, drop = FALSE]
    suffix   <- ""
  } else {
    plot_mat <- jaccard_mat
    suffix   <- "_full"
  }
  
  fname <- file.path(output_dir,
                     paste0("label_jaccard_heatmap", suffix, "_", pair_tag, ".pdf"))
  
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    library(pheatmap)
    
    pdf(fname,
        width  = max(14, ncol(plot_mat) * 0.25),
        height = max(8, nrow(plot_mat) * 0.5))
    
    pheatmap(plot_mat,
             color           = colorRampPalette(c("white", "cornflowerblue",
                                                  "darkblue"))(100),
             cluster_rows    = TRUE,
             cluster_cols    = TRUE,
             clustering_method = "ward.D2",
             display_numbers = filtered,
             number_format   = "%.2f",
             number_color    = if (filtered) ifelse(plot_mat > 0.15, "white", "grey30") else "grey30",
             fontsize_number = 7,
             fontsize_row    = 10,
             fontsize_col    = if (filtered) 8 else 6,
             angle_col       = 45,
             main            = paste0("Jaccard: ", nm_row, " (rows) vs ",
                                      nm_col, " (cols)",
                                      if (filtered) paste0("\nTop ", n_top_markers, " markers") else ""),
             border_color    = NA)
    dev.off()
    
  } else {
    plot_df <- as.data.frame(as.table(plot_mat))
    colnames(plot_df) <- c("type_row", "type_col", "Jaccard")
    
    p <- ggplot(plot_df, aes(x = type_col, y = type_row, fill = Jaccard)) +
      geom_tile(color = "grey90") +
      geom_text(aes(label = ifelse(Jaccard >= 0.05,
                                   sprintf("%.2f", Jaccard), "")),
                size = 2.5) +
      scale_fill_gradient(low = "white", high = "darkblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y = element_text(size = 10)) +
      labs(title = paste0("Jaccard: ", nm_row, " vs ", nm_col),
           x = nm_col, y = nm_row)
    
    ggsave(fname, p,
           width = max(14, ncol(plot_mat) * 0.3), height = 10)
  }
  
  cat("  Saved:", basename(fname), "\n")
}
