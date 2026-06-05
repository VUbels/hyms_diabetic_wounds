################################################################################
#
#  Two-tier manual annotation of merged Xenium (proseg) data
#
#  Tier 1 — broad types:  fibroblast, keratinocyte, immune, endothelial, etc.
#  Tier 2 — subtypes:     e.g. FIB_papillary, KC_basal, KC_spinous, …
#
#  Workflow:
#    1. Load & merge 6 proseg RDS objects
#    2. SCTransform → PCA → UMAP → clustering (low resolution)
#    3. Inspect markers → assign broad labels  (TIER 1 — you edit the mapping)
#    4. Sub-cluster each broad type at higher resolution
#    5. Inspect sub-markers → assign fine labels (TIER 2 — you edit the mapping)
#    6. Combine everything back into a single annotated object and save
#
################################################################################

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

set.seed(42)

# ── Paths ────────────────────────────────────────────────────────────────────
proseg_dir <- "./proseg_results"
output_dir <- "./manual_annotation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Discover and load the 6 RDS files ────────────────────────────────────────

rds_files <- list.files(
  proseg_dir,
  pattern    = "_proseg_seurat\\.rds$",
  recursive  = TRUE,
  full.names = TRUE
)
cat("Found", length(rds_files), "RDS files:\n")
print(basename(rds_files))

obj_list    <- list()
sample_ids  <- c()

for (f in rds_files) {
  nm <- sub("_proseg_seurat\\.rds$", "", basename(f))
  cat("Loading:", nm, "\n")
  obj <- readRDS(f)
  DefaultAssay(obj) <- "RNA"
  
  # Proseg produces fractional expected counts.
  # SCTransform's NB model expects integers — round them.
  obj[["RNA"]]$counts <- round(obj[["RNA"]]$counts)
  
  obj$sample <- nm
  obj_list[[nm]]  <- obj
  sample_ids      <- c(sample_ids, nm)
  rm(obj); gc(verbose = FALSE)
}

# ── Merge (not integrate) ────────────────────────────────────────────────────
# Same Xenium panel across all samples → minimal technical batch.
# SCTransform handles per-sample depth via split layers automatically.
# Integration (Harmony) would remove biological L-vs-D differences.

merged <- merge(
  x           = obj_list[[1]],
  y           = obj_list[-1],
  add.cell.ids = sample_ids
)
rm(obj_list); gc()

cat("Merged object:", ncol(merged), "cells,", nrow(merged), "genes\n")
cat("Samples:", paste(unique(merged$sample), collapse = ", "), "\n")


################################################################################
# STEP 1 — SCTransform, PCA, UMAP, clustering
################################################################################

merged <- SCTransform(merged, assay = "RNA", verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, features = rownames(merged), verbose = FALSE)

# Quick elbow plot — adjust dims if needed
ElbowPlot(merged, ndims = 30)
ggsave(file.path(output_dir, "elbow_plot.pdf"), width = 6, height = 4)

merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.8, verbose = FALSE, algorithm = 4)

# ── Visual batch check ───────────────────────────────────────────────────────
# If the UMAP splits by sample rather than by biology, consider Harmony.
p1 <- DimPlot(merged, group.by = "sample", raster = FALSE) + ggtitle("By sample")
p2 <- DimPlot(merged, label = TRUE, raster = FALSE)       + ggtitle("Clusters")
p1 + p2
ggsave(file.path(output_dir, "01_umap_batch_check.pdf"), width = 14, height = 6)


# Visualise one FOV with the merged clusters

ImageDimPlot(
  merged,
  fov          = "proseg",
  group.by     = "seurat_clusters",
  boundaries   = "segmentation",
  border.color = NA,
  border.size  = 0,
  dark.background = TRUE
)

cat("\n=== CHECK 01_umap_batch_check.pdf ===\n")
cat("If samples mix well → proceed as-is.\n")
cat("If they separate by sample → uncomment the Harmony block below and re-run.\n\n")

# ── Optional: Harmony integration (uncomment if needed) ──────────────────────
# library(harmony)
# merged <- RunHarmony(merged, group.by.vars = "sample", reduction = "pca",
#                      assay.use = "SCT", reduction.save = "harmony")
# merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
# merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
# merged <- FindClusters(merged, resolution = 0.3, verbose = FALSE)


################################################################################
# STEP 2 — Explore markers for TIER 1 (broad annotation)
################################################################################

# Find markers for every cluster
all_markers <- FindAllMarkers(
  merged,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.3,
  verbose         = FALSE
)

write.csv(all_markers, file.path(output_dir, "02_all_cluster_markers.csv"),
          row.names = FALSE)
cat("Cluster markers saved to 02_all_cluster_markers.csv\n")

# Top 5 per cluster for quick inspection
top5 <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(avg_log2FC, n = 5)
print(top5, n = 60)

# ── Dot plots of known broad-type markers ────────────────────────────────────
# Adjust these gene lists to match your panel.
broad_markers <- list(
  Keratinocyte = c("KRT14", "KRT5", "KRT10", "KRT1", "LORICRIN", "IVL"),
  Fibroblast   = c("COL1A1", "COL1A2", "DCN", "LUM", "VIM", "PDGFRA"),
  Immune       = c("PTPRC", "CD3E", "CD8A", "CD68", "CD163", "CD14"),
  Endothelial  = c("PECAM1", "VWF", "CDH5", "FLT1"),
  Melanocyte   = c("MLANA", "PMEL", "TYRP1", "DCT"),
  Smooth_Musc  = c("ACTA2", "MYH11", "TAGLN", "DES"),
  Neural_Crest = c("SOX10", "S100B", "PLP1", "MPZ")
)

FeaturePlot(merged, features = "Pecam1")

# Keep only genes that exist in the panel
broad_markers_present <- lapply(broad_markers, function(g) {
  g[g %in% rownames(merged)]
})
broad_markers_present <- broad_markers_present[lengths(broad_markers_present) > 0]

DotPlot(merged, features = broad_markers_present, cluster.idents = TRUE) +
  RotatedAxis() +
  ggtitle("Broad-type markers by cluster")
ggsave(file.path(output_dir, "03_broad_marker_dotplot.pdf"), width = 14, height = 7)

# Feature plots for a few key markers
key_genes <- c("KRT14", "COL1A1", "PTPRC", "PECAM1", "MLANA", "SOX10")
key_genes <- key_genes[key_genes %in% rownames(merged)]
if (length(key_genes) > 0) {
  FeaturePlot(merged, features = key_genes, ncol = 3, raster = FALSE)
  ggsave(file.path(output_dir, "04_key_feature_plots.pdf"), width = 14, height = 9)
}


################################################################################
# STEP 3 — Assign TIER 1 labels
#
#  >>> EDIT THE MAPPING BELOW <<<
#  Look at 02_all_cluster_markers.csv + 03_broad_marker_dotplot.pdf
#  and map each cluster number to a broad cell type.
#
################################################################################

tier1_map <- c(
  "0"  = "Fibroblast",
  "1"  = "Keratinocyte",
  "2"  = "Keratinocyte",
  "3"  = "Fibroblast",
  "4"  = "Immune",
  "5"  = "Endothelial",
  "6"  = "Fibroblast",
  "7"  = "Melanocyte",
  "8"  = "Smooth_Muscle",
  "9"  = "Neural_Crest"
  # ... add all your clusters
)

merged$broad_type <- tier1_map[as.character(merged$seurat_clusters)]

DimPlot(merged, group.by = "broad_type", label = TRUE, raster = FALSE) +
  ggtitle("Tier 1 — Broad cell types")
ggsave(file.path(output_dir, "05_tier1_broad_types.pdf"), width = 9, height = 7)

# Save checkpoint
saveRDS(merged, file.path(output_dir, "merged_tier1.rds"))
cat("Tier 1 annotated object saved.\n")


################################################################################
# STEP 4 — TIER 2: Sub-cluster a broad type
#
#  Two helper functions:
#    subcluster()    — subsets, re-processes, clusters at higher resolution
#    assign_subtypes() — writes the fine labels back into the main object
#
################################################################################

subcluster <- function(obj,
                       broad_label,
                       resolution  = 0.5,
                       dims        = 1:20,
                       output_dir  = "./manual_annotation") {
  
  cat("\n── Sub-clustering:", broad_label, "──\n")
  sub <- subset(obj, broad_type == broad_label)
  cat("   Cells:", ncol(sub), "\n")
  
  sub <- SCTransform(sub, assay = "RNA", verbose = FALSE)
  sub <- RunPCA(sub, npcs = 30, features = rownames(sub), verbose = FALSE)
  sub <- RunUMAP(sub, dims = dims, verbose = FALSE)
  sub <- FindNeighbors(sub, reduction = "pca", dims = dims, verbose = FALSE)
  sub <- FindClusters(sub, resolution = resolution, verbose = FALSE)
  
  # Markers
  sub_markers <- FindAllMarkers(
    sub, only.pos = TRUE, min.pct = 0.25,
    logfc.threshold = 0.3, verbose = FALSE
  )
  fname <- paste0("06_subcluster_markers_", broad_label, ".csv")
  write.csv(sub_markers, file.path(output_dir, fname), row.names = FALSE)
  
  top5 <- sub_markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = 5)
  cat("   Top markers per sub-cluster:\n")
  print(top5, n = 60)
  
  p <- DimPlot(sub, label = TRUE, raster = FALSE) +
    ggtitle(paste(broad_label, "sub-clusters"))
  print(p)
  ggsave(file.path(output_dir, paste0("07_subcluster_umap_", broad_label, ".pdf")),
         width = 8, height = 6)
  
  return(sub)
}


assign_subtypes <- function(main_obj, sub_obj, subtype_map) {
  # subtype_map: named character vector, cluster -> fine label
  # e.g. c("0" = "FIB_papillary", "1" = "FIB_reticular", ...)
  
  fine_labels <- subtype_map[as.character(sub_obj$seurat_clusters)]
  names(fine_labels) <- colnames(sub_obj)
  
  # Write into the main object
  idx <- match(names(fine_labels), colnames(main_obj))
  if (is.null(main_obj$cell_type)) {
    main_obj$cell_type <- main_obj$broad_type
  }
  main_obj$cell_type[idx] <- fine_labels
  
  cat("Assigned", length(fine_labels), "cells to subtypes:",
      paste(unique(fine_labels), collapse = ", "), "\n")
  return(main_obj)
}


# ──────────────────────────────────────────────────────────────────────────────
#  EXAMPLE: Sub-cluster fibroblasts
# ──────────────────────────────────────────────────────────────────────────────

fib_sub <- subcluster(merged, "Fibroblast", resolution = 0.5)

# Inspect known fibroblast sub-type markers from your panel:
fib_markers <- c("PDGFRA", "COL1A1", "COL3A1", "DCN", "LUM",
                 "POSTN", "COMP", "PRG4", "APCDD1", "WIF1",
                 "SFRP2", "SFRP4", "RSPO1", "DPP4", "PDPN")
fib_markers <- fib_markers[fib_markers %in% rownames(fib_sub)]
if (length(fib_markers) > 0) {
  DotPlot(fib_sub, features = fib_markers) + RotatedAxis()
  ggsave(file.path(output_dir, "08_fibroblast_subtype_dotplot.pdf"),
         width = 12, height = 5)
}

# >>> EDIT THIS MAPPING after inspecting the sub-cluster markers <<<
fib_subtype_map <- c(
  "0" = "FIB_upper",
  "1" = "FIB_lower",
  "2" = "FIB_origin",
  "3" = "FIB_deep"
  # ... add your sub-clusters
)
merged <- assign_subtypes(merged, fib_sub, fib_subtype_map)
rm(fib_sub); gc()


# ──────────────────────────────────────────────────────────────────────────────
#  EXAMPLE: Sub-cluster keratinocytes
# ──────────────────────────────────────────────────────────────────────────────

kc_sub <- subcluster(merged, "Keratinocyte", resolution = 0.4)

kc_markers <- c("KRT14", "KRT5", "KRT15", "KRT10", "KRT1",
                "LORICRIN", "FLG", "IVL", "LOR", "KRT6A",
                "KRT16", "KRT17", "MKI67", "TOP2A", "COL17A1")
kc_markers <- kc_markers[kc_markers %in% rownames(kc_sub)]
if (length(kc_markers) > 0) {
  DotPlot(kc_sub, features = kc_markers) + RotatedAxis()
  ggsave(file.path(output_dir, "09_keratinocyte_subtype_dotplot.pdf"),
         width = 12, height = 5)
}

# >>> EDIT THIS MAPPING <<<
kc_subtype_map <- c(
  "0" = "KC_basal",
  "1" = "KC_spinous",
  "2" = "KC_granular",
  "3" = "KC_proliferating"
  # ... add your sub-clusters
)
merged <- assign_subtypes(merged, kc_sub, kc_subtype_map)
rm(kc_sub); gc()


# ──────────────────────────────────────────────────────────────────────────────
#  Repeat for any other broad types you want to subdivide:
#
#  immune_sub <- subcluster(merged, "Immune", resolution = 0.6)
#  immune_map <- c("0" = "T_cell", "1" = "Macrophage", ...)
#  merged <- assign_subtypes(merged, immune_sub, immune_map)
# ──────────────────────────────────────────────────────────────────────────────


################################################################################
# STEP 5 — Final overview + save
################################################################################

# Any cells not sub-clustered retain their broad_type label
if (is.null(merged$cell_type)) {
  merged$cell_type <- merged$broad_type
} else {
  still_broad <- is.na(merged$cell_type) | merged$cell_type == merged$broad_type
  merged$cell_type[still_broad] <- merged$broad_type[still_broad]
}

DimPlot(merged, group.by = "cell_type", label = TRUE,
        repel = TRUE, raster = FALSE) +
  ggtitle("Final annotation") +
  NoLegend()
ggsave(file.path(output_dir, "10_final_annotation.pdf"), width = 10, height = 8)

# Composition table
comp <- table(merged$cell_type, merged$sample)
write.csv(as.data.frame.matrix(comp),
          file.path(output_dir, "11_cell_type_composition.csv"))
cat("\nCell-type composition saved.\n")
print(comp)

# Save final object
saveRDS(merged, file.path(output_dir, "merged_annotated_final.rds"))
cat("\nDone. Final object saved to", file.path(output_dir, "merged_annotated_final.rds"), "\n")