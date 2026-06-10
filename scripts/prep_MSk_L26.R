#!/usr/bin/env Rscript
suppressMessages({library(Seurat); library(ggplot2); library(schard)})

H5AD    <- "/mnt/d/scRNA_datasets/HYMS_metal_diabetes/Reference/reconstructed_skin_reference.h5ad"
RDS_OUT <- sub("\\.h5ad$", ".rds", H5AD)
CLUSTER_COL <- "leiden_scVI_0.9"

# full-gene log-norm matrix (adata.raw); save Seurat object
seu  <- schard::h5ad2seurat(H5AD, use.raw = TRUE)
expr <- LayerData(seu, layer = "counts")
saveRDS(seu, RDS_OUT)

markers <- c("Crabp1","Twist2","Pdgfra","Dpp4","Wif1","Apcdd1","Dlk1","Agtr2",
             "Fabp4","Mfap5","Ebf2","Meox2","Igfbp7","Col11a1","Cd200","Acan",
             "Sox18","Cxcr4","Sox2","Alpl","Alx4","Col23a1","Actg2","Itga8",
             "Coch","Matn4","Dlx5","Trp63","Lhx2","Edar","Barx2","Sox9",
             "Il11ra1","Krt79","Sostdc1","Apoe","Msx2","Krt71","Krt5","Krt14",
             "Krt1","Krt10","Cpa3","Mcpt4","Cxcr2","Itgam","Ptprc","Cd3g",
             "Cd68","Cd86","Mrc1","Cd163","Krt8","Krt18","Msc","Ttn","Pax7",
             "Rgs5","Sox10","Dct","Tyr","Pecam1","Cdh5","Lyve1")
markers <- markers[markers %in% rownames(expr)]

grp <- as.character(seu@meta.data[[CLUSTER_COL]])
g0  <- factor(grp, levels = sort(unique(grp)))

# per-gene z-score across all cells (matches sc.pp.scale), then per-cluster mean
m <- as.matrix(expr[markers, , drop = FALSE])
z <- t(scale(t(m))); z[is.na(z)] <- 0
color_mat <- sapply(levels(g0), function(g) rowMeans(z[, g0 == g, drop = FALSE]))
pct_mat   <- sapply(levels(g0), function(g) rowMeans(m[, g0 == g, drop = FALSE] > 0)) * 100

# order clusters the way the paper does: by the marker-axis center-of-mass of their
# positive enrichment (left-peaking clusters on top -> diagonal), no labels needed
w     <- pmax(color_mat, 0)
denom <- colSums(w)
pos   <- ifelse(denom > 0, colSums(w * seq_along(markers)) / denom,
                apply(color_mat, 2, which.max))
lev   <- names(sort(pos))
color_mat <- color_mat[, lev]; pct_mat <- pct_mat[, lev]

df <- expand.grid(gene = factor(markers, levels = markers),
                  group = factor(lev, levels = lev), stringsAsFactors = FALSE)
df$meanz <- as.vector(color_mat)
df$pct   <- as.vector(pct_mat)

write.csv(df, file = "ref_markers.csv")

rdbu_r <- c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7",
            "#fddbc7","#f4a582","#d6604d","#b2182b","#67001f")
print(
  ggplot(df, aes(gene, group)) +
    geom_point(aes(color = meanz, size = pct)) +
    scale_color_gradientn(colours = rdbu_r, limits = c(-2.5, 2.5),
                          oob = scales::squish, name = "mean z-score") +
    scale_size(range = c(0, 6), name = "% expressing") +
    scale_y_discrete(limits = rev(lev)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
          axis.title = element_blank())
)

cluster_names <- c(
  "0" = "BK-2", "1" = "FIB-1", "2" = "FIB-2", "3" = "emFIB-2",
  "4" = "BK-3", "5" = "HP/HG", "6" = "BK-1", "7" = "DC/DP",
  "8" = "Spinous", "9" = "FIB-3", "10" = "FIB-4", "11" = "emFIB-1",
  "12" = "Pre-adipo", "13" = "DP", "14" = "FIB-2", "15" = "Muscle",
  "16" = "emK", "17" = "emK", "18" = "MAC-2", "19" = "Mast-cell",
  "20" = "APM", "21" = "Chond-Fib", "22" = "Lymphocyte", "23" = "Pericyte",
  "24" = "Schwann", "25" = "MAC-1", "26" = "Melanocyte", "27" = "Coch-Fib",
  "28" = "VEC", "29" = "LEC", "30" = "Neutrophil", "31" = "Merkel"
)

seu$cell_type <- unname(cluster_names[as.character(seu@meta.data[[CLUSTER_COL]])])
DimPlot(seu, group.by = "cell_type", label = TRUE)
FeaturePlot(seu, features = "Tagln")

saveRDS(seu, file = "./reference/Lee26_reference.rds")

