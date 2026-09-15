################################################################################
# panel_design_functions.R
#
# helpers for design_panel.R
#
# membership -> new genes -> population detection -> anchor genes -> export
################################################################################

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})


#######################################################################
# candidate list -> long table                                        #
#######################################################################

#one row per gene, first block kept when a gene is listed twice
candidates_to_table <- function(candidates) {
  tab <- data.frame(
    gene  = unlist(candidates, use.names = FALSE),
    block = rep(names(candidates), lengths(candidates)),
    stringsAsFactors = FALSE
  )
  tab[!duplicated(tab$gene), ]
}


#######################################################################
# 1. on / off the 5K panel                                            #
#######################################################################

check_panel_membership <- function(obj, candidates, assay = "RNA") {
  panel_genes <- rownames(obj[[assay]])
  tab <- candidates_to_table(candidates)
  tab$on_5k <- tab$gene %in% panel_genes
  
  cat("on 5K     (", sum(tab$on_5k),  "): ",
      paste(tab$gene[tab$on_5k],  collapse = ", "), "\n\n", sep = "")
  cat("not on 5K (", sum(!tab$on_5k), "): ",
      paste(tab$gene[!tab$on_5k], collapse = ", "), "\n\n", sep = "")
  
  tab
}


#######################################################################
# 2. new genes                                                        #
#######################################################################

#off-5K candidates only. slots per block are proportional to how many of
#that block's genes are off the panel; within a block, list order wins
select_new_genes <- function(membership, block_reason, n_new = 100) {
  off <- membership[!membership$on_5k, ]
  n_new <- min(n_new, nrow(off))
  
  blocks <- unique(off$block)
  share  <- n_new * as.numeric(table(factor(off$block, levels = blocks))) / nrow(off)
  quota  <- floor(share)
  extra  <- n_new - sum(quota)
  if (extra > 0) {
    top <- order(share - quota, decreasing = TRUE)[seq_len(extra)]
    quota[top] <- quota[top] + 1
  }
  names(quota) <- blocks
  
  picked <- do.call(rbind, lapply(blocks, function(b) {
    head(off[off$block == b, ], quota[[b]])
  }))
  
  data.frame(
    gene       = picked$gene,
    set        = "New",
    tag        = "Required",
    population = picked$block,
    pct_in     = NA_real_,
    pct_other  = NA_real_,
    reason     = unname(block_reason[picked$block]),
    stringsAsFactors = FALSE
  )
}


#######################################################################
# 3. per-population detection from the split counts layers            #
#######################################################################

#% of cells in each population with >= 1 count, summed over every
#counts.<sample> layer. populations missing from a sample are kept as
#zero columns instead of being dropped
population_detection <- function(obj,
                                 label_col = "cell_type",
                                 assay     = "RNA",
                                 exclude   = "Rem\\.|LowQual|Unk|Undet|Unassigned|Rescue") {
  
  labels <- setNames(as.character(obj[[]][[label_col]]), colnames(obj))
  keep   <- names(labels)[!is.na(labels) & !grepl(exclude, labels)]
  types  <- sort(unique(labels[keep]))
  genes  <- rownames(obj[[assay]])
  
  hits  <- matrix(0, nrow = length(genes), ncol = length(types),
                  dimnames = list(genes, types))
  ncell <- setNames(numeric(length(types)), types)
  
  layers <- grep("^counts", Layers(obj[[assay]]), value = TRUE)
  for (layer in layers) {
    counts <- LayerData(obj, assay = assay, layer = layer)
    cells  <- intersect(colnames(counts), keep)
    if (!length(cells)) next
    
    group <- fac2sparse(factor(labels[cells], levels = types),
                        drop.unused.levels = FALSE)          # types x cells
    detected <- counts[, cells, drop = FALSE] > 0            # explicit zeros ignored
    layer_hits <- as.matrix(detected %*% t(group))           # genes x types
    
    hits[rownames(layer_hits), colnames(layer_hits)] <-
      hits[rownames(layer_hits), colnames(layer_hits)] + layer_hits
    ncell[rownames(group)] <- ncell[rownames(group)] + Matrix::rowSums(group)
    
    cat(layer, ":", length(cells), "cells\n")
  }
  
  pct <- 100 * sweep(hits, 2, ncell, "/")
  
  #highest % among all OTHER populations, per gene and population
  first  <- apply(pct, 1, max)
  second <- apply(pct, 1, function(x) sort(x, decreasing = TRUE)[2])
  other  <- matrix(first, nrow(pct), ncol(pct), dimnames = dimnames(pct))
  is_top <- pct == first
  other[is_top] <- matrix(second, nrow(pct), ncol(pct))[is_top]
  
  list(pct = pct, other = other, specificity = pct - other, ncell = ncell)
}


#######################################################################
# 4. anchor genes                                                     #
#######################################################################

#on-5K genes, taken round-robin across populations by specificity
#(% in population minus % in next-highest population). on-5K candidate
#genes are always kept. the n_optional least specific are Optional
select_anchor_genes <- function(detection,
                                membership,
                                n_anchor   = 150,
                                n_optional = 50,
                                depth      = 30,
                                exclude    = character(0)) {
  
  pct   <- detection$pct
  other <- detection$other
  spec  <- detection$specificity
  types <- colnames(pct)
  
  anchor_row <- function(gene, type) {
    data.frame(gene = gene, population = type,
               pct_in = pct[gene, type], pct_other = other[gene, type],
               specificity = spec[gene, type], stringsAsFactors = FALSE)
  }
  
  forced_genes <- membership$gene[membership$on_5k]
  forced <- do.call(rbind, lapply(forced_genes, function(g) {
    anchor_row(g, types[which.max(spec[g, ])])
  }))
  
  ranked <- do.call(rbind, lapply(types, function(type) {
    top <- head(order(spec[, type], decreasing = TRUE), depth)
    cbind(anchor_row(rownames(spec)[top], type), rank = seq_along(top))
  }))
  ranked <- ranked[ranked$specificity > 0, ]
  ranked <- ranked[order(ranked$rank, -ranked$specificity), names(forced)]
  
  anchors <- rbind(forced, ranked)
  anchors <- anchors[!duplicated(anchors$gene) & !anchors$gene %in% exclude, ]
  anchors <- head(anchors, n_anchor)
  
  optional <- rank(anchors$specificity, ties.method = "first") <= n_optional
  
  data.frame(
    gene       = anchors$gene,
    set        = "Anchor",
    tag        = ifelse(optional, "Optional", "Required"),
    population = anchors$population,
    pct_in     = round(anchors$pct_in, 2),
    pct_other  = round(anchors$pct_other, 2),
    reason     = sprintf("%s%s marker: %.1f%% of cells vs %.1f%% in next-highest population",
                         ifelse(anchors$gene %in% forced_genes, "candidate on 5K; ", ""),
                         anchors$population, anchors$pct_in, anchors$pct_other),
    stringsAsFactors = FALSE
  )
}


#######################################################################
# 5. export                                                           #
#######################################################################

write_panel_design <- function(panel,
                               membership,
                               detection,
                               dir,
                               format = c("csv", "xlsx")) {
  
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    proposed_panel    = panel,
    panel_membership  = membership,
    population_pct    = data.frame(gene = rownames(detection$pct),
                                   round(detection$pct, 3), check.names = FALSE)
  )
  
  if ("csv" %in% format) {
    for (nm in names(tables))
      write.csv(tables[[nm]], file.path(dir, paste0(nm, ".csv")), row.names = FALSE)
  }
  if ("xlsx" %in% format) {
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(tables, file.path(dir, "panel_design.xlsx"))
    } else {
      warning("openxlsx not installed; xlsx skipped")
    }
  }
  
  print(table(panel$set, panel$tag))
  invisible(tables)
}

################################################################################
# design_panel.R
#
# candidate genes -> on/off 5K -> new genes -> population detection
#                 -> anchor genes -> export
#
# 250-gene panel = ~100 new genes (not on the 5K) + up to 150 anchors (on the
# 5K). the 50 least specific anchors are Optional, i.e. the slots that can go
# to metal-interaction genes while still keeping 100 anchors.
################################################################################

output_dir <- "./annotated_data"
sketch_dir <- file.path(output_dir, "_sketch")
design_dir <- file.path(output_dir, "_panel_design")


#######################################################################
#######################################################################
#                                                                     #
#                          1. LOAD OBJECT                             #
#                                                                     #
#######################################################################
#######################################################################

obj <- readRDS(file.path(sketch_dir, "merged_annotated.rds"))
DefaultAssay(obj) <- "RNA"


#######################################################################
#######################################################################
#                                                                     #
#                        2. CANDIDATE GENES                           #
#                                                                     #
#######################################################################
#######################################################################

#block order and within-block order = priority for the new-gene slots
candidates <- list(
  A_Keratins     = c("Krt5","Krt14","Krt15","Krt1","Krt10","Krtdap","Krt77",
                     "Krt6a","Krt6b","Krt16","Krt17","Krt75",
                     "Krt71","Krt25","Krt26","Krt27","Krt28","Krt73",
                     "Krt32","Krt35","Krt82","Krt84",
                     "Krt31","Krt33a","Krt33b","Krt34","Krt81","Krt83","Krt85","Krt86",
                     "Krt80","Krt7","Krt8","Krt18","Krt19","Krt20",
                     "Krtap3-1","Krtap15","Krtap16-3"),
  B_Barrier      = c("Lor","Flg","Flg2","Ivl","Sbsn","Cnfn","Tgm1","Tgm3","Cdsn",
                     "Asprv1","Rptn","Hrnr","Lce1a1","Klk7","Klk8","Spink5",
                     "Casp14","Abca12","Nipal4","Cers3"),
  C_HairFollicle = c("Shh","Lef1","Msx2","Dlx3","Hoxc13","Foxn1","Tchh","Sostdc1",
                     "Bmp2","Bmp4","Dkk1","Rspo3","Wnt5a","Corin","Sox2","Wif1",
                     "Alx4","Enpp2","Nfatc1","Gli1","Edar","Foxi3","Shisa2",
                     "Prdm1","Plet1"),
  D_Melanocyte   = c("Mitf","Tyr","Pmel","Mlana","Slc45a2","Slc24a5","Oca2",
                     "Trpm1","Mc1r","Pax3","Kitl"),
  E_NerveMerkel  = c("Plp1","Mpz","Mbp","Pmp22","Prx","Cdh19","Scn7a","Foxd3",
                     "Egr2","Ncmap","Cldn19","Pllp","Apod","Ptn","Gfra3","Nrxn1",
                     "Prph","Nefl","Nefm","Calca","Tac1","Ntrk1","Scn9a","Trpv1",
                     "Mrgprd","Snap25","Atoh1","Isl1","Chga","Syp"),
  F_Lymphocyte   = c("Cd3d","Cd3g","Cd4","Cd8a","Cd8b1","Foxp3","Ikzf2","Ctla4",
                     "Il2ra","Tnfrsf4","Trdc","Tcrg-C1","Trgv5","Rorc","Il17a",
                     "Ifng","Nkg7","Prf1","Gzmb","Klrb1c","Eomes","Cd79a","Cd79b",
                     "Ms4a1","Jchain"),
  G_Granulocyte  = c("S100a8","S100a9","Ly6g","Retnlg","Mpo","Elane","Ltf","Camp",
                     "Ngp","Padi4"),
  H_DC           = c("Cd207","Xcr1","Clec9a","Zbtb46","Mgl2","H2-Ab1","H2-Eb1","Cd74"),
  I_Mural        = c("Rgs5","Higd1b","Anpep","Acta2","Cnn1","Lmod1","Mylk","Pln",
                     "Steap4","Art3","Foxs1","Ndufa4l2"),
  J_Fibroblast   = c("Pi16","Ly6a","En1","Dlk1","Ptgds","Fmo1","Mgp","Gpx3",
                     "Col6a5","Coch","Ndnf","Wnt2","Comp","Tnmd"),
  K_Adipose      = c("Cfd","Fabp4","Cidec","Car3","Ebf2","Wt1","Zfp423","Ucp1"),
  L_Open         = c("Cthrc1","Wisp1","Nrg1","Areg","Tgfb2","Fgf2","Pax7","Myf5",
                     "Calcr","Foxa1","Best2","Ano1","Dcpp1","Oxtr")
)

block_reason <- c(
  A_Keratins     = "Keratin: basal/suprabasal/wound/follicle/Merkel identity (panel is keratin-excluded)",
  B_Barrier      = "Granular/cornified strata, currently inferred from lipid enzymes",
  C_HairFollicle = "Follicle compartments; dermal papilla vs sheath (Fib.HFDermal)",
  D_Melanocyte   = "Separates melanocytes from matrix progenitors in HF.Bulb",
  E_NerveMerkel  = "Schwann myelin/Remak, sensory axons, Merkel; reinnervation readout",
  F_Lymphocyte   = "Splits Lymph.TNK into T subsets, gdT/DETC, NK, B/plasma",
  G_Granulocyte  = "Confirms Neut.Sell / Neut.Mmp9",
  H_DC           = "cDC1 vs cDC2; confirms Langerhans cells",
  I_Mural        = "Pericyte vs vSMC vs arrector pili; Mural.Myofibro vs Fib.Activated",
  J_Fibroblast   = "Papillary vs reticular vs fascia lineage",
  K_Adipose      = "Preadipocyte / adipocyte states",
  L_Open         = "Wound / satellite cell / gland markers"
)


#######################################################################
#######################################################################
#                                                                     #
#                     3. ON / OFF THE 5K PANEL                        #
#                                                                     #
#######################################################################
#######################################################################

membership <- check_panel_membership(obj, candidates, assay = "RNA")


#######################################################################
#######################################################################
#                                                                     #
#                           4. NEW GENES                              #
#                                                                     #
#######################################################################
#######################################################################

new_genes <- select_new_genes(membership,
                              block_reason = block_reason,
                              n_new        = 100)


#######################################################################
#######################################################################
#                                                                     #
#                          5. ANCHOR GENES                            #
#                                                                     #
#######################################################################
#######################################################################

#% detecting cells per cell_type, summed over the counts.<sample> layers.
#Rem./LowQual/Unk/Undet/Unassigned/Rescue labels are left out
detection <- population_detection(obj,
                                  label_col = "cell_type",
                                  assay     = "RNA",
                                  exclude   = "Rem\\.|LowQual|Unk|Undet|Unassigned|Rescue")

anchor_genes <- select_anchor_genes(detection,
                                    membership = membership,
                                    n_anchor   = 150,
                                    n_optional = 50,
                                    depth      = 30,
                                    exclude    = new_genes$gene)


#######################################################################
#######################################################################
#                                                                     #
#                            6. EXPORT                                #
#                                                                     #
#######################################################################
#######################################################################

panel <- rbind(new_genes, anchor_genes)

write_panel_design(panel,
                   membership = membership,
                   detection  = detection,
                   dir        = design_dir,
                   format     = c("csv", "xlsx"))