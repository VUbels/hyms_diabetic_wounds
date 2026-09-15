library(SeuratObject); library(Matrix)
# obj: Seurat v5, assay "RNA" (split counts.<sample> layers), meta column "cell_type"
# supp_panel: named list of §11 blocks (block order and within-block order = priority)

n_new <- 100; n_anchor <- 150; n_optional <- 50
panel_genes <- rownames(obj[["RNA"]])

supp_panel <- list(
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
                     "Trpm1","Mc1r","Pax3","Kitl","Pmel"),
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
  A_Keratins = "Keratin: basal/suprabasal/wound/follicle/Merkel identity (panel is keratin-excluded)",
  B_Barrier = "Granular/cornified strata, currently inferred from lipid enzymes",
  C_HairFollicle = "Follicle compartments; dermal papilla vs sheath (Fib.HFDermal)",
  D_Melanocyte = "Separates melanocytes from matrix progenitors in HF.Bulb",
  E_NerveMerkel = "Schwann myelin/Remak, sensory axons, Merkel; reinnervation readout",
  F_Lymphocyte = "Splits Lymph.TNK into T subsets, gdT/DETC, NK, B/plasma",
  G_Granulocyte = "Confirms Neut.Sell / Neut.Mmp9",
  H_DC = "cDC1 vs cDC2; confirms Langerhans cells",
  I_Mural = "Pericyte vs vSMC vs arrector pili; Mural.Myofibro vs Fib.Activated",
  J_Fibroblast = "Papillary vs reticular vs fascia lineage",
  K_Adipose = "Preadipocyte / adipocyte states",
  L_Open = "Wound / satellite cell / gland markers")

# ---- 1. Sub-panel genes on vs not on the 5K --------------------------------
cand <- unique(unlist(supp_panel))
on_5k  <- cand[cand %in% panel_genes]
off_5k <- cand[!cand %in% panel_genes]

# ---- 2. New genes: off-5K candidates, quota per block proportional to size --
nb <- data.frame(gene = unlist(supp_panel, use.names = FALSE),
                 block = rep(names(supp_panel), lengths(supp_panel)))
nb <- nb[!duplicated(nb$gene) & !nb$gene %in% panel_genes, ]
n_new <- min(n_new, nrow(nb))
q <- n_new * table(factor(nb$block, unique(nb$block))) / nrow(nb)
quota <- floor(q)
extra <- n_new - sum(quota)
if (extra > 0) { i <- order(q - quota, decreasing = TRUE)[seq_len(extra)]; quota[i] <- quota[i] + 1 }
new_df <- do.call(rbind, lapply(names(quota), function(b) head(nb[nb$block == b, ], quota[[b]])))
new_df <- data.frame(gene = new_df$gene, set = "New", tag = "Required",
                     population = new_df$block, reason = block_reason[new_df$block])

# ---- 3. Anchors: on-5K genes most specific to each annotated population -----
ct   <- setNames(obj$cell_type, colnames(obj))
keep <- names(ct)[!is.na(ct) & !grepl("Rem\\.|LowQual|Unk|Rescue", ct)]
types <- sort(unique(ct[keep]))
hits  <- matrix(0, length(panel_genes), length(types), dimnames = list(panel_genes, types))
ncell <- setNames(numeric(length(types)), types)
for (l in grep("^counts", Layers(obj[["RNA"]]), value = TRUE)) {
  m  <- LayerData(obj, assay = "RNA", layer = l)
  cl <- intersect(colnames(m), keep)
  f  <- fac2sparse(factor(ct[cl], levels = types))          # types x cells
  b  <- m[, cl, drop = FALSE]; b@x[] <- 1
  h  <- as.matrix(b %*% t(f))
  hits[rownames(h), ] <- hits[rownames(h), ] + h            # match by gene name
  ncell <- ncell + Matrix::rowSums(f)
}
pct   <- 100 * sweep(hits, 2, ncell, "/")
top1  <- apply(pct, 1, max)
top2  <- apply(pct, 1, function(x) sort(x, decreasing = TRUE)[2])
other <- matrix(top1, nrow(pct), ncol(pct), dimnames = dimnames(pct))
is1   <- pct == top1
other[is1] <- matrix(top2, nrow(pct), ncol(pct))[is1]
spec  <- pct - other                                        # % in type minus best other type

row_for <- function(g, t) data.frame(gene = g, population = t, pct_in = pct[g, t],
                                     pct_other = other[g, t], spec = spec[g, t])
forced <- do.call(rbind, lapply(on_5k, function(g) row_for(g, types[which.max(spec[g, ])])))
ranked <- do.call(rbind, lapply(types, function(t) {
  g <- names(sort(spec[, t], decreasing = TRUE))[1:30]
  cbind(row_for(g, t), rank = 1:30)
}))
ranked <- ranked[ranked$spec > 0, ]
ranked <- ranked[order(ranked$rank, -ranked$spec), names(forced)]
anc <- rbind(forced, ranked)
anc <- head(anc[!duplicated(anc$gene) & !anc$gene %in% new_df$gene, ], n_anchor)

anchor_df <- data.frame(
  gene = anc$gene, set = "Anchor (on 5K)",
  tag = ifelse(rank(anc$spec, ties.method = "first") <= n_optional, "Optional", "Required"),
  population = anc$population,
  reason = sprintf("%s%s marker: %.1f%% of cells vs %.1f%% in next-highest population",
                   ifelse(anc$gene %in% on_5k, "Sub-panel gene; ", ""),
                   anc$population, anc$pct_in, anc$pct_other))

# ---- Output ------------------------------------------------------------------
final_panel <- rbind(new_df, anchor_df); rownames(final_panel) <- NULL
cat("On 5K (", length(on_5k), "):", paste(on_5k, collapse = ", "), "\n\n")
cat("Not on 5K (", length(off_5k), "):", paste(off_5k, collapse = ", "), "\n\n")
print(table(final_panel$set, final_panel$tag))
write.csv(final_panel, "proposed_250_panel.csv", row.names = FALSE)