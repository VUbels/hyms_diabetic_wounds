################################################################################
# prepare_objects.R
#
# proseg -> merge -> sketch -> broad labels -> subcluster -> fine labels
#        -> rescue Rem. clusters -> projection -> export
#
# Sections 3 and 5 each stop for manual labelling. Both offer the same two
# routes: option A labels in this script, option B labels in a csv.
################################################################################

source("./scripts/setup_py_environment.R")
setup_py_env("hyms_metal_diabetes", "/home/uvictor/miniconda3/condabin/conda")

source("./scripts/helper_functions.R")
source("./scripts/convert_to_h5ad.R")

reticulate::py_run_string("
import sys
sys.path.append('./scripts')
")

#zarr converter used by build_proseg_seurat
reticulate::py_run_string("import ssl")
zconv <- reticulate::import_from_path("zarr_h5ad_conversion", path = "./scripts")

output_dir <- "./annotated_data"
sketch_dir <- file.path(output_dir, "_sketch")
dir.create(sketch_dir, recursive = TRUE, showWarnings = FALSE)


#######################################################################
#######################################################################
#                                                                     #
#                       1. BUILD REGION OBJECTS                       #
#                                                                     #
#######################################################################
#######################################################################

# build_proseg_seurat(
#   proseg_dir = "./proseg_results_mask",
#   xenium_dir = "/mnt/d/HYMS/metal_diabetes",
#   overwrite  = FALSE
# )

sample_sheet <- read_sample_sheet("./samples.csv")
regions      <- regions_from_sheet(sample_sheet)

check_proseg_regions(regions)


#######################################################################
#######################################################################
#                                                                     #
#                        2. MERGE AND SKETCH                          #
#                                                                     #
#######################################################################
#######################################################################

#on_disk = TRUE keeps the full counts matrix in BPCells instead of memory
obj <- merge_proseg_regions(
  regions      = regions,
  sample_sheet = sample_sheet,
  min_counts   = 10,
  on_disk      = FALSE
)

obj <- sketch_proseg_object(
  obj,
  sketch_cells = 20000,
  method       = "leverage",
  n_features   = 2500
)

obj <- cluster_sketch(
  obj,
  dims        = 1:30,
  resolutions = c(0.3, 0.5, 0.8, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3),
  resolution  = 0.5,
  integrate   = FALSE
)

saveRDS(obj, file.path(sketch_dir, "merged_sketched.rds"))


#######################################################################
#######################################################################
#                                                                     #
#                    3. BROAD COMPARTMENT LABELS                      #
#                                                                     #
#######################################################################
#######################################################################

markers_broad <- sketch_cluster_markers(
  obj,
  cluster_col  = "sketch_cluster",
  sketch_assay = "sketch",
  output_dir   = output_dir,
  format       = c("csv", "xlsx")
)


top_cluster_markers(markers_broad, 15)

#######################################################################
# OPTION A: MANUAL BROAD CLUSTER ANNOTATION                           #
#######################################################################

broad_labels <- c(
  "0"  = "Stromal",        # Col1a1, Postn, Pdgfra        (activated fibroblast)
  "1"  = "Stromal",        # Col1a1, Pdgfra, Col14a1      (reticular/ECM fibroblast)
  "2"  = "Epithelial",     # Dsg1a, Dsc1, Gjb2            (suprabasal keratinocyte)
  "3"  = "Immune",         # Cd68, Adgre1, Itgam          (macrophage)
  "4"  = "Endothelial",    # Pecam1, Cdh5, Egfl7          (blood EC)
  "5"  = "Epithelial",     # Sox9, Epcam, Trp63           (hair follicle / JZ)
  "6"  = "Immune",         # Cd14, Itgb2, Csf3r           (inflammatory myeloid)
  "7"  = "Immune",         # Ptprc, Itgb2, Cd53           (DC / lymphoid)
  "8"  = "Immune",         # Csf1r, Adgre1, Mrc1          (resident macrophage)
  "9"  = "Mural_Muscle",   # Pdgfrb, Notch3, Tagln        (pericyte/SMC)
  "10" = "Epithelial",     # Itgb4, Itga6, Trp63          (basal keratinocyte)
  "11" = "Stromal",        # Col1a1, Col1a2, Mmp2         (fibro/macrophage doublets)
  "12" = "Stromal",        # Adipoq, Plin1, Lpl           (dermal adipocyte)
  "13" = "Mural_Muscle",   # Des, Ttn, Actn2              (skeletal myofibre)
  "14" = "Epithelial",     # Cidea, Far2, Elovl4          (sebocyte)
  "15" = "NeuralCrest",    # Sox10, Ngfr, Erbb3           (Schwann/glia)
  "16" = "Endothelial",    # Pecam1, Egfl7, Prox1         (lymphatic EC)
  "17" = "Immune",         # Mki67, Top2a, Cd84           (cycling, provisional)
  "18" = "Mural_Muscle"    # Myh3, Mymk, Actc1            (myoblast/fusing fibre)
)

obj <- label_sketch_clusters(obj, broad_labels,
                             cluster_col = "sketch_cluster",
                             label_col   = "broad_type")

#######################################################################
# OPTION B: LABEL IN ./annotated_data/_sketch/cluster_annotation.csv  #
#######################################################################
# option B overwrites broad_type with whatever is in the csv, so leave it
# commented out while option A above is in use

# write_cluster_template(obj, markers_broad,
#                        cluster_col = "sketch_cluster",
#                        output_dir  = output_dir)
#
# obj <- label_sketch_clusters(obj,
#                              cluster_col = "sketch_cluster",
#                              label_col   = "broad_type",
#                              output_dir  = output_dir)


#builds pca.full and full.umap, then spreads broad_type to every cell
obj <- project_sketch_labels(obj, label_col = "broad_type", dims = 1:30)

saveRDS(obj, file.path(sketch_dir, "merged_broad.rds"))


#######################################################################
#######################################################################
#                                                                     #
#                4. SUBCLUSTER EACH COMPARTMENT                       #
#                                                                     #
#######################################################################
#######################################################################

subclusters <- run_subclustering(
  obj,
  broad_col    = "broad_type",
  output_dir   = output_dir,
  skip_types   = NULL,
  min_cells    = 500,
  dims         = 1:20,
  resolutions  = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
  resolution   = 0.6,
  n_features   = 2000,
  sketch_cells = 20000,
  direct_max   = 75000,
  markers      = FALSE,   #deferred until the resolution is chosen
  overwrite    = FALSE
)


#######################################################################
#######################################################################
#                                                                     #
#                        5. FINE CELL TYPES                           #
#                                                                     #
#######################################################################
#######################################################################

###----------------------------------------------------------------------------
### 5a. one compartment at a time: choose the resolution, then get markers
###----------------------------------------------------------------------------

#--- Epithelial ---------------------------------------------------------------

#prints the resolutions available and their cluster counts
epi <- load_subcluster("Epithelial", output_dir)

#compare shapes before committing
DimPlot(epi$obj, group.by = epi$resolutions$column, label = TRUE)
DimPlot(epi$obj, group.by = "RNA_snn_res.1", label = TRUE)

#try a candidate resolution. dir = NULL means nothing is written to disk,
#so this is throwaway. give it a dir plus format to get a file out
trial <- subcluster_markers(epi$obj,
                            cluster_col = "RNA_snn_res.1",
                            dir         = epi$dir,
                            format      = c("csv", "xlsx"),
                            only_pos    = TRUE,
                            logfc       = 0.2,
                            min_pct     = 0.2)

table(trial$cluster)              #clusters at 0 have nothing separating them
top_cluster_markers(trial, 15)

#trial is one vs rest, so two similar clusters share genes against the
#pooled remainder. test suspicious pairs against each other directly
compare_clusters(epi$obj, 11, 17, cluster_col = "RNA_snn_res.1")

#or do it systematically: every cluster against its closest neighbour
local_markers <- neighbour_markers(epi$obj,
                                   cluster_col = "RNA_snn_res.1",
                                   reduction   = "pca",
                                   dims        = 1:20,
                                   k           = 1,
                                   logfc       = 0.2,
                                   min_pct     = 0.1,
                                   dir         = epi$dir,
                                   format      = c("csv", "xlsx"))

attr(local_markers, "pairs")              #who was compared with whom
subset(local_markers, cluster == "11")

#lock it in: sets subcluster, recomputes markers at that resolution,
#rewrites annotation.csv and saves the object back to disk
epi <- finalise_subcluster(epi, resolution = 1, format = c("csv", "xlsx"))
epi$top


#--- Stromal ------------------------------------------------------------------

#res 0.4 was under-resolved, 0.8 splits the inflammatory/proliferative
#composite and separates adipocytes from the adipocyte-macrophage mixture
str <- load_subcluster("Stromal", output_dir)
DimPlot(str$obj, group.by = "RNA_snn_res.0.8", label = TRUE)

str <- finalise_subcluster(str, resolution = 0.8, format = c("csv", "xlsx"))
str_local <- neighbour_markers(str$obj,
                               cluster_col = "subcluster",
                               reduction   = "pca",
                               dims        = 1:20,
                               k           = 1,
                               logfc       = 0.2,
                               min_pct     = 0.1,
                               dir         = str$dir,
                               format      = c("csv", "xlsx"))
str$top


#--- Immune -------------------------------------------------------------------

#22 clusters (0-21)
imm <- load_subcluster("Immune", output_dir)
DimPlot(imm$obj, group.by = imm$resolutions$column, label = TRUE)

#resolution not recorded in the chats: pick the one giving 22 clusters
imm <- finalise_subcluster(imm, resolution = 1.2, format = c("csv", "xlsx"))
imm_local <- neighbour_markers(imm$obj,
                               cluster_col = "subcluster",
                               reduction   = "pca",
                               dims        = 1:20,
                               k           = 1,
                               logfc       = 0.2,
                               min_pct     = 0.1,
                               dir         = imm$dir,
                               format      = c("csv", "xlsx"))
imm$top


#--- Endothelial --------------------------------------------------------------

#8 clusters (0-7)
endo <- load_subcluster("Endothelial", output_dir)
DimPlot(endo$obj, group.by = endo$resolutions$column, label = TRUE)

#resolution not recorded in the chats: pick the one giving 8 clusters
endo <- finalise_subcluster(endo, resolution = 0.6, format = c("csv", "xlsx"))
endo_local <- neighbour_markers(endo$obj,
                                cluster_col = "subcluster",
                                reduction   = "pca",
                                dims        = 1:20,
                                k           = 1,
                                logfc       = 0.2,
                                min_pct     = 0.1,
                                dir         = endo$dir,
                                format      = c("csv", "xlsx"))
endo$top


#--- Mural_Muscle -------------------------------------------------------------

#10 clusters (0-9), rebuilt after the Mural-Muscle projection split
mus <- load_subcluster("Mural_Muscle", output_dir)
DimPlot(mus$obj, group.by = mus$resolutions$column, label = TRUE)

#resolution not recorded in the chats: pick the one giving 10 clusters
mus <- finalise_subcluster(mus, resolution = 0.6, format = c("csv", "xlsx"))
mus_local <- neighbour_markers(mus$obj,
                               cluster_col = "subcluster",
                               reduction   = "pca",
                               dims        = 1:20,
                               k           = 1,
                               logfc       = 0.2,
                               min_pct     = 0.1,
                               dir         = mus$dir,
                               format      = c("csv", "xlsx"))
mus$top


#--- NeuralCrest --------------------------------------------------------------

#7 clusters (0-6)
nc <- load_subcluster("NeuralCrest", output_dir)
DimPlot(nc$obj, group.by = nc$resolutions$column, label = TRUE)

#resolution not recorded in the chats: pick the one giving 7 clusters
nc <- finalise_subcluster(nc, resolution = 0.6, format = c("csv", "xlsx"))
nc_local <- neighbour_markers(nc$obj,
                              cluster_col = "subcluster",
                              reduction   = "pca",
                              dims        = 1:20,
                              k           = 1,
                              logfc       = 0.2,
                              min_pct     = 0.1,
                              dir         = nc$dir,
                              format      = c("csv", "xlsx"))
nc$top


###----------------------------------------------------------------------------
### 5b. apply the fine labels
###----------------------------------------------------------------------------

#######################################################################
# OPTION A: MANUAL SUBCLUSTER ANNOTATION                              #
#######################################################################
#one named vector per compartment. compartments left out of this list fall
#back to their annotation.csv, so the two routes can be mixed.
#
#naming: Lineage.Type_Suffix. Rem.* = belongs to another compartment and is
#moved to the Rescue pool in 5c.

sub_labels <- list(

  "Epithelial" = c(
    "0"  = "Unk.LowQual",          # Cxcl12, Neat1, Rsrp1
    "1"  = "IFE.Spinous_L",        # Aqp3, Gja1, Dsc3, Klf5, Tdh, Il6st, Stom, Mafb
    "2"  = "KC.Basal",             # Trp63, Col7a1, Itgb4, Lama3, Fgfr2, Dst, Dll1, Hlf
    "3"  = "HF.Isthmus",           # Lrig1, Krt79, Alox12e, Sox9, Fzd7, Fst, Sostdc1
    "4"  = "IFE.Spinous",          # Dsc1, Dsg1a, Gjb2, Gjb6, Klf4, Hal, Csta1
    "5"  = "KC.Basal_ORS",         # Moxd1, Wnt10a, Sox21, Irx2, Dsg2, Lamb3, Itga3
    "6"  = "HF.Bulge",             # Fgf18, Cd34, Dkk3, Nt5e, Cyp26b1, Sox9, Lgr5, Cxcl14
    "7"  = "IFE.Spinous_U",        # Krt2, Tjp3, Cgn, Ptgs1, Krt78, Lsr, Exph5
    "8"  = "KC.Basal_Prol",        # Mki67, Slc2a1, Il33, Ldha, Pgk1, Pthlh, Bnc1, Cdh3
    "9"  = "KC.Basal_Cyc",         # Cdk1, Top2a, Ccna2, Cenpa, Cdc20, Kif20a, Birc5
    "10" = "Rem.Fibro",            # Col1a1, Col6a2, Fn1, Fstl1, Fbn1, Mmp2
    "11" = "SG.Sebocyte",          # Fasn, Acly, Hmgcr, Msmo1, Idi1, Lss, Lpl, Gpam
    "12" = "IFE.Granular",         # Smpd3, Alox12b, Ovol1, Prdm1, Klk5, Klk6, Il1f5
    "13" = "Rem.Immune_LC",        # Ptprc, Cd53, Ciita, Runx3, Aif1, Id2, Selplg
    "14" = "IFE.Cornified",        # Cpa4, Gba2, Aqp5, Aif1l, Gpsm1, Scnn1b, Lipe
    "15" = "KC.Basal_Mig",         # Serpine1, Mmp13, Snai2, Mmp9, Itga5, Pdpn, Inhba, Has2
    "16" = "IFE.Suprabasal_Infl",  # Saa1, Chil1, Serpinb3b, Il1rn, Osmr, Socs3, Cd274
    "17" = "SG.Sebocyte_Mat",      # Plin2, Cd36, Tgm2, Hopx, Cers5, Krt79, Rab27a
    "18" = "IFE.Spinous",          # Hal, Krt78, Dsg1a, Ppl, Evpl, Cldn1 (low count)
    "19" = "HF.Bulb",              # Dct, Tyrp1, Lgr5, Ptch1, Ptch2, Mycn, Lhx2, Tcf7
    "20" = "IFE.Cornified_Act"     # Sprr2f, Arc, Egr1, Maff, Hk2, Osmr, Hcar2
  ),

  "Stromal" = c(
    "0"  = "Rem.Immune_Mac",       # Cd68, Adgre1, C3ar1, Ccr5, Ms4a7, Csf1r (glycolytic)
    "1"  = "Fib.Dpp4",             # Dpp4, Sfrp2, Cd34, Sfrp4, Sema3c, Cd55, Tnxb, Ace, C3
    "2"  = "Fib.Activated",        # Lrrc15, Mmp13, Postn, Col7a1, Crabp1, Twist1, Itga11, Tgfb1i1
    "3"  = "Fib.Reticular",        # Fmod, Thbs4, Igf2, Eln, Ogn, Mest, H19, Piezo2, Kcnma1
    "4"  = "Rem.LowQual",          # 29 genes, all pct.1 < pct.2
    "5"  = "Fib.Inflam_Gly",       # Cxcl5, Ereg, Serpine1, Cxcl1, Has2, Slc2a1, Hk2, Ldha
    "6"  = "Rem.LowQual",          # Col12a1, Lrrc15, Tnc (max FC 1.2; papillary-flavoured)
    "7"  = "Fib.Col15a1",          # Col15a1, Col4a1, Hmcn2, Col18a1, Lama2, Abca8a, Cxcl14, Lpl
    "8"  = "Fib.Activated_Cyc",    # Mki67, Birc5, Cdk1, Aurkb, Top2a, Foxm1 (+Lrrc15/Crabp1)
    "9"  = "Adip.Mature",          # Adipoq, Plin1, Retn, Lep, Lipe, Pparg, Srebf1, Mlxipl
    "10" = "Fib.HFDermal",         # Cyp26b1, Ntn1, Nkd2, Igfbp2, Lgr5, Ltbp1, Slit2, Fzd2
    "11" = "Rem.LowQual",          # Neat1, Rsrp1, Dnm1, Hnrnph1, Xist (nuclear-biased)
    "12" = "Fib.Inflam",           # Il33, Il1rl1, Cxcl5, Timp1, Inhba, Ccl2, Fgf7, Vegfa, Pdpn
    "13" = "Rem.Immune_Mac",       # Cd163, Folr2, Clec10a, Ccl8, F13a1, Mrc1, Siglec1
    "14" = "Rem.Endo",             # Pecam1, Cdh5, Egfl7, Plvap, Aplnr, Sox18, Tie1, Kdr
    "15" = "Rem.Mixed_AdipMac"     # Cd163, Clec10a, Adipoq, Folr2, F13a1, Plin1, Lep
  ),

  "Immune" = c(
    "0"  = "Mono.Hypoxic",         # Arg1, Vegfa, Hif1a, Cd14, Itgb2, Egln3, Bnip3, Ldha
    "1"  = "Rem.Fibro",            # Col1a2, Col1a1, Col6a2 (+weak Mrc1/Adgre1 spillover)
    "2"  = "Mac.Dermal",           # Mrc1, Adgre1, Csf1r, Ms4a7, Cd68, Sirpa, Mertk, Trem2
    "3"  = "Mono.Infl",            # Ccl3, Cxcl2, Cd14, Acod1, Il1rn, Cxcl3
    "4"  = "Mono.Infl",            # Nlrp3, Acod1, Cd14, Cd53, Cxcr2, Trem1, C5ar1, Ptgs2
    "5"  = "Mac.Lipid",            # Trem2, Gpnmb, Lpl, Cd36, Plin2, Mmp12, Cd68, Itgb2
    "6"  = "DC.Migratory",         # Flt3, Ccr7, Ciita, Fscn1, Cd209a, Batf3, Naaa, Irf8
    "7"  = "Mac.Resident",         # Cd163, Folr2, Csf1r, Mrc1, C4b, Clec10a, Siglec1, Gas6
    "8"  = "Mac.Undet",            # Mrc1, Pltp, Ms4a7, Grn, Hexa, F13a1
    "9"  = "Neut.Sell",            # Sell, Csf3r, Cd14, Mxd1, Marcksl1
    "10" = "Mac.Perivasc",         # Lyve1, Csf1r, Mrc1, Cd163, C6, Ednrb, Siglec1, Fcrls
    "11" = "Lymph.TNK",            # Cd3e, Lck, Ptprc, Il2rb, Lat, Gzma, Ncr1, Il7r, Tcf7
    "12" = "Mono.Infl",            # Ccl3, Il1rn, Cd14, Thbs1, Marcksl1 (low conf.)
    "13" = "Mono.Classical",       # Chil3, Vcan, Ccr2, Fn1, Sell, Itgam, Fcgr1, Cd244a
    "14" = "Rem.Mural",            # Tagln, Cav1, Cavin1, Cttn, Mprip, Fscn1
    "15" = "Mac.Dermal_IFN",       # Ifit3, Zbp1, Ddx58, Stat1, Stat2, Isg15, Oasl2, Irf1
    "16" = "Neut.Sell_IFN",        # Isg15, Oasl2, Cd14, Sell, Csf3r, Mxd1, Cxcl2
    "17" = "Mac.Lipid_Cyc",        # Mki67, Top2a, Mrc1, Ccna2, Cdk1, Birc5 (+Lpl/Cd36/Gpnmb)
    "18" = "Neut.Mmp9",            # Mmp9, Csf3r, Cd14, Mxd1, Il1rn
    "19" = "Rem.Fibro",            # Col12a1, Lrrc15, Postn, Tnc, Crabp1, Lum, Pdgfrb, Prrx1
    "20" = "Mast.Cell",            # Tpsb2, Kit, Ms4a2, Slc18a2, Il1rl1, Gata2, Tpsab1
    "21" = "Mac.Osteoclastic"      # Acp5, Ctsk, Csf1r, Nfatc1, Ocstamp, Atp6v0d2, Itgb3
  ),

  "Endothelial" = c(
    "0"  = "Endo.Capillary",       # Cd36, Tek, Flt1, Esam, Ephb4, Notch4, Tjp1
    "1"  = "Endo.Angiogenic",      # Apln, Igfbp3, Cd276, Nid2, Lama4, Mest, Fn1, Mmp14, Cxcr4
    "2"  = "Rem.Fibro",            # Col1a1, Lum, Pdgfra, Col14a1, Postn, Sfrp2, Smoc2, Ogn
    "3"  = "Endo.Venous",          # Ackr1, Selp, Sele, Vwf, Nr2f2, Icam1, Jam2, Lifr
    "4"  = "Endo.Angiogenic_Cyc",  # Mki67, Cdk1, Birc5, Aurkb, Top2a, Foxm1 (+Apln/Angpt2)
    "5"  = "Rem.Immune_Mac",       # Cd163, Mrc1, Csf1r, Adgre1, Cd68, Itgb2, Cd53, F13a1
    "6"  = "Endo.Arterial",        # Gja5, Gja4, Hey1, Efnb2, Dll4, Sox17, Jag1 (+Myh11 spillover)
    "7"  = "Endo.Lymphatic"        # Prox1, Lyve1, Ccl21a, Ccl21b, Mmrn1, Flt4, Reln
  ),

  "Mural_Muscle" = c(
    "0"  = "Mural.Myofibro",       # Lrrc15, Tnc, Col12a1, Notch3, Pdgfrb, Mcam, Ednra, Cd248, Tgfb1i1
    "1"  = "Mural.SMC",            # Myh11, Tagln, Smtn, Synpo2, Actg2, Gucy1a1 (low depth)
    "2"  = "Rem.Endo",             # Pecam1, Cdh5, Vwf, Tie1, Kdr, Plvap, Sox18, Clec14a, Tek
    "3"  = "Musc.Myofibre",        # Myh1, Ttn, Ryr1, Actn2, Cacna1s, Casq1, Trdn, Des (low depth)
    "4"  = "Musc.Myofibre",        # Myh4, Myh2, Myf6, Klhl41, Igfn1, Six1, Trim63, Pdk4, Sgca
    "5"  = "Mural.Pericyte_Cyc",   # Mki67, Cdk1, Aurka, Cspg4, Cd248, Mcam, Col15a1
    "6"  = "Musc.Regenerating",    # Mymk, Myog, Myod1, Myh3, Myh8, Chrng, Chrnd, Chrna1, Musk
    "7"  = "Rem.Immune_Mac",       # Cd163, Mrc1, Csf1r, Adgre1, Cd68, F13a1, Trem2, Fcrls
    "8"  = "Rem.Fibro",            # Thbs4, Col14a1, Ogn, Eln, Mfap5, Tnxb, Lum, Smoc2, Hmcn2
    "9"  = "Mural.Pericyte"        # Kcnj8, Abcc9, Vtn, Pdgfrb, Notch3, Sept4, Agtr1a, Gipr, Il34
  ),

  "NeuralCrest" = c(
    "0"  = "Rem.LowQual",          # Gfra1, Gap43, Ngfr (14 genes, all pct.1 < pct.2)
    "1"  = "Rem.Fibro",            # Col12a1, Postn, Col6a1, Lrrc15, Spon1, Mest, Mki67
    "2"  = "NC.Schwann_Rep",       # Ngfr, Nes, Gap43, Zeb2, Cadm1, Nrcam, Met, Mcam (provisional)
    "3"  = "Rem.Fibro",            # Col1a2, Mmp2, Col6a2, Slit3, Col1a1, Fbn1, Thbs2
    "4"  = "NC.Schwann_Myel",      # Dusp15, Mag, Ugt8a, Dhh, Cnp, S100b, Cntf, Gpr37l1, Sema6c
    "5"  = "Nerve.Perineurial",    # Slc2a1, Cldn1, Gjb2, Foxc2, Cxadr, Klf5, Itgb4 (provisional)
    "6"  = "Rem.Immune_Mac"        # Cxcl16, Fcgr3, Ly86, Il10ra, Ccr5, Syk, Ms4a7, Itgb2
  )
)

obj <- apply_subcluster_labels(obj, labels = sub_labels,
                              broad_col  = "broad_type",
                              label_col  = "cell_type",
                              output_dir = output_dir,
                              skip_types = NULL,
                              prefix     = TRUE)

#######################################################################
# OPTION B: LABEL IN EACH _subclusters/<compartment>/annotation.csv    #
#######################################################################
# as in section 3, option B overwrites cell_type, so leave it commented
# out while option A above is in use

# obj <- apply_subcluster_labels(obj,
#                               broad_col  = "broad_type",
#                               label_col  = "cell_type",
#                               output_dir = output_dir,
#                               skip_types = NULL,
#                               prefix     = TRUE)


#######################################################################
#######################################################################
#                                                                     #
#                     5c. RECLUSTER THE Rem. CLUSTERS                 #
#                                                                     #
#######################################################################
#######################################################################

#move every Rem.* cell out of its compartment into one Rescue pool.
#include_unlabelled also collects cells whose cell_type is NA or still the
#bare compartment name
obj <- stage_rescue_compartment(obj,
                                label_col          = "cell_type",
                                broad_col          = "broad_type",
                                prefix             = "Rem.",
                                rescue_as          = "Rescue",
                                include_unlabelled = TRUE)

#overwrite = FALSE, so only the new Rescue compartment gets built
subclusters <- run_subclustering(
  obj,
  broad_col    = "broad_type",
  output_dir   = output_dir,
  skip_types   = NULL,
  min_cells    = 500,
  dims         = 1:20,
  resolutions  = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
  resolution   = 0.6,
  n_features   = 2000,
  sketch_cells = 20000,
  direct_max   = 75000,
  markers      = FALSE,
  overwrite    = FALSE
)

resc <- load_subcluster("Rescue", output_dir)

DimPlot(resc$obj, group.by = "RNA_snn_res.0.6", label = TRUE)

resc_markers <- subcluster_markers(resc$obj,
                                   cluster_col = "RNA_snn_res.0.6",
                                   only_pos    = TRUE,
                                   logfc       = 0.2,
                                   min_pct     = 0.2,
                                   dir         = resc$dir,
                                   format      = c("csv", "xlsx"))

local_markers <- neighbour_markers(resc$obj,
                                   cluster_col = "RNA_snn_res.0.6",
                                   reduction   = "pca",
                                   dims        = 1:20,
                                   k           = 1,
                                   logfc       = 0.2,
                                   min_pct     = 0.1,
                                   dir         = resc$dir,
                                   format      = c("csv", "xlsx"))

resc <- finalise_subcluster(resc, resolution = 0.6)

#add the Rescue labels without touching the sub_labels list above
sub_labels[["Rescue"]] <- c(
  "0"  = "Fib.Activated",          # Col12a1, Tnc, Postn      (+Lrrc15, Thbs2, Fzd1, Timp1, Itga5)
  "1"  = "Mac.Lipid",              # Trem2, Gpnmb, Plin2      (+Cd68, C3ar1, Ms4a7, Cd36, Mrc1, Csf1r)
  "2"  = "Rem.LowQual",            # Col1a2, Col1a1, Fbn1     (21 genes, 19 with pct.1 < pct.2)
  "3"  = "Fib.Reticular",          # Col14a1, Ogn, Tnxb       (+Sfrp2, Smoc2, Eln, Lum, Mfap5, C4b)
  "4"  = "Endo.Vascular",          # Pecam1, Cdh5, Egfl7      (+Plvap, Tie1, Sox18, Kdr, Aplnr, Cd93)
  "5"  = "Rem.LowQual",            # Cxcl12, Neat1, Rsrp1     (3 genes, all pct.1 < pct.2)
  "6"  = "Rem.LowQual",            # Neat1, App, Cxcl12       (5 genes; Dsp/Trim29 vs c5, KC-flavoured)
  "7"  = "Mac.Hypoxic",            # Arg1, Cd14, Hif1a        (+Hmox1, Slc2a1, Ldha, Itgam)
  "8"  = "Imm.DCLymph",            # Flt3, Ciita, Irf8        (+Batf3, Il7r, Klrk1, Itgb7; merged)
  "9"  = "Adip.Mature",            # Adipoq, Plin1, Retn      (+Lep, Gpd1, Lipe, Cebpa, Pnpla2)
  "10" = "Mac.Resident",           # Cd163, Folr2, Clec10a    (+Ccl8, F13a1, Siglec1, Ednrb, C6)
  "11" = "Mural.PericyteSMC",      # Notch3, Kcnj8, Abcc9     (+Myh11, Mcam, Pdgfrb, Il34, Gucy1a1)
  "12" = "Cyc.Unassigned",         # Mki67, Cdk1, Top2a       (+Cdca8, Birc5; no lineage)
  "13" = "NC.Schwann",             # Sox10, Ngfr, Erbb3       (+Gldn, Gap43, Gfra1, L1cam, Cadm1)
  "14" = "Musc.Myofibre",          # Ttn, Actn2, Ryr1         (+Myh1, Myh2, Myh4, Des, Cacna1s)
  "15" = "KC.Epidermal",           # Dsp, Dsc3, Dsg1a         (+Cdh1, Epcam, Trp63, Aqp3, Csta1)
  "16" = "Rem.LowQual"             # Cxcl12, Rsrp1            (2 genes, both pct.1 < pct.2)
)

#prefix takes the compartments to prefix, so Rescue labels stay unqualified
obj <- apply_subcluster_labels(obj, labels = sub_labels,
                              broad_col  = "broad_type",
                              label_col  = "cell_type",
                              output_dir = output_dir,
                              skip_types = NULL,
                              prefix     = setdiff(names(sub_labels), "Rescue"))

#every cell should now carry a fine label: both should be 0
sum(is.na(obj$cell_type))
sum(obj$cell_type == obj$broad_type, na.rm = TRUE)


#######################################################################
#######################################################################
#                                                                     #
#                    5d. COLOURS AND FULL OBJECT                      #
#                                                                     #
#######################################################################
#######################################################################

colours <- build_colour_map(obj, label_col = "cell_type",
                            output_dir = output_dir)

#image names are sanitised by Seurat, so the fov is dotted not hyphenated
ImageDimPlot(obj, fov = "XE765.L", group.by = "cell_type",
             cols = colours, size = 1)

plot_cohort_spatial(obj, group_col = "cell_type", colours = colours)

DefaultAssay(obj) <- "RNA"

#full cohort, both formats, before the per region split
saveRDS(obj, file.path(sketch_dir, "merged_annotated.rds"))

full <- obj
full[["sketch"]] <- NULL
full[["RNA"]]    <- JoinLayers(full[["RNA"]])
export_h5ad(full, file.path(sketch_dir, "merged_annotated.h5ad"))
rm(full); gc(verbose = FALSE)


#######################################################################
#######################################################################
#                                                                     #
#                    6. EXPORT PER REGION OUTPUTS                     #
#                                                                     #
#######################################################################
#######################################################################

DefaultAssay(obj) <- "RNA"

export_annotated_regions(
  obj,
  output_dir = output_dir,
  label_col  = "cell_type",
  colours    = colours,
  write_rds  = TRUE,
  write_h5ad = TRUE,
  plots      = TRUE
)
