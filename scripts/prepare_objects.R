source("./scripts/setup_py_environment.R")
setup_py_env("hyms_metal_diabetes", "/home/uvictor/miniconda3/condabin/conda")


source("./scripts/helper_functions.R")
source("./scripts/convert_to_h5ad.R")

reticulate::py_run_string("
import sys
sys.path.append('./scripts')
")

zconv <- reticulate::import("zarr_h5ad_conversion", convert = TRUE)

build_proseg_seurat(
  proseg_dir = "./proseg_results_mask",
  xenium_dir = "/mnt/d/HYMS/metal_diabetes",
  overwrite = TRUE
)

annotate_proseg_seurat_10x(
  proseg_dir = "./proseg_results_mask",
  reference_dir = "./reference",
  ref_label_col = "subclustering",
  output_dir = "./annotated_data",
  ref_method = "pca",
  dims = 1:30,
  k.weight = 20,
  sketch_ncells = 45000,
  plots = TRUE,
  overwrite = TRUE,
  replot_only = FALSE
)

# annotate_proseg_seurat(
#     proseg_dir         = "./proseg_results",
#     reference_dir      = "./reference",
#     ref_label_col      = "CellType",
#     annotation_method  = "both",
#     cluster_resolution = 5,
#     rctd_mode          = "doublet", 
#     rctd_max_cores     = 32,
#     singler_delta_floor = 0.05,
#     singler_delta_pctile_threshold = 0.5,  # was 0.95 — much too strict for singler-only
#     overwrite          = TRUE
# )

conversion_log <- convert_annotated_rds(
  base_dir     = "./annotated_data",
  main_layer   = "counts",
  overwrite = TRUE,
  other_layers = NULL
)


print(conversion_log)

