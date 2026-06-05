setup_py_env <- function(py_env_name, py_location) {
  options(reticulate.conda_binary = py_location)
  
  existing_envs <- system(paste(py_location, "env list"), intern = TRUE)
  env_exists <- any(grepl(py_env_name, existing_envs))
  
  # Define env_path FIRST
  env_path <- paste0("/home/uvictor/miniconda3/envs/", py_env_name)
  
  # THEN set LD_LIBRARY_PATH
  Sys.setenv(LD_LIBRARY_PATH = paste0(
    env_path, "/lib:",
    Sys.getenv("LD_LIBRARY_PATH")
  ))
  
  if (!env_exists) {
    system(paste(
      py_location, "create --yes --name", py_env_name,
      "'python=3.12' pip umap-learn --quiet -c conda-forge"
    ))
  }
  
  reticulate::use_condaenv(env_path, required = TRUE)
  
  if (!env_exists) {
    reticulate::py_install(
      packages = c("torch", "torchvision"),
      pip = TRUE,
      extra_options = c("--index-url", "https://download.pytorch.org/whl/rocm6.4")
    )
    
    reticulate::py_install(
      packages = "git+https://github.com/broadinstitute/CellBender.git@refs/pull/420/head",
      pip = TRUE
    )
    
    reticulate::py_install(
      packages = c("anndata", "scipy", "h5py", "spatialdata"),
      pip = TRUE
    )
  }
  
  reticulate::py_config()
}