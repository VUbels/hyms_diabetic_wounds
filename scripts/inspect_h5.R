#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(hdf5r))

h5 <- hdf5r::H5File$new(
  "/mnt/d/scrna_datasets/hyms_metal_diabetes/Jacob2023-object_All.h5",
  mode = "r"
)

inspect <- function(node, prefix = "", max_depth = 3, depth = 0) {
  if (depth > max_depth) { cat(sprintf("%s... (max depth)\n", prefix)); return() }
  children <- node$ls()
  for (i in seq_len(nrow(children))) {
    nm   <- children$name[i]
    child <- tryCatch(node[[nm]], error = function(e) NULL)
    if (is.null(child)) { cat(sprintf("%s[ERR]  %s — could not open\n", prefix, nm)); next }
    
    attrs <- tryCatch(hdf5r::list.attributes(child), error = function(e) character(0))
    attr_str <- ""
    if (length(attrs) > 0) {
      attr_vals <- sapply(attrs, function(a) {
        val <- tryCatch({
          v <- child$attr_open(a)$read()
          paste(head(v, 5), collapse=",")
        }, error = function(e) "<err>")
        paste0(a, "=", val)
      })
      attr_str <- paste0("  {", paste(attr_vals, collapse = "; "), "}")
    }
    
    if (inherits(child, "H5Group")) {
      cat(sprintf("%s[GROUP] %s%s\n", prefix, nm, attr_str))
      inspect(child, prefix = paste0(prefix, "  "), max_depth = max_depth, depth = depth + 1)
    } else {
      dims <- tryCatch(child$dims, error = function(e) NA)
      dtype <- tryCatch(child$get_type()$to_text(), error = function(e) "unknown")
      cat(sprintf("%s[DSET]  %s  dtype=%s  dims=%s%s\n",
                  prefix, nm, dtype, paste(dims, collapse="x"), attr_str))
      
      # Only try to read values for safe types
      if (length(dims) == 1 && !grepl("H5T_REFERENCE|H5T_VLEN|H5T_COMPOUND", dtype)) {
        tryCatch({
          vals <- if (dims[1] <= 20) child$read() else child$read()[1:5]
          label <- if (dims[1] <= 20) "values" else "first 5"
          cat(sprintf("%s  -> %s: %s\n", prefix, label, paste(head(vals, 20), collapse=", ")))
        }, error = function(e) {
          cat(sprintf("%s  -> <read error: %s>\n", prefix, conditionMessage(e)))
        })
      }
    }
  }
}

cat("=== ROOT ATTRIBUTES ===\n")
tryCatch({
  for (a in hdf5r::list.attributes(h5)) {
    val <- tryCatch(h5$attr_open(a)$read(), error = function(e) "<err>")
    cat(sprintf("  %s = %s\n", a, paste(head(val, 3), collapse=", ")))
  }
}, error = function(e) cat("  <error reading root attrs>\n"))

cat("\n=== FULL TREE ===\n")
inspect(h5)

# Specifically list obs children with their types
cat("\n=== OBS CHILDREN (detailed) ===\n")
obs <- h5[["obs"]]
obs_children <- obs$ls()
for (i in seq_len(nrow(obs_children))) {
  nm <- obs_children$name[i]
  child <- tryCatch(obs[[nm]], error = function(e) NULL)
  if (is.null(child)) { cat(sprintf("  [ERR] %s\n", nm)); next }
  cls <- class(child)[1]
  dtype <- tryCatch(child$get_type()$to_text(), error = function(e) "N/A")
  dims <- tryCatch(child$dims, error = function(e) NA)
  cat(sprintf("  %s : class=%s  dtype=%s  dims=%s\n", nm, cls, dtype, paste(dims, collapse="x")))
}

h5$close_all()
cat("\nDone.\n")