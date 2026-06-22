# One-command TNBC1 -> Cottrazm BRCA1 reference workflow.
#
# The core script builds cell labels, ranks sparse-expression markers, creates
# sig_exp, and writes QC metadata. This wrapper then reorders each marker list
# by cell-type specificity in sig_exp so Cottrazm's first 25 markers are the
# most discriminative ones. This is particularly important because malignant
# and normal epithelial cells share several supplied annotation genes.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(
  sub("^--file=", "", script_arg[[1]]),
  winslash = "/",
  mustWork = TRUE
)
script_dir <- dirname(script_file)
core_script <- file.path(script_dir, "prepare_tnbc_cottrazm_reference_core.R")
if (!file.exists(core_script)) {
  stop("Core preprocessing script not found: ", core_script, call. = FALSE)
}

source(core_script, local = FALSE)

sig_exp <- readRDS(sig_file)
clustermarkers_list <- readRDS(marker_file)
for (cell_type in colnames(sig_exp)) {
  genes <- intersect(clustermarkers_list[[cell_type]], rownames(sig_exp))
  other_types <- setdiff(colnames(sig_exp), cell_type)
  specificity <- sig_exp[genes, cell_type] -
    rowMeans(sig_exp[genes, other_types, drop = FALSE])
  clustermarkers_list[[cell_type]] <- genes[
    order(specificity, decreasing = TRUE, na.last = NA)
  ]
}
saveRDS(clustermarkers_list, marker_file, compress = "gzip")

marker_check <- readRDS(marker_file)
stopifnot(
  identical(colnames(sig_exp), names(marker_check)),
  all(lengths(marker_check) > 0L)
)
message("Marker lists reordered by cell-type specificity.")
