# One-command TNBC1 -> Cottrazm BRCA1 reference workflow.
#
# The core script builds the reference. This wrapper orders markers by
# cell-type specificity and pushes mitochondrial/ribosomal/ubiquitous genes
# behind biological markers because Cottrazm prioritizes the first 25 genes.

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
technical_pattern <- "^(MT-|RPL|RPS)|^(MALAT1|NEAT1)$"

for (cell_type in colnames(sig_exp)) {
  genes <- intersect(clustermarkers_list[[cell_type]], rownames(sig_exp))
  other_types <- setdiff(colnames(sig_exp), cell_type)
  specificity <- sig_exp[genes, cell_type] -
    rowMeans(sig_exp[genes, other_types, drop = FALSE])
  ordered_genes <- genes[order(specificity, decreasing = TRUE, na.last = NA)]
  is_technical <- grepl(technical_pattern, ordered_genes, ignore.case = TRUE)
  clustermarkers_list[[cell_type]] <- c(
    ordered_genes[!is_technical],
    ordered_genes[is_technical]
  )
}
saveRDS(clustermarkers_list, marker_file, compress = "gzip")

marker_check <- readRDS(marker_file)
stopifnot(
  identical(colnames(sig_exp), names(marker_check)),
  all(lengths(marker_check) > 0L),
  all(vapply(
    marker_check,
    function(x) !any(grepl(technical_pattern, head(x, 25), ignore.case = TRUE)),
    FUN.VALUE = logical(1)
  ))
)
message("Marker lists reordered; technical genes excluded from each top 25.")
