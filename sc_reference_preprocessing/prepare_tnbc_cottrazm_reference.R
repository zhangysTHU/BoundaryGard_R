# Convert an annotated TNBC Seurat object into the two single-cell reference
# files consumed by the Cottrazm R pipeline:
#   sig_exp.rds.gz
#   clustermarkers_list.rds.gz
#
# Usage:
#   Rscript --vanilla prepare_tnbc_cottrazm_reference.R \
#     path/to/tnbc_scRNA.rds \
#     path/to/output/single_cell \
#     cottrazm_celltype \
#     RNA
#
# The label column should contain the final deconvolution labels. To work with
# the defaults in 00_config.R, include these exact labels:
#   Malignant epithelial cells
#   Epithelial cells
#   Fibroblast cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage:",
      "Rscript --vanilla prepare_tnbc_cottrazm_reference.R",
      "<seurat_rds> <output_single_cell_dir> <label_column> [assay]"
    ),
    call. = FALSE
  )
}

input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
label_column <- args[[3]]
assay_name <- if (length(args) >= 4L) args[[4]] else "RNA"

# Match the original Cottrazm vignette defaults.
logfc_threshold <- 0.25
min_pct <- 0.10
only_positive <- TRUE
min_cells_per_type <- 20L
required_pipeline_labels <- c(
  "Malignant epithelial cells",
  "Epithelial cells",
  "Fibroblast cells"
)

required_packages <- c("Seurat", "Matrix")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

message("Reading: ", input_rds)
sc_obj <- readRDS(input_rds)
if (!inherits(sc_obj, "Seurat")) {
  stop("Input object is not a Seurat object.", call. = FALSE)
}
if (!assay_name %in% names(sc_obj@assays)) {
  stop("Assay not found: ", assay_name, call. = FALSE)
}
if (!label_column %in% colnames(sc_obj@meta.data)) {
  stop("Metadata column not found: ", label_column, call. = FALSE)
}

labels <- trimws(as.character(sc_obj@meta.data[[label_column]]))
bad_labels <- is.na(labels) | labels == ""
if (any(bad_labels)) {
  stop(
    sum(bad_labels),
    " cells have missing/empty labels in metadata column ",
    label_column,
    ".",
    call. = FALSE
  )
}
sc_obj@meta.data[[label_column]] <- factor(labels, levels = unique(labels))
Seurat::Idents(sc_obj) <- label_column
Seurat::DefaultAssay(sc_obj) <- assay_name

cell_counts <- sort(table(labels), decreasing = TRUE)
message("Cells per reference type:")
print(cell_counts)
if (any(cell_counts < min_cells_per_type)) {
  warning(
    "Reference types with fewer than ",
    min_cells_per_type,
    " cells: ",
    paste(names(cell_counts)[cell_counts < min_cells_per_type], collapse = ", ")
  )
}

missing_pipeline_labels <- setdiff(required_pipeline_labels, names(cell_counts))
if (length(missing_pipeline_labels) > 0L) {
  warning(
    "These labels used by the defaults in scripts_format_R/00_config.R are absent: ",
    paste(missing_pipeline_labels, collapse = ", "),
    ". Either rename your labels or update params$decon_*_cluster in 00_config.R."
  )
}

# FindAllMarkers expects normalized expression in the assay data layer. If it
# is absent, create it from counts with Seurat's standard LogNormalize method.
get_assay_data <- function(object, assay, layer) {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}

norm_exp <- get_assay_data(sc_obj, assay_name, "data")
if (nrow(norm_exp) == 0L || ncol(norm_exp) == 0L) {
  message("Normalized data layer is empty; running NormalizeData().")
  sc_obj <- Seurat::NormalizeData(
    sc_obj,
    assay = assay_name,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  norm_exp <- get_assay_data(sc_obj, assay_name, "data")
}

message("Finding positive markers for ", length(cell_counts), " reference types...")
markers <- Seurat::FindAllMarkers(
  object = sc_obj,
  assay = assay_name,
  only.pos = only_positive,
  logfc.threshold = logfc_threshold,
  min.pct = min_pct,
  return.thresh = 0.05,
  verbose = FALSE
)
if (nrow(markers) == 0L) {
  stop("FindAllMarkers returned no markers.", call. = FALSE)
}

gene_column <- if ("gene" %in% colnames(markers)) {
  "gene"
} else {
  stop("FindAllMarkers output does not contain a gene column.", call. = FALSE)
}
fc_column <- intersect(
  c("avg_log2FC", "avg_logFC", "avg_diff"),
  colnames(markers)
)
if (length(fc_column) == 0L) {
  stop("Cannot identify a fold-change column in FindAllMarkers output.", call. = FALSE)
}
fc_column <- fc_column[[1]]

# Preserve the identity order used in sig_exp and sort each marker vector so
# Cottrazm's use of the first 25 markers selects the strongest markers.
type_order <- levels(sc_obj@meta.data[[label_column]])
markers$cluster <- as.character(markers$cluster)
markers <- markers[
  !is.na(markers[[gene_column]]) &
    markers[[gene_column]] != "" &
    markers$cluster %in% type_order,
  ,
  drop = FALSE
]

clustermarkers_list <- setNames(vector("list", length(type_order)), type_order)
for (cell_type in type_order) {
  sub_markers <- markers[markers$cluster == cell_type, , drop = FALSE]
  sub_markers <- sub_markers[
    order(sub_markers[[fc_column]], decreasing = TRUE, na.last = NA),
    ,
    drop = FALSE
  ]
  clustermarkers_list[[cell_type]] <- unique(as.character(sub_markers[[gene_column]]))
}

empty_marker_types <- names(clustermarkers_list)[lengths(clustermarkers_list) == 0L]
if (length(empty_marker_types) > 0L) {
  stop(
    "No markers were found for: ",
    paste(empty_marker_types, collapse = ", "),
    call. = FALSE
  )
}

signature_genes <- unique(unlist(clustermarkers_list, use.names = FALSE))
signature_genes <- intersect(signature_genes, rownames(norm_exp))
if (length(signature_genes) == 0L) {
  stop("No marker genes overlap the normalized expression matrix.", call. = FALSE)
}

# Reproduce the original Cottrazm get_sig_exp() transformation exactly:
#   norm_exp <- 2^(RNA@data) - 1
# This is kept for compatibility with the CRC1 reference files.
norm_exp <- norm_exp[signature_genes, , drop = FALSE]
linear_exp <- 2^norm_exp - 1

sig_exp <- vapply(
  type_order,
  function(cell_type) {
    selected_cells <- which(labels == cell_type)
    Matrix::rowMeans(linear_exp[, selected_cells, drop = FALSE])
  },
  FUN.VALUE = numeric(length(signature_genes))
)
rownames(sig_exp) <- signature_genes
colnames(sig_exp) <- type_order

if (!identical(colnames(sig_exp), names(clustermarkers_list))) {
  stop("Internal error: sig_exp columns and marker-list names differ.", call. = FALSE)
}
if (anyDuplicated(rownames(sig_exp))) {
  stop("sig_exp contains duplicated gene names.", call. = FALSE)
}
if (any(!is.finite(sig_exp)) || any(sig_exp < 0)) {
  stop("sig_exp contains invalid values.", call. = FALSE)
}

marker_coverage <- vapply(
  clustermarkers_list,
  function(x) mean(x %in% rownames(sig_exp)),
  FUN.VALUE = numeric(1)
)

sig_file <- file.path(output_dir, "sig_exp.rds.gz")
marker_file <- file.path(output_dir, "clustermarkers_list.rds.gz")
qc_file <- file.path(output_dir, "reference_qc.tsv")

saveRDS(sig_exp, sig_file, compress = "gzip")
saveRDS(clustermarkers_list, marker_file, compress = "gzip")

qc <- data.frame(
  cell_type = type_order,
  n_cells = as.integer(cell_counts[type_order]),
  n_markers = lengths(clustermarkers_list[type_order]),
  marker_coverage_in_sig_exp = marker_coverage[type_order],
  stringsAsFactors = FALSE
)
utils::write.table(
  qc,
  qc_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Round-trip validation of the actual files written to disk.
sig_check <- readRDS(sig_file)
marker_check <- readRDS(marker_file)
stopifnot(
  is.matrix(sig_check),
  is.list(marker_check),
  identical(colnames(sig_check), names(marker_check)),
  identical(rownames(sig_check), signature_genes)
)

message("Cottrazm reference written successfully:")
message("  ", sig_file, " [", nrow(sig_exp), " genes x ", ncol(sig_exp), " cell types]")
message("  ", marker_file)
message("  ", qc_file)
