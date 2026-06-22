# Build the Cottrazm single-cell reference for BRCA1 from TNBC1.
#
# Fixed input:
#   TNBC1/sample_aggr_combined.rds
#   TNBC1/clustermarkers_list.csv
#
# Fixed output:
#   ../input/BRCA1/single_cell/sig_exp.rds.gz
#   ../input/BRCA1/single_cell/clustermarkers_list.rds.gz
#   ../input/BRCA1/single_cell/reference_qc.tsv
#
# Epithelial and basal/myoepithelial cells in *_ca samples are assigned to
# "Malignant epithelial cells"; their counterparts in *_n samples are assigned
# to "Epithelial cells". Other labels come from clustermarkers_list.csv.

options(stringsAsFactors = FALSE)
set.seed(666)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("prepare_tnbc_cottrazm_reference.R", winslash = "/", mustWork = FALSE)
}
script_dir <- dirname(script_file)
scripts_r_dir <- dirname(script_dir)

input_dir <- file.path(script_dir, "TNBC1")
input_rds <- file.path(input_dir, "sample_aggr_combined.rds")
annotation_csv <- file.path(input_dir, "clustermarkers_list.csv")
output_dir <- file.path(scripts_r_dir, "input", "BRCA1", "single_cell")

assay_name <- "RNA"
cluster_column <- "seurat_clusters"
sample_column <- "orig.ident"
marker_logfc_threshold <- 0.25
marker_min_pct <- 0.10
marker_max_cells_per_type <- 5000L
min_cells_per_type <- 20L

required_packages <- c("Seurat", "Matrix")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Use the scripts_format_R renv library.",
    call. = FALSE
  )
}

if (!file.exists(input_rds)) {
  stop("Seurat RDS not found: ", input_rds, call. = FALSE)
}
if (!file.exists(annotation_csv)) {
  stop("Cluster annotation CSV not found: ", annotation_csv, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading cluster annotations: ", annotation_csv)
annotation <- utils::read.csv(
  annotation_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
required_annotation_columns <- c("Cluster", "Celltype", "Marker")
if (!all(required_annotation_columns %in% colnames(annotation))) {
  stop(
    "Annotation CSV must contain: ",
    paste(required_annotation_columns, collapse = ", "),
    call. = FALSE
  )
}
annotation$Cluster <- as.character(annotation$Cluster)
annotation$Celltype <- trimws(annotation$Celltype)
if (anyDuplicated(annotation$Cluster)) {
  stop("Annotation CSV contains duplicated Cluster values.", call. = FALSE)
}

message("Reading large Seurat object: ", input_rds)
sc_obj <- readRDS(input_rds)
if (!inherits(sc_obj, "Seurat")) {
  stop("sample_aggr_combined.rds is not a Seurat object.", call. = FALSE)
}
if (!assay_name %in% names(sc_obj@assays)) {
  stop("RNA assay is absent from the Seurat object.", call. = FALSE)
}
if (!all(c(cluster_column, sample_column) %in% colnames(sc_obj@meta.data))) {
  stop(
    "Required metadata columns are absent: ",
    paste(setdiff(c(cluster_column, sample_column), colnames(sc_obj@meta.data)), collapse = ", "),
    call. = FALSE
  )
}

cluster_id <- as.character(sc_obj@meta.data[[cluster_column]])
sample_id <- as.character(sc_obj@meta.data[[sample_column]])
cluster_to_type <- stats::setNames(annotation$Celltype, annotation$Cluster)
unannotated_clusters <- setdiff(unique(cluster_id), names(cluster_to_type))
if (length(unannotated_clusters) > 0L) {
  stop(
    "Clusters absent from clustermarkers_list.csv: ",
    paste(sort(unannotated_clusters), collapse = ", "),
    call. = FALSE
  )
}

broad_type <- unname(cluster_to_type[cluster_id])
is_tumor_sample <- grepl("_ca$", sample_id, ignore.case = TRUE)
is_normal_sample <- grepl("_n$", sample_id, ignore.case = TRUE)
unknown_sample_group <- !(is_tumor_sample | is_normal_sample)
if (any(unknown_sample_group)) {
  stop(
    "Cannot infer tumor/normal status for samples: ",
    paste(sort(unique(sample_id[unknown_sample_group])), collapse = ", "),
    ". Expected sample IDs ending in _ca or _n.",
    call. = FALSE
  )
}

epithelial_types <- c("Epithelial_cells", "Basal_myoepithelial_cells")
final_type <- broad_type
final_type[broad_type %in% epithelial_types & is_tumor_sample] <- "Malignant epithelial cells"
final_type[broad_type %in% epithelial_types & is_normal_sample] <- "Epithelial cells"

other_type_map <- c(
  Fibroblasts = "Fibroblast cells",
  Macrophages = "Macrophage cells",
  T_NK_cells = "T/NK cells",
  B_cells = "B cells",
  Endothelial_cells = "Endothelial cells",
  Plasma_cells = "Plasma cells",
  Mast_cells = "Mast cells"
)
for (source_name in names(other_type_map)) {
  final_type[broad_type == source_name] <- unname(other_type_map[[source_name]])
}

preferred_order <- c(
  "Malignant epithelial cells",
  "Epithelial cells",
  "Fibroblast cells",
  "Macrophage cells",
  "T/NK cells",
  "B cells",
  "Plasma cells",
  "Endothelial cells",
  "Mast cells"
)
unmapped_types <- setdiff(unique(final_type), preferred_order)
if (length(unmapped_types) > 0L) {
  stop(
    "Cell types lack a Cottrazm mapping: ",
    paste(sort(unmapped_types), collapse = ", "),
    call. = FALSE
  )
}

type_order <- preferred_order[preferred_order %in% unique(final_type)]
sc_obj$cottrazm_celltype <- factor(final_type, levels = type_order)
Seurat::Idents(sc_obj) <- "cottrazm_celltype"
Seurat::DefaultAssay(sc_obj) <- assay_name

cell_counts <- table(sc_obj$cottrazm_celltype)
message("Cells per Cottrazm reference type:")
print(cell_counts)
if (any(cell_counts < min_cells_per_type)) {
  warning(
    "Reference types with fewer than ",
    min_cells_per_type,
    " cells: ",
    paste(names(cell_counts)[cell_counts < min_cells_per_type], collapse = ", ")
  )
}

get_assay_data <- function(object, assay, layer) {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}

norm_exp <- get_assay_data(sc_obj, assay_name, "data")
if (nrow(norm_exp) == 0L || ncol(norm_exp) == 0L) {
  message("RNA data layer is empty; running NormalizeData().")
  sc_obj <- Seurat::NormalizeData(
    sc_obj,
    assay = assay_name,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  norm_exp <- get_assay_data(sc_obj, assay_name, "data")
}

message(
  "Finding positive markers using at most ",
  marker_max_cells_per_type,
  " cells per reference type..."
)
markers <- Seurat::FindAllMarkers(
  object = sc_obj,
  assay = assay_name,
  only.pos = TRUE,
  logfc.threshold = marker_logfc_threshold,
  min.pct = marker_min_pct,
  return.thresh = 0.05,
  max.cells.per.ident = marker_max_cells_per_type,
  random.seed = 666,
  verbose = FALSE
)
if (nrow(markers) == 0L) {
  stop("FindAllMarkers returned no marker genes.", call. = FALSE)
}

fc_column <- intersect(c("avg_log2FC", "avg_logFC", "avg_diff"), colnames(markers))
if (!"gene" %in% colnames(markers) || length(fc_column) == 0L) {
  stop("Unexpected FindAllMarkers output columns.", call. = FALSE)
}
fc_column <- fc_column[[1]]
markers$cluster <- as.character(markers$cluster)

split_marker_string <- function(x) {
  genes <- trimws(unlist(strsplit(ifelse(is.na(x), "", x), ",", fixed = TRUE)))
  unique(genes[nzchar(genes)])
}
annotation_gene_list <- lapply(annotation$Marker, split_marker_string)
names(annotation_gene_list) <- annotation$Cluster

seed_markers <- setNames(vector("list", length(type_order)), type_order)
for (cell_type in type_order) {
  if (cell_type %in% c("Malignant epithelial cells", "Epithelial cells")) {
    source_clusters <- annotation$Cluster[annotation$Celltype %in% epithelial_types]
  } else {
    reverse_map <- names(other_type_map)[other_type_map == cell_type]
    source_clusters <- annotation$Cluster[annotation$Celltype %in% reverse_map]
  }
  seed_markers[[cell_type]] <- unique(
    unlist(annotation_gene_list[source_clusters], use.names = FALSE)
  )
}

clustermarkers_list <- setNames(vector("list", length(type_order)), type_order)
for (cell_type in type_order) {
  sub_markers <- markers[markers$cluster == cell_type, , drop = FALSE]
  sub_markers <- sub_markers[
    order(sub_markers[[fc_column]], decreasing = TRUE, na.last = NA),
    ,
    drop = FALSE
  ]
  supplied_markers <- intersect(seed_markers[[cell_type]], rownames(norm_exp))
  de_markers <- unique(as.character(sub_markers$gene))
  clustermarkers_list[[cell_type]] <- unique(c(supplied_markers, de_markers))
}

empty_types <- names(clustermarkers_list)[lengths(clustermarkers_list) == 0L]
if (length(empty_types) > 0L) {
  stop("No markers found for: ", paste(empty_types, collapse = ", "), call. = FALSE)
}

signature_genes <- unique(unlist(clustermarkers_list, use.names = FALSE))
signature_genes <- intersect(signature_genes, rownames(norm_exp))
clustermarkers_list <- lapply(
  clustermarkers_list,
  function(x) intersect(unique(x), signature_genes)
)

# Match the transformation in Cottrazm::get_sig_exp() and the CRC1 reference.
norm_exp <- norm_exp[signature_genes, , drop = FALSE]
linear_exp <- 2^norm_exp - 1
sig_exp <- vapply(
  type_order,
  function(cell_type) {
    cells <- which(sc_obj$cottrazm_celltype == cell_type)
    Matrix::rowMeans(linear_exp[, cells, drop = FALSE])
  },
  FUN.VALUE = numeric(length(signature_genes))
)
rownames(sig_exp) <- signature_genes
colnames(sig_exp) <- type_order

if (!identical(colnames(sig_exp), names(clustermarkers_list))) {
  stop("sig_exp columns and clustermarkers_list names differ.", call. = FALSE)
}
if (anyDuplicated(rownames(sig_exp)) || any(!is.finite(sig_exp)) || any(sig_exp < 0)) {
  stop("sig_exp failed matrix validation.", call. = FALSE)
}

sig_file <- file.path(output_dir, "sig_exp.rds.gz")
marker_file <- file.path(output_dir, "clustermarkers_list.rds.gz")
qc_file <- file.path(output_dir, "reference_qc.tsv")
cell_metadata_file <- file.path(output_dir, "reference_cell_metadata.tsv.gz")

saveRDS(sig_exp, sig_file, compress = "gzip")
saveRDS(clustermarkers_list, marker_file, compress = "gzip")

marker_coverage <- vapply(
  clustermarkers_list,
  function(x) mean(x %in% rownames(sig_exp)),
  FUN.VALUE = numeric(1)
)
qc <- data.frame(
  cell_type = type_order,
  n_cells = as.integer(cell_counts[type_order]),
  n_markers = lengths(clustermarkers_list[type_order]),
  marker_coverage_in_sig_exp = marker_coverage[type_order],
  stringsAsFactors = FALSE
)
utils::write.table(qc, qc_file, sep = "\t", quote = FALSE, row.names = FALSE)

cell_metadata <- data.frame(
  cell_barcode = colnames(sc_obj),
  sample_id = sample_id,
  sample_group = ifelse(is_tumor_sample, "TNBC", "Normal"),
  seurat_cluster = cluster_id,
  broad_annotation = broad_type,
  cottrazm_celltype = final_type,
  stringsAsFactors = FALSE
)
metadata_connection <- gzfile(cell_metadata_file, open = "wt")
utils::write.table(
  cell_metadata,
  metadata_connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
close(metadata_connection)

sig_check <- readRDS(sig_file)
marker_check <- readRDS(marker_file)
stopifnot(
  is.matrix(sig_check),
  is.list(marker_check),
  identical(colnames(sig_check), names(marker_check)),
  identical(rownames(sig_check), signature_genes),
  all(lengths(marker_check) > 0L)
)

message("Cottrazm BRCA1 reference completed:")
message("  ", sig_file, " [", nrow(sig_exp), " genes x ", ncol(sig_exp), " cell types]")
message("  ", marker_file)
message("  ", qc_file)
message("  ", cell_metadata_file)
