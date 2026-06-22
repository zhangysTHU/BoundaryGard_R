# Build the BRCA1 Cottrazm reference from the four TNBC (_ca) samples.
#
# Inputs:
#   TNBC1/sample_aggr_combined.rds
#   TNBC1/clustermarkers_list.csv
#
# Outputs:
#   ../input/BRCA1/single_cell/sig_exp.rds.gz
#   ../input/BRCA1/single_cell/clustermarkers_list.rds.gz
#   ../input/BRCA1/single_cell/reference_qc.tsv
#   ../input/BRCA1/single_cell/reference_cell_metadata.tsv.gz
#
# Cell identities are assigned strictly by the CSV Cluster -> Celltype mapping.
# No malignant-cell label is inferred. Normal (_n) samples are excluded.

options(stringsAsFactors = FALSE)
set.seed(666)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0L) {
  stop("Run this file with Rscript.", call. = FALSE)
}
script_file <- normalizePath(
  sub("^--file=", "", script_arg[[1]]),
  winslash = "/",
  mustWork = TRUE
)
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
technical_pattern <- "^(MT-|RPL|RPS)|^(MALAT1|NEAT1)$"

required_packages <- c("Seurat", "Matrix")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
if (!file.exists(input_rds) || !file.exists(annotation_csv)) {
  stop("TNBC1 RDS or annotation CSV is missing.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation <- utils::read.csv(
  annotation_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
required_columns <- c("Cluster", "Celltype", "Marker")
if (!all(required_columns %in% colnames(annotation))) {
  stop("Annotation CSV must contain Cluster, Celltype and Marker.", call. = FALSE)
}
annotation$Cluster <- as.character(annotation$Cluster)
annotation$Celltype <- trimws(annotation$Celltype)
if (anyDuplicated(annotation$Cluster)) {
  stop("Annotation CSV contains duplicated cluster IDs.", call. = FALSE)
}

message("Reading large Seurat object: ", input_rds)
sc_obj <- readRDS(input_rds)
if (!inherits(sc_obj, "Seurat") || !assay_name %in% names(sc_obj@assays)) {
  stop("Input must be a Seurat object with an RNA assay.", call. = FALSE)
}
if (!all(c(cluster_column, sample_column) %in% colnames(sc_obj@meta.data))) {
  stop("Required metadata columns are missing.", call. = FALSE)
}

# Final Cottrazm reference uses only the four cancer-tissue samples.
all_sample_id <- as.character(sc_obj@meta.data[[sample_column]])
ca_cells <- grepl("_ca$", all_sample_id, ignore.case = TRUE)
if (!any(ca_cells)) {
  stop("No sample IDs ending in _ca were found.", call. = FALSE)
}
sc_obj <- subset(sc_obj, cells = colnames(sc_obj)[ca_cells])
sample_id <- as.character(sc_obj@meta.data[[sample_column]])
ca_samples <- sort(unique(sample_id))
if (length(ca_samples) != 4L) {
  warning("Expected four _ca samples but found ", length(ca_samples), ": ",
          paste(ca_samples, collapse = ", "))
}
message("Included _ca samples: ", paste(ca_samples, collapse = ", "))
message("Excluded _n cells: ", sum(!ca_cells))

cluster_id <- as.character(sc_obj@meta.data[[cluster_column]])
cluster_to_type <- stats::setNames(annotation$Celltype, annotation$Cluster)
missing_clusters <- setdiff(unique(cluster_id), names(cluster_to_type))
if (length(missing_clusters) > 0L) {
  stop("Clusters absent from annotation CSV: ",
       paste(sort(missing_clusters), collapse = ", "), call. = FALSE)
}
csv_celltype <- unname(cluster_to_type[cluster_id])

# This changes spelling only, not biological grouping. Epithelial and
# fibroblast names intentionally match 00_config.R.
display_name_map <- c(
  Epithelial_cells = "Epithelial cells",
  Basal_myoepithelial_cells = "Basal myoepithelial cells",
  Fibroblasts = "Fibroblast cells",
  Macrophages = "Macrophage cells",
  T_NK_cells = "T/NK cells",
  B_cells = "B cells",
  Endothelial_cells = "Endothelial cells",
  Plasma_cells = "Plasma cells",
  Mast_cells = "Mast cells"
)
unmapped <- setdiff(unique(csv_celltype), names(display_name_map))
if (length(unmapped) > 0L) {
  stop("CSV cell types lack a display-name mapping: ",
       paste(sort(unmapped), collapse = ", "), call. = FALSE)
}
final_type <- unname(display_name_map[csv_celltype])
type_order <- unname(display_name_map[
  names(display_name_map) %in% unique(csv_celltype)
])

sc_obj$cottrazm_celltype <- factor(final_type, levels = type_order)
Seurat::Idents(sc_obj) <- "cottrazm_celltype"
Seurat::DefaultAssay(sc_obj) <- assay_name
cell_counts <- table(sc_obj$cottrazm_celltype)
message("Cells per Cottrazm type:")
print(cell_counts)

get_assay_data <- function(object, assay, layer) {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}
norm_exp <- get_assay_data(sc_obj, assay_name, "data")
if (nrow(norm_exp) == 0L || ncol(norm_exp) == 0L) {
  sc_obj <- Seurat::NormalizeData(sc_obj, assay = assay_name, verbose = FALSE)
  norm_exp <- get_assay_data(sc_obj, assay_name, "data")
}

# Cottrazm needs ranked marker identities, not marker-test p-values. Rank by
# group-vs-rest mean log-expression difference and within-group detection.
message("Ranking cell-type markers from sparse RNA expression...")
marker_cells <- unlist(
  lapply(type_order, function(cell_type) {
    cells <- which(sc_obj$cottrazm_celltype == cell_type)
    if (length(cells) > marker_max_cells_per_type) {
      sample(cells, marker_max_cells_per_type)
    } else {
      cells
    }
  }),
  use.names = FALSE
)
marker_exp <- norm_exp[, marker_cells, drop = FALSE]
marker_labels <- droplevels(sc_obj$cottrazm_celltype[marker_cells])
de_marker_list <- setNames(vector("list", length(type_order)), type_order)
for (cell_type in type_order) {
  in_group <- which(marker_labels == cell_type)
  out_group <- which(marker_labels != cell_type)
  mean_in <- Matrix::rowMeans(marker_exp[, in_group, drop = FALSE])
  mean_out <- Matrix::rowMeans(marker_exp[, out_group, drop = FALSE])
  pct_in <- Matrix::rowMeans(marker_exp[, in_group, drop = FALSE] > 0)
  specificity <- mean_in - mean_out
  keep <- which(
    is.finite(specificity) &
      specificity >= marker_logfc_threshold &
      pct_in >= marker_min_pct
  )
  ranked <- keep[order(specificity[keep], pct_in[keep], decreasing = TRUE)]
  genes <- rownames(marker_exp)[ranked]
  is_technical <- grepl(technical_pattern, genes, ignore.case = TRUE)
  de_marker_list[[cell_type]] <- c(genes[!is_technical], genes[is_technical])
}

# Retain the user-supplied annotation genes, but append them after the
# data-driven markers so shared annotation genes do not dominate the top 25.
split_marker_string <- function(x) {
  genes <- trimws(unlist(strsplit(ifelse(is.na(x), "", x), ",", fixed = TRUE)))
  unique(genes[nzchar(genes)])
}
annotation_gene_list <- lapply(annotation$Marker, split_marker_string)
names(annotation_gene_list) <- annotation$Cluster

clustermarkers_list <- setNames(vector("list", length(type_order)), type_order)
for (csv_type in names(display_name_map)) {
  cell_type <- unname(display_name_map[[csv_type]])
  if (!cell_type %in% type_order) {
    next
  }
  source_clusters <- annotation$Cluster[annotation$Celltype == csv_type]
  supplied <- unique(unlist(annotation_gene_list[source_clusters], use.names = FALSE))
  supplied <- intersect(supplied, rownames(norm_exp))
  clustermarkers_list[[cell_type]] <- unique(c(de_marker_list[[cell_type]], supplied))
}
if (any(lengths(clustermarkers_list) == 0L)) {
  stop("At least one reference type has no marker genes.", call. = FALSE)
}

signature_genes <- intersect(
  unique(unlist(clustermarkers_list, use.names = FALSE)),
  rownames(norm_exp)
)
clustermarkers_list <- lapply(
  clustermarkers_list,
  function(x) intersect(unique(x), signature_genes)
)

# Match Cottrazm::get_sig_exp(): 2^(RNA@data)-1, then mean by cell type.
# Transform only sparse nonzero entries to avoid allocating a multi-GB dense
# gene-by-cell matrix.
norm_signature <- norm_exp[signature_genes, , drop = FALSE]
linear_signature <- norm_signature
linear_signature@x <- 2^linear_signature@x - 1
sig_exp <- vapply(
  type_order,
  function(cell_type) {
    cells <- which(sc_obj$cottrazm_celltype == cell_type)
    Matrix::rowMeans(linear_signature[, cells, drop = FALSE])
  },
  FUN.VALUE = numeric(length(signature_genes))
)
rownames(sig_exp) <- signature_genes
colnames(sig_exp) <- type_order

if (
  !identical(colnames(sig_exp), names(clustermarkers_list)) ||
    anyDuplicated(rownames(sig_exp)) ||
    any(!is.finite(sig_exp)) ||
    any(sig_exp < 0)
) {
  stop("Generated Cottrazm reference failed validation.", call. = FALSE)
}

sig_file <- file.path(output_dir, "sig_exp.rds.gz")
marker_file <- file.path(output_dir, "clustermarkers_list.rds.gz")
qc_file <- file.path(output_dir, "reference_qc.tsv")
metadata_file <- file.path(output_dir, "reference_cell_metadata.tsv.gz")
saveRDS(sig_exp, sig_file, compress = "gzip")
saveRDS(clustermarkers_list, marker_file, compress = "gzip")

qc <- data.frame(
  cell_type = type_order,
  n_cells = as.integer(cell_counts[type_order]),
  n_markers = lengths(clustermarkers_list[type_order]),
  marker_coverage_in_sig_exp = 1,
  stringsAsFactors = FALSE
)
utils::write.table(qc, qc_file, sep = "\t", quote = FALSE, row.names = FALSE)

cell_metadata <- data.frame(
  cell_barcode = colnames(sc_obj),
  sample_id = sample_id,
  seurat_cluster = cluster_id,
  csv_celltype = csv_celltype,
  cottrazm_celltype = final_type,
  stringsAsFactors = FALSE
)
metadata_connection <- gzfile(metadata_file, open = "wt")
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
  !"Malignant epithelial cells" %in% colnames(sig_check),
  all(lengths(marker_check) > 0L)
)
message("Completed: ", sig_file)
message("Completed: ", marker_file)
