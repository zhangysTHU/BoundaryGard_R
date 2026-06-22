# Build BRCA1 Cottrazm reference from TNBC1/sample_aggr_combined.rds.
options(stringsAsFactors = FALSE)
set.seed(666)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
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

required_packages <- c("Seurat", "Matrix")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
if (!file.exists(input_rds) || !file.exists(annotation_csv)) {
  stop("TNBC1 input RDS or annotation CSV is missing.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation <- utils::read.csv(
  annotation_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
if (!all(c("Cluster", "Celltype", "Marker") %in% colnames(annotation))) {
  stop("Annotation CSV needs Cluster, Celltype and Marker columns.", call. = FALSE)
}
annotation$Cluster <- as.character(annotation$Cluster)
annotation$Celltype <- trimws(annotation$Celltype)
if (anyDuplicated(annotation$Cluster)) {
  stop("Duplicated Cluster values in annotation CSV.", call. = FALSE)
}

message("Reading large Seurat object: ", input_rds)
sc_obj <- readRDS(input_rds)
if (!inherits(sc_obj, "Seurat") || !assay_name %in% names(sc_obj@assays)) {
  stop("Input is not a Seurat object with an RNA assay.", call. = FALSE)
}
if (!all(c(cluster_column, sample_column) %in% colnames(sc_obj@meta.data))) {
  stop("Required Seurat metadata columns are missing.", call. = FALSE)
}

cluster_id <- as.character(sc_obj@meta.data[[cluster_column]])
sample_id <- as.character(sc_obj@meta.data[[sample_column]])
cluster_to_type <- stats::setNames(annotation$Celltype, annotation$Cluster)
missing_clusters <- setdiff(unique(cluster_id), names(cluster_to_type))
if (length(missing_clusters) > 0L) {
  stop("Unannotated clusters: ", paste(missing_clusters, collapse = ", "), call. = FALSE)
}
broad_type <- unname(cluster_to_type[cluster_id])

is_tumor_sample <- grepl("_ca$", sample_id, ignore.case = TRUE)
is_normal_sample <- grepl("_n$", sample_id, ignore.case = TRUE)
if (any(!(is_tumor_sample | is_normal_sample))) {
  stop("Sample IDs must end in _ca or _n.", call. = FALSE)
}

epithelial_types <- c("Epithelial_cells", "Basal_myoepithelial_cells")
final_type <- broad_type
final_type[broad_type %in% epithelial_types & is_tumor_sample] <-
  "Malignant epithelial cells"
final_type[broad_type %in% epithelial_types & is_normal_sample] <-
  "Epithelial cells"

other_type_map <- c(
  Fibroblasts = "Fibroblast cells",
  Macrophages = "Macrophage cells",
  T_NK_cells = "T/NK cells",
  B_cells = "B cells",
  Plasma_cells = "Plasma cells",
  Endothelial_cells = "Endothelial cells",
  Mast_cells = "Mast cells"
)
for (source_type in names(other_type_map)) {
  final_type[broad_type == source_type] <- unname(other_type_map[[source_type]])
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
unmapped <- setdiff(unique(final_type), preferred_order)
if (length(unmapped) > 0L) {
  stop("Unmapped cell types: ", paste(unmapped, collapse = ", "), call. = FALSE)
}
type_order <- preferred_order[preferred_order %in% unique(final_type)]

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

# Cottrazm needs ranked marker identities, not test p-values. Rank genes by
# sparse group-vs-rest mean difference and within-group detection frequency.
message("Ranking markers from sparse expression...")
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
  log_fc <- mean_in - mean_out
  keep <- which(
    is.finite(log_fc) &
      log_fc >= marker_logfc_threshold &
      pct_in >= marker_min_pct
  )
  ranked <- keep[order(log_fc[keep], pct_in[keep], decreasing = TRUE)]
  de_marker_list[[cell_type]] <- rownames(marker_exp)[ranked]
}

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
    source_type <- names(other_type_map)[other_type_map == cell_type]
    source_clusters <- annotation$Cluster[annotation$Celltype %in% source_type]
  }
  seed_markers[[cell_type]] <- unique(
    unlist(annotation_gene_list[source_clusters], use.names = FALSE)
  )
}

clustermarkers_list <- setNames(vector("list", length(type_order)), type_order)
for (cell_type in type_order) {
  supplied <- intersect(seed_markers[[cell_type]], rownames(norm_exp))
  clustermarkers_list[[cell_type]] <- unique(c(supplied, de_marker_list[[cell_type]]))
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

# Match Cottrazm::get_sig_exp() and the CRC1 reference transformation.
norm_signature <- norm_exp[signature_genes, , drop = FALSE]
linear_signature <- 2^norm_signature - 1
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
  stop("Generated reference failed validation.", call. = FALSE)
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
  sample_group = ifelse(is_tumor_sample, "TNBC", "Normal"),
  seurat_cluster = cluster_id,
  broad_annotation = broad_type,
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
  all(lengths(marker_check) > 0L)
)
message("Completed: ", sig_file)
message("Completed: ", marker_file)
