# 03：运行 R inferCNV 并生成统一格式 CNV calls。
# 输入：
# - intermediate/02_TumorST_clustered.rds.gz：含 Spatial counts、seurat_clusters、NormalScore。
# - intermediate/InferCNV/CellAnnotation.txt：CellID 和 cluster，无表头。
# - resources/gencode_v38_gene_pos.txt：gene、chromosome、start、end。
# 输出：
# - intermediate/03_cnv_calls.tsv：三列 cell_ID、CNVLabel、cnv_score，04 只依赖这个文件。
# - intermediate/03_infercnv_obj.rds.gz：完整 inferCNV 对象，供复核。
# - output/03_infercnv/output_Spatial/：inferCNV checkpoint、HMM observations、dendrogram 等。
# 备注：COTTRAZM_FAST_CNV=1 只生成 smoke-test CNV calls，不代表真实 CNV。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
options(scipen = 100)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "infercnv", "readr", "ape", "dendextend"))

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1) {
    stop(name, " must be an integer >= 1; got ", sQuote(value), call. = FALSE)
  }
  parsed
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) default else value
}

env_flag <- function(name, default = TRUE) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (identical(value, "")) {
    return(default)
  }
  value %in% c("1", "true", "yes", "y")
}

write_cnv_calls <- function(out_dir, TumorST, normal_cluster) {
  calls_file <- file.path(paths$intermediate, "03_cnv_calls.tsv")
  obs_files <- list.files(
    out_dir,
    pattern = "observations\\.txt$",
    full.names = TRUE
  )
  obs_files <- obs_files[!grepl("dendrogram", basename(obs_files), ignore.case = TRUE)]
  if (length(obs_files) < 1) {
    group_file <- list.files(
      out_dir,
      pattern = "cell_groupings$",
      full.names = TRUE
    )
    gene_file <- list.files(
      out_dir,
      pattern = "pred_cnv_genes\\.dat$",
      full.names = TRUE
    )
    if (length(group_file) < 1 || length(gene_file) < 1) {
      log_step("no HMM observations table or HMM grouping/gene files found; skipping ", calls_file)
      return(invisible(FALSE))
    }
    group_file <- group_file[which.max(file.info(group_file)$mtime)]
    gene_file <- gene_file[which.max(file.info(gene_file)$mtime)]
    log_step("building CNV calls from ", basename(group_file), " and ", basename(gene_file))

    cell_groups <- read.delim(group_file, stringsAsFactors = FALSE, check.names = FALSE)
    cnv_genes <- read.delim(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
    group_scores <- tapply(abs(as.numeric(cnv_genes$state) - 3), cnv_genes$cell_group_name, sum, na.rm = TRUE)
    group_scores <- group_scores[!is.na(names(group_scores))]
    cell_groups$cnv_score <- as.numeric(group_scores[cell_groups$cell_group_name])
    cell_groups$cnv_score[is.na(cell_groups$cnv_score)] <- 0

    calls <- data.frame(
      cell_ID = rownames(TumorST@meta.data),
      CNVLabel = "Normal",
      cnv_score = 0,
      stringsAsFactors = FALSE
    )
    normal_cells <- rownames(TumorST@meta.data)[as.character(TumorST$seurat_clusters) == as.character(normal_cluster)]
    observation_groups <- unique(cell_groups$cell_group_name[!cell_groups$cell %in% normal_cells])
    observation_scores <- tapply(
      cell_groups$cnv_score[cell_groups$cell_group_name %in% observation_groups],
      cell_groups$cell_group_name[cell_groups$cell_group_name %in% observation_groups],
      mean
    )
    if (length(observation_scores) > 0) {
      centers <- min(params$cnv_k, length(unique(observation_scores)))
      if (centers <= 1) {
        group_labels <- setNames(rep("1", length(observation_scores)), names(observation_scores))
      } else {
        set.seed(666)
        group_labels <- setNames(as.character(stats::kmeans(as.matrix(observation_scores), centers = centers)$cluster), names(observation_scores))
      }
      cell_groups$CNVLabel <- as.character(group_labels[cell_groups$cell_group_name])
      cell_groups$CNVLabel[is.na(cell_groups$CNVLabel)] <- "Normal"
      matched <- match(calls$cell_ID, cell_groups$cell)
      calls$CNVLabel[!is.na(matched)] <- cell_groups$CNVLabel[matched[!is.na(matched)]]
      calls$cnv_score[!is.na(matched)] <- cell_groups$cnv_score[matched[!is.na(matched)]]
    }
    calls$CNVLabel[calls$cell_ID %in% normal_cells] <- "Normal"
    calls$cnv_score[calls$CNVLabel == "Normal"] <- 0
    write.table(calls, calls_file, sep = "\t", row.names = FALSE, quote = FALSE)
    log_step("wrote ", calls_file)
    return(invisible(TRUE))
  }

  obs_file <- obs_files[which.max(file.info(obs_files)$mtime)]
  log_step("building CNV calls from ", basename(obs_file))
  cnv_table <- read.table(obs_file, header = TRUE, check.names = FALSE)
  cnv_score_table <- abs(as.matrix(cnv_table) - 3)
  cell_scores <- colSums(cnv_score_table)
  names(cell_scores) <- gsub("\\.", "-", names(cell_scores))

  calls <- data.frame(
    cell_ID = rownames(TumorST@meta.data),
    CNVLabel = "Normal",
    cnv_score = 0,
    stringsAsFactors = FALSE
  )

  tree_files <- list.files(
    out_dir,
    pattern = "observations_dendrogram\\.txt$",
    full.names = TRUE
  )
  tree_files <- tree_files[grepl("HMM", basename(tree_files), ignore.case = TRUE)]
  observation_cells <- intersect(calls$cell_ID, names(cell_scores))
  if (length(tree_files) > 0) {
    tree_file <- tree_files[which.max(file.info(tree_files)$mtime)]
    log_step("assigning CNV labels from ", basename(tree_file))
    infercnv_label <- as.data.frame(dendextend::cutree(ape::read.tree(file = tree_file), k = params$cnv_k))
    colnames(infercnv_label) <- "CNVLabel"
    rownames(infercnv_label) <- gsub("\\.", "-", rownames(infercnv_label))
    matched <- match(calls$cell_ID, rownames(infercnv_label))
    calls$CNVLabel[!is.na(matched)] <- as.character(infercnv_label$CNVLabel[matched[!is.na(matched)]])
  } else if (length(observation_cells) > 0) {
    log_step("no dendrogram found; assigning CNV labels by kmeans on HMM scores")
    set.seed(666)
    km <- stats::kmeans(cell_scores[observation_cells], centers = min(params$cnv_k, length(unique(cell_scores[observation_cells]))))
    calls$CNVLabel[match(observation_cells, calls$cell_ID)] <- as.character(km$cluster)
  }

  score_match <- match(calls$cell_ID, names(cell_scores))
  calls$cnv_score[!is.na(score_match)] <- as.numeric(cell_scores[score_match[!is.na(score_match)]])
  calls$CNVLabel[calls$cell_ID %in% rownames(TumorST@meta.data)[as.character(TumorST$seurat_clusters) == as.character(normal_cluster)]] <- "Normal"
  calls$cnv_score[calls$CNVLabel == "Normal"] <- 0
  write.table(calls, calls_file, sep = "\t", row.names = FALSE, quote = FALSE)
  log_step("wrote ", calls_file)
  invisible(TRUE)
}

log_step("loading clustered Seurat object")
TumorST <- readr::read_rds(file.path(paths$intermediate, "02_TumorST_clustered.rds.gz"))
assay <- params$infercnv_assay
out_dir <- file.path(paths$output, "03_infercnv", paste0("output_", assay))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("loaded object: ", nrow(TumorST@meta.data), " spots; assay=", assay)

if (identical(Sys.getenv("COTTRAZM_FAST_CNV"), "1")) {
  log_step("COTTRAZM_FAST_CNV=1 detected; writing smoke-test CNV calls")
  normal_cluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[, c("seurat_clusters", "NormalScore")], TumorST@meta.data$seurat_clusters), function(x) mean(x$NormalScore))), decreasing = TRUE)[1]]
  calls <- data.frame(
    cell_ID = rownames(TumorST@meta.data),
    CNVLabel = as.character(TumorST@meta.data$seurat_clusters),
    cnv_score = as.numeric(TumorST@meta.data$nCount_Spatial)
  )
  calls$CNVLabel[calls$CNVLabel == as.character(normal_cluster)] <- "Normal"
  calls$cnv_score <- calls$cnv_score / max(calls$cnv_score, na.rm = TRUE)
  calls$cnv_score[calls$CNVLabel == "Normal"] <- 0
  write.table(calls, file.path(paths$intermediate, "03_cnv_calls.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  log_step("wrote ", file.path(paths$intermediate, "03_cnv_calls.tsv"))
  quit(save = "no", status = 0)
}

annotation_file <- file.path(paths$intermediate, "InferCNV", "CellAnnotation.txt")
normal_cluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[, c("seurat_clusters", "NormalScore")], TumorST@meta.data$seurat_clusters), function(x) mean(x$NormalScore))), decreasing = TRUE)[1]]
gene_order_file <- file.path(paths$resources, "gencode_v38_gene_pos.txt")
log_step("reference cluster selected: ", normal_cluster)

if (identical(Sys.getenv("COTTRAZM_INFERCNV_POSTPROCESS_ONLY"), "1")) {
  log_step("COTTRAZM_INFERCNV_POSTPROCESS_ONLY=1 detected; writing CNV calls from existing inferCNV outputs")
  ok <- write_cnv_calls(out_dir, TumorST, normal_cluster)
  quit(save = "no", status = if (isTRUE(ok)) 0 else 1)
}

checkpoint <- file.path(out_dir, "04_logtransformedHMMi6.infercnv_obj")
use_checkpoint <- file.exists(checkpoint) && env_flag("COTTRAZM_INFERCNV_USE_CHECKPOINT", TRUE)
if (use_checkpoint) {
  log_step("loading inferCNV checkpoint: ", checkpoint)
  infercnv_obj <- readRDS(checkpoint)
  log_step("checkpoint loaded: ", nrow(infercnv_obj@expr.data), " genes x ", ncol(infercnv_obj@expr.data), " cells")
  gc()
} else {
  log_step("extracting dense count matrix for inferCNV")
  matrix <- as.matrix(Seurat::GetAssayData(TumorST, layer = "counts", assay = assay))
  log_step("matrix ready: ", nrow(matrix), " genes x ", ncol(matrix), " spots")
  log_step("creating inferCNV object")
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = matrix,
    annotations_file = annotation_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = normal_cluster
  )
  log_step("inferCNV object created")
}

infercnv_threads <- env_int("COTTRAZM_INFERCNV_THREADS", min(params$infercnv_threads, 2))
partition_method <- env_chr("COTTRAZM_INFERCNV_PARTITION", "qnorm")
analysis_mode <- env_chr("COTTRAZM_INFERCNV_ANALYSIS_MODE", "subclusters")
Sys.setenv(
  OMP_NUM_THREADS = infercnv_threads,
  MKL_NUM_THREADS = infercnv_threads,
  OPENBLAS_NUM_THREADS = infercnv_threads,
  VECLIB_MAXIMUM_THREADS = infercnv_threads,
  NUMEXPR_NUM_THREADS = infercnv_threads
)
options(mc.cores = infercnv_threads, scipen = 100)

log_step(
  "starting infercnv::run; output dir=", out_dir,
  "; threads=", infercnv_threads,
  "; analysis_mode=", analysis_mode,
  "; partition=", partition_method,
  "; HMM=TRUE; plots=off"
)
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = out_dir,
  cluster_by_groups = FALSE,
  analysis_mode = analysis_mode,
  denoise = TRUE,
  HMM = TRUE,
  tumor_subcluster_partition_method = partition_method,
  HMM_type = "i6",
  BayesMaxPNormal = 0,
  num_threads = infercnv_threads,
  plot_steps = FALSE,
  inspect_subclusters = FALSE,
  resume_mode = TRUE,
  plot_probabilities = FALSE,
  no_prelim_plot = TRUE,
  no_plot = TRUE,
  save_rds = TRUE,
  save_final_rds = TRUE
)
log_step("infercnv::run completed")

log_step("saving inferCNV object")
readr::write_rds(infercnv_obj, file.path(paths$intermediate, "03_infercnv_obj.rds.gz"), compress = "gz")

log_step("building CNV calls from inferCNV outputs")
ok <- write_cnv_calls(out_dir, TumorST, normal_cluster)
if (!isTRUE(ok)) {
  stop("inferCNV finished, but no CNV calls could be reconstructed from output files in ", out_dir, call. = FALSE)
}
