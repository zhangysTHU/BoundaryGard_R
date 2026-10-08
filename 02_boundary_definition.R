# 02: inferCNV, CNV scoring, and Mal/Bdy/nMal boundary definition.

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
options(scipen = 100)
source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "R", "cottrazm_boundary_core.R"))
source(file.path(script_dir, "R", "visium_hex_neighbors.R"))
source(file.path(script_dir, "R", "spatial_plot_core.R"))
load_required_packages(c("Seurat", "infercnv", "readr", "ape", "dendextend", "ggplot2", "ggpubr", "magrittr", "dplyr", "purrr", "tibble", "assertthat", "png", "grid", "jsonlite", "patchwork"))

# Spatial-neighbor and boundary-expansion helpers.
# 空间邻接和边界扩展 helper。
# 输入通常是 Space Ranger 坐标、UMAP 坐标和已标注的 Mal/Bdy/Normal spot 集合。

nbrs <- function(df_j, MalCellIDAdd, CellIDRaw) {
  # 对一组恶性边缘 spot，找出尚未出现在 CellIDRaw 中的空间邻居，作为下一轮候选。
  nbrs_of_Mal <- lapply(MalCellIDAdd, function(id) {
    nbs <- df_j[[id]]
    nbs[!nbs %in% CellIDRaw]
  })
  names(nbrs_of_Mal) <- MalCellIDAdd
  nbrs_of_Mal
}

load_spaceranger_positions <- function(spaceranger_dir, cells) {
  # 兼容 Space Ranger 新旧坐标文件：
  # tissue_positions.csv 有 header，tissue_positions_list.csv 通常无 header。
  # 统一输出 row/col/imagerow/imagecol，并按 Seurat 对象中的 cells 排序。
  spatial_dir <- file.path(spaceranger_dir, "spatial")
  pos_file <- file.path(spatial_dir, "tissue_positions.csv")
  if (!file.exists(pos_file)) {
    pos_file <- file.path(spatial_dir, "tissue_positions_list.csv")
  }
  if (!file.exists(pos_file)) {
    stop("Cannot find tissue_positions.csv or tissue_positions_list.csv in ", spatial_dir, call. = FALSE)
  }

  first_line <- readLines(pos_file, n = 1, warn = FALSE)
  has_header <- grepl("barcode|array_row|pxl_row", first_line, ignore.case = TRUE)
  pos_raw <- utils::read.csv(pos_file, header = has_header, stringsAsFactors = FALSE)
  if (!has_header) {
    colnames(pos_raw) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  }

  required <- c("barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  missing_cols <- setdiff(required, colnames(pos_raw))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", pos_file, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  pos <- data.frame(
    row = pos_raw$array_row,
    col = pos_raw$array_col,
    imagerow = pos_raw$pxl_row_in_fullres,
    imagecol = pos_raw$pxl_col_in_fullres,
    row.names = pos_raw$barcode
  )
  pos <- pos[cells, , drop = FALSE]
  if (any(!stats::complete.cases(pos))) {
    missing_cells <- cells[!stats::complete.cases(pos)]
    stop("Missing spaceranger coordinates for ", length(missing_cells), " cells; first missing: ", missing_cells[[1]], call. = FALSE)
  }
  pos
}

ClusterUpdate <- function(MalCellIDN, position, df_j, UMAPembeddings, NormalCellID, BdyCellID, MalCellID, x,
                          expand_mal_radius = 0.8) {
  # 第 x 轮边界扩展：
  # 对上一轮新 Mal spot 的未标注邻居，比较其到 Mal 中心和 Bdy 中心的 UMAP 距离，
  # 决定标为 Malx 还是 Bdy。
  P <- lapply(MalCellIDN, function(i) {
    ncellID <- df_j[[i]]
    cpos <- UMAPembeddings[i, ]
    npos <- UMAPembeddings[ncellID, , drop = FALSE]
    nMalID <- ncellID[ncellID %in% MalCellID]
    CiMal <- data.frame(t(apply(rbind(npos[nMalID, , drop = FALSE], cpos), 2, mean)))
    rMal <- lapply(c(nMalID, i), function(id) sqrt(sum((rbind(npos, cpos)[id, ] - CiMal)^2))) |> unlist()
    nBdyID <- ncellID[ncellID %in% BdyCellID]
    nbrsID <- ncellID[!ncellID %in% MalCellID & !ncellID %in% BdyCellID & !ncellID %in% NormalCellID]

    if (length(nbrsID) == 0) {
      return(data.frame(cellID = character(), p1 = numeric(), p2 = numeric(), cluster = character()))
    }

    if (length(nBdyID) <= 1) {
      p <- lapply(nbrsID, function(id) sqrt(sum((npos[id, ] - CiMal)^2)))
      d <- data.frame(cellID = nbrsID, p1 = unlist(p), p2 = unlist(p))
      d$cluster <- ifelse(d$p1 <= expand_mal_radius * max(rMal), "Mal", "Bdy")
    } else {
      CiBdy <- data.frame(t(apply(npos[nBdyID, , drop = FALSE], 2, mean)))
      rBdy <- lapply(nBdyID, function(id) sqrt(sum((npos[id, ] - CiBdy)^2))) |> unlist()
      p1 <- lapply(nbrsID, function(id) sqrt(sum((npos[id, ] - CiMal)^2)))
      p2 <- lapply(nbrsID, function(id) sqrt(sum((npos[id, ] - CiBdy)^2)))
      d <- data.frame(cellID = nbrsID, p1 = unlist(p1), p2 = unlist(p2))
      d$cluster <- ifelse(d$p1 <= expand_mal_radius * max(rMal) & d$p2 > max(rBdy), "Mal", "Bdy")
    }
    d
  })

  PDF <- do.call(rbind, lapply(P, data.frame))
  if (is.null(PDF) || nrow(PDF) == 0) {
    return(data.frame(CellID = character(), Location = character()))
  }
  rownames(PDF) <- NULL
  Cluster <- as.data.frame.array(table(PDF$cellID, PDF$cluster))
  if ("Mal" %in% colnames(Cluster)) {
    Cluster$Location <- ifelse(Cluster$Mal > 0, paste("Mal", x, sep = ""), "Bdy")
  } else {
    Cluster$Location <- rep("Bdy", nrow(Cluster))
  }
  data.frame(CellID = rownames(Cluster), Location = as.character(Cluster$Location))
}

# 02A helpers: inferCNV calls.
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

# 02A: run inferCNV or build CNV calls.
log_step("loading clustered Seurat object")
TumorST <- readr::read_rds(file.path(paths$intermediate, "02_TumorST_clustered.rds.gz"))
assay <- params$infercnv_assay
out_dir <- file.path(paths$output, "03_infercnv", paste0("output_", assay))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_step("loaded object: ", nrow(TumorST@meta.data), " spots; assay=", assay)

normal_cluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[, c("seurat_clusters", "NormalScore")], TumorST@meta.data$seurat_clusters), function(x) mean(x$NormalScore))), decreasing = TRUE)[1]]
annotation_file <- file.path(paths$intermediate, "InferCNV", "CellAnnotation.txt")
gene_order_file <- file.path(paths$resources, "gencode_v38_gene_pos.txt")
log_step("reference cluster selected: ", normal_cluster)

if (identical(Sys.getenv("COTTRAZM_FAST_CNV"), "1")) {
  log_step("COTTRAZM_FAST_CNV=1 detected; writing smoke-test CNV calls")
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
} else if (identical(Sys.getenv("COTTRAZM_INFERCNV_POSTPROCESS_ONLY"), "1")) {
  log_step("COTTRAZM_INFERCNV_POSTPROCESS_ONLY=1 detected; writing CNV calls from existing inferCNV outputs")
  ok <- write_cnv_calls(out_dir, TumorST, normal_cluster)
  if (!isTRUE(ok)) {
    stop("No CNV calls could be reconstructed from existing inferCNV outputs in ", out_dir, call. = FALSE)
  }
} else {
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
}

# 02B: score and visualize CNV.
TumorST <- readr::read_rds(file.path(paths$intermediate, "02_TumorST_clustered.rds.gz"))
assay <- params$infercnv_assay
cnv_outdir <- file.path(paths$output, "03_infercnv", paste0("output_", assay))
out_dir <- file.path(paths$output, "04_cnv_score")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tree_file <- file.path(cnv_outdir, "infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt")
calls_file <- file.path(paths$intermediate, "03_cnv_calls.tsv")
if (file.exists(calls_file)) {
  calls <- read.delim(calls_file, stringsAsFactors = FALSE)
  TumorST@meta.data$CNVLabel <- calls$CNVLabel[match(rownames(TumorST@meta.data), calls$cell_ID)]
} else {
  cell_groupings <- ape::read.tree(file = tree_file)
  infercnv_label <- as.data.frame(dendextend::cutree(cell_groupings, k = params$cnv_k))
  colnames(infercnv_label) <- "infercnv.label"
  missing_cells <- rownames(TumorST@meta.data)[!rownames(TumorST@meta.data) %in% rownames(infercnv_label)]
  infercnv_label <- rbind(infercnv_label, data.frame(row.names = missing_cells, infercnv.label = rep("Normal", length(missing_cells))))
  TumorST@meta.data$CNVLabel <- infercnv_label$infercnv.label[match(rownames(TumorST@meta.data), rownames(infercnv_label))]
}

pdf(file.path(out_dir, paste0(sample_name, "_cnv_label.pdf")), width = 7, height = 7)
print(Seurat::SpatialDimPlot(TumorST, group.by = "CNVLabel", cols = cluster_cols) +
  ggplot2::scale_fill_manual(values = cluster_cols))
dev.off()

pdf(file.path(out_dir, paste0(sample_name, "_reduction_cnvlabel.pdf")), width = 7, height = 7)
print(Seurat::DimPlot(TumorST, group.by = "CNVLabel", cols = cluster_cols) +
  ggplot2::scale_fill_manual(values = cluster_cols))
dev.off()

if (exists("calls")) {
  TumorST@meta.data$cnv_score <- calls$cnv_score[match(rownames(TumorST@meta.data), calls$cell_ID)]
} else {
  cnv_table <- read.table(file.path(cnv_outdir, "infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt"), header = TRUE)
  cnv_score_table <- abs(as.matrix(cnv_table) - 3)
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table))
  colnames(cell_scores_CNV) <- "cnv_score"
  rownames(cell_scores_CNV) <- gsub("\\.", "-", rownames(cell_scores_CNV))
  TumorST@meta.data$cnv_score <- cell_scores_CNV$cnv_score[match(rownames(TumorST@meta.data), rownames(cell_scores_CNV))]
}
TumorST@meta.data$cnv_score <- ifelse(TumorST@meta.data$CNVLabel == "Normal", 0, TumorST@meta.data$cnv_score)

pdf(file.path(out_dir, paste0(sample_name, "_cnv_observation_vlnplot.pdf")), width = 6, height = 4)
print(ggplot2::ggplot(TumorST@meta.data[, c("CNVLabel", "cnv_score")], ggplot2::aes(x = CNVLabel, y = cnv_score, fill = CNVLabel)) +
  ggplot2::geom_violin(alpha = 0.5) +
  ggplot2::geom_boxplot(stat = "boxplot", alpha = 1, width = .5, outlier.size = 0.5) +
  ggpubr::stat_compare_means() +
  ggplot2::scale_fill_manual(values = cluster_cols) +
  ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
  ggplot2::labs(title = "CNV Scores", y = "CNV_scores") +
  Seurat::NoLegend())
dev.off()

readr::write_rds(TumorST, file.path(paths$intermediate, "04_TumorST_cnv_scored.rds.gz"), compress = "gz")

# Batch preparation for boundary-model training stops here: inferCNV labels and
# scores are retained, while no Mal/Bdy/nMal label is generated.  Keeping this
# switch inside the canonical script avoids maintaining a second inferCNV path.
if (env_flag("COTTRAZM_STOP_AFTER_CNV_SCORE", FALSE)) {
  log_step("COTTRAZM_STOP_AFTER_CNV_SCORE=1; CNV-scored object is ready and boundary definition is skipped")
  quit(save = "no", status = 0L, runLast = FALSE)
}

# 02C: define boundary.
TumorST <- readr::read_rds(file.path(paths$intermediate, "04_TumorST_cnv_scored.rds.gz"))
boundary_run_id <- params$boundary_run_id %||% ""
if (nzchar(boundary_run_id) && grepl("[/\\\\]", boundary_run_id)) {
  stop("boundary_run_id must not contain path separators: ", boundary_run_id, call. = FALSE)
}
out_dir_root <- file.path(paths$output, "05_boundary")
out_dir <- if (nzchar(boundary_run_id)) file.path(out_dir_root, boundary_run_id) else out_dir_root
intermediate_out_dir <- if (nzchar(boundary_run_id)) file.path(paths$intermediate, boundary_run_id) else paths$intermediate
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_out_dir, recursive = TRUE, showWarnings = FALSE)

boundary_mal_cluster_fraction <- params$boundary_mal_cluster_fraction
boundary_umap_mal_ratio <- params$boundary_umap_mal_ratio
boundary_expand_mal_radius <- params$boundary_expand_mal_radius
boundary_max_rounds <- params$boundary_max_rounds
if (!is.finite(boundary_mal_cluster_fraction) || boundary_mal_cluster_fraction <= 0 || boundary_mal_cluster_fraction > 1) {
  stop("boundary_mal_cluster_fraction must be in (0, 1]; got ", boundary_mal_cluster_fraction, call. = FALSE)
}
if (!is.finite(boundary_umap_mal_ratio) || boundary_umap_mal_ratio <= 0) {
  stop("boundary_umap_mal_ratio must be > 0; got ", boundary_umap_mal_ratio, call. = FALSE)
}
if (!is.finite(boundary_expand_mal_radius) || boundary_expand_mal_radius <= 0) {
  stop("boundary_expand_mal_radius must be > 0; got ", boundary_expand_mal_radius, call. = FALSE)
}
if (!is.finite(boundary_max_rounds) || boundary_max_rounds < 1) {
  stop("boundary_max_rounds must be >= 1; got ", boundary_max_rounds, call. = FALSE)
}

# UMAP 用来衡量表达/形态相似性；空间邻接表用来限制边界只能沿相邻 spot 扩展。
UMAPembeddings <- as.data.frame(TumorST@reductions$umap@cell.embeddings)
colnames(UMAPembeddings) <- c("x", "y")
position <- load_spaceranger_positions(paths$spaceranger, rownames(TumorST@meta.data))
position$spot.ids <- seq_len(nrow(position))
df_j <- vh_build_neighbors(position)
spatial_neighbor_qc <- vh_neighbor_qc(df_j)
vh_assert_neighbor_qc(spatial_neighbor_qc)
readr::write_tsv(spatial_neighbor_qc, file.path(out_dir, "neighbor_qc.tsv"))
readr::write_tsv(spatial_neighbor_qc, file.path(intermediate_out_dir, "neighbor_qc.tsv"))

plot_context <- sp_read_context(paths$spaceranger, load_image = TRUE)
scale_factor <- plot_context$tissue_lowres_scalef
plot_position <- position %>%
  tibble::rownames_to_column("cell_ID") %>%
  dplyr::mutate(
    X = imagecol * scale_factor,
    Y = imagerow * scale_factor
  )

make_spatial_boundary_plot <- function(label_df, label_col, label_cols, title = NULL) {
  point_df <- plot_position %>%
    dplyr::inner_join(label_df[, c("cell_ID", label_col), drop = FALSE], by = "cell_ID") %>%
    dplyr::transmute(cell_ID, X, Y, Location = factor(.data[[label_col]], levels = names(label_cols)))
  spot_polygons <- sp_make_spot_polygons(
    point_df,
    plot_context,
    segments = params$spatial_plot_spot_segments %||% 24L
  )
  sp_boundary_base(
    spot_polygons,
    plot_context,
    label_cols,
    title = title,
    show_image = plot_context$has_image,
    legend_name = label_col
  )
}

# The pure core is also called by the benchmark adapter. Reference/Filtered
# labels remain ineligible, and no analytic boundary is passed to this call.
configured_malignant_labels <- params$boundary_malignant_cnv_labels %||%
  params$malignant_cnv_labels
core_role <- if ("infercnv_role" %in% colnames(TumorST@meta.data)) {
  as.character(TumorST@meta.data$infercnv_role)
} else {
  rep("Observation", nrow(TumorST@meta.data))
}
boundary_attempts <- data.frame(
  malignant_label_n = 2L,
  malignant_cluster_fraction = boundary_mal_cluster_fraction,
  umap_malignant_ratio = boundary_umap_mal_ratio,
  stringsAsFactors = FALSE
)
if (is.null(configured_malignant_labels)) {
  fallback_attempts <- expand.grid(
    malignant_label_n = c(2L, 1L),
    malignant_cluster_fraction = c(0.30, 0.20, 0.10, 0.05),
    umap_malignant_ratio = c(0.50, 0.40, 0.30, 0.70, 1.00),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  fallback_attempts$priority <-
    ifelse(fallback_attempts$malignant_label_n == 2L, 0, 10) +
    abs(fallback_attempts$malignant_cluster_fraction - boundary_mal_cluster_fraction) * 10 +
    abs(fallback_attempts$umap_malignant_ratio - boundary_umap_mal_ratio)
  fallback_attempts <- fallback_attempts[
    order(fallback_attempts$priority, seq_len(nrow(fallback_attempts))),
    c("malignant_label_n", "malignant_cluster_fraction", "umap_malignant_ratio"),
    drop = FALSE
  ]
  boundary_attempts <- unique(rbind(boundary_attempts, fallback_attempts))
}
boundary_attempts$attempt_id <- sprintf("attempt_%02d", seq_len(nrow(boundary_attempts)))

boundary_core <- NULL
selected_boundary_attempt <- NULL
boundary_attempt_log <- data.frame(
  attempt_id = character(), malignant_label_n = integer(),
  malignant_cluster_fraction = numeric(), umap_malignant_ratio = numeric(),
  status = character(), detail = character(), stringsAsFactors = FALSE
)
for (attempt_index in seq_len(nrow(boundary_attempts))) {
  attempt <- boundary_attempts[attempt_index, , drop = FALSE]
  result <- tryCatch(
    ct_define_boundary(
      spot_ids = rownames(TumorST@meta.data),
      cnv_label = TumorST@meta.data$CNVLabel,
      cnv_score = TumorST@meta.data$cnv_score,
      infercnv_role = core_role,
      cluster = TumorST@meta.data$seurat_clusters,
      normal_score = TumorST@meta.data$NormalScore,
      embedding = as.matrix(UMAPembeddings[, c("x", "y"), drop = FALSE]),
      neighbors = df_j,
      malignant_labels = configured_malignant_labels,
      malignant_label_n = attempt$malignant_label_n,
      malignant_cluster_fraction = attempt$malignant_cluster_fraction,
      umap_malignant_ratio = attempt$umap_malignant_ratio,
      expand_malignant_radius = boundary_expand_mal_radius,
      maximum_rounds = boundary_max_rounds
    ),
    error = function(error) error
  )
  succeeded <- !inherits(result, "error")
  boundary_attempt_log <- rbind(
    boundary_attempt_log,
    data.frame(
      attempt_id = attempt$attempt_id,
      malignant_label_n = attempt$malignant_label_n,
      malignant_cluster_fraction = attempt$malignant_cluster_fraction,
      umap_malignant_ratio = attempt$umap_malignant_ratio,
      status = if (succeeded) "selected" else "failed",
      detail = if (succeeded) "" else conditionMessage(result),
      stringsAsFactors = FALSE
    )
  )
  if (succeeded) {
    boundary_core <- result
    selected_boundary_attempt <- attempt
    break
  }
}
readr::write_tsv(boundary_attempt_log, file.path(out_dir, "parameter_attempts.tsv"))
readr::write_tsv(boundary_attempt_log, file.path(intermediate_out_dir, "parameter_attempts.tsv"))
if (is.null(boundary_core)) {
  stop(
    "Cottrazm boundary definition failed for all ", nrow(boundary_attempts),
    " parameter attempts; see ", file.path(out_dir, "parameter_attempts.tsv"),
    call. = FALSE
  )
}
MalLabel <- boundary_core$selected_malignant_labels

boundary_params <- data.frame(
  parameter = c(
    "sample_name",
    "boundary_run_id",
    "neighbor_method",
    "boundary_parameter_attempt",
    "boundary_malignant_cnv_labels",
    "boundary_mal_cluster_fraction",
    "boundary_umap_mal_ratio",
    "boundary_expand_mal_radius",
    "boundary_max_rounds"
  ),
  value = c(
    sample_name,
    if (nzchar(boundary_run_id)) boundary_run_id else "legacy",
    "visium_array_hex_6_neighbor",
    selected_boundary_attempt$attempt_id,
    paste(MalLabel, collapse = ","),
    as.character(selected_boundary_attempt$malignant_cluster_fraction),
    as.character(selected_boundary_attempt$umap_malignant_ratio),
    as.character(boundary_expand_mal_radius),
    as.character(boundary_max_rounds)
  )
)
readr::write_tsv(boundary_params, file.path(out_dir, "params.tsv"))
readr::write_tsv(boundary_params, file.path(intermediate_out_dir, "params.tsv"))
TumorST@meta.data$Location <- boundary_core$location
TumorST@meta.data$Location <- factor(TumorST@meta.data$Location, levels = c("Mal", "Bdy", "nMal"))
TumorST@misc$boundary_params <- boundary_params
subset_ids <- unique(c(boundary_core$malignant_ids, boundary_core$boundary_ids,
                       boundary_core$normal_ids))
TumorSTn <- subset(TumorST, cells = subset_ids)
TumorSTn@meta.data$LabelNew <- TumorSTn@meta.data$Location

boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
final_df <- TumorST@meta.data %>%
  tibble::rownames_to_column("cell_ID") %>%
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))
boundary_plot <- make_spatial_boundary_plot(final_df, "Location", boundary_cols)
sp_save_plot(
  boundary_plot,
  file.path(out_dir, paste0(sample_name, "_BoundaryDefine.pdf")),
  width = params$spatial_plot_width %||% 8,
  height = params$spatial_plot_height %||% 7,
  fixed_layout = TRUE,
  map_width = params$spatial_plot_map_width %||% 6.45,
  legend_width = params$spatial_plot_legend_width %||% 1.55
)

location_counts <- as.data.frame(table(TumorST@meta.data$Location, useNA = "ifany"))
colnames(location_counts) <- c("Location", "n")
readr::write_tsv(location_counts, file.path(out_dir, "location_counts.tsv"))
readr::write_tsv(location_counts, file.path(intermediate_out_dir, "location_counts.tsv"))
readr::write_rds(TumorSTn, file.path(intermediate_out_dir, "05_TumorST_boundary_subset.rds.gz"), compress = "gz")
readr::write_rds(TumorST, file.path(intermediate_out_dir, "05_TumorST_boundary_defined.rds.gz"), compress = "gz")
