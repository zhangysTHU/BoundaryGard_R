# 02: inferCNV, CNV scoring, and Mal/Bdy/nMal boundary definition.

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
options(scipen = 100)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "infercnv", "readr", "ape", "dendextend", "ggplot2", "ggpubr", "magrittr", "dplyr", "purrr", "tibble", "assertthat", "png", "grid"))

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

compute_interspot_distances <- function(position, scale.factor = 1.05) {
  # 用阵列 row/col 与图像像素 imagerow/imagecol 的线性关系估计相邻 spot 半径。
  cols <- c("row", "col", "imagerow", "imagecol")
  assertthat::assert_that(all(cols %in% colnames(position)))
  list(
    xdist = coef(lm(position$imagecol ~ position$col))[2],
    ydist = coef(lm(position$imagerow ~ position$row))[2]
  ) |>
    within(radius <- (abs(xdist) + abs(ydist)) * scale.factor)
}

find_neighbors <- function(position, radius, method = c("manhattan", "euclidean")) {
  # 根据像素坐标计算所有 spot 的邻居表；返回 named list，名字是 spot barcode。
  method <- match.arg(method)
  pdist <- as.matrix(stats::dist(as.matrix(position[, c("imagecol", "imagerow")]), method = method))
  neighbors <- (pdist <= radius & pdist > 0)
  df_j2 <- sapply(seq_len(nrow(position)), function(x) as.vector(names(neighbors[x, ])[which(neighbors[x, ])]))
  names(df_j2) <- rownames(position)
  df_j2
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
dists <- compute_interspot_distances(position = position, scale.factor = 1.05)
df_j <- find_neighbors(position = position, radius = dists$radius, method = "manhattan")

slice <- names(TumorST@images)[1]
scale_factor <- if (!is.na(slice) && nzchar(slice)) {
  TumorST@images[[slice]]@scale.factors$lowres %||% 1
} else {
  1
}
plot_position <- position %>%
  tibble::rownames_to_column("cell_ID") %>%
  dplyr::mutate(
    X = imagecol * scale_factor,
    Y = imagerow * scale_factor
  )
img_path <- file.path(paths$spaceranger, "spatial", "tissue_lowres_image.png")
has_image <- file.exists(img_path)
if (has_image) {
  img <- png::readPNG(img_path)
  img_grob <- grid::rasterGrob(
    img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc")
  )
}

make_spatial_boundary_plot <- function(label_df, label_col, label_cols, title = NULL) {
  point_df <- plot_position %>%
    dplyr::inner_join(label_df[, c("cell_ID", label_col), drop = FALSE], by = "cell_ID") %>%
    dplyr::mutate(label = factor(.data[[label_col]], levels = names(label_cols)))

  if (has_image) {
    point_df <- point_df %>%
      dplyr::mutate(Y = -Y)
    p <- ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grob = img_grob,
        xmin = 0,
        xmax = ncol(img),
        ymin = -nrow(img),
        ymax = 0
      ) +
      ggplot2::geom_point(
        data = point_df,
        ggplot2::aes(x = X, y = Y, fill = label),
        shape = 21,
        size = 1.8,
        stroke = 0.1,
        color = "grey20",
        alpha = 0.82,
        na.rm = TRUE
      ) +
      ggplot2::coord_fixed(
        ratio = 1,
        xlim = c(0, ncol(img)),
        ylim = c(-nrow(img), 0),
        expand = FALSE,
        clip = "on"
      )
  } else {
    p <- ggplot2::ggplot(point_df, ggplot2::aes(x = X, y = Y, fill = label)) +
      ggplot2::geom_point(shape = 21, size = 1.8, stroke = 0.1, color = "grey25", alpha = 0.9, na.rm = TRUE) +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_fixed()
  }

  p <- p +
    ggplot2::scale_fill_manual(values = label_cols, drop = FALSE, name = label_col) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}

# 若未手动指定恶性 CNV 标签，就从有效 Observation 标签中取 CNV score 中位数最高的两个。
# Reference/Filtered 绝不能成为恶性种子；全 NA 的低质量 spot 也不参与排序。
MalLabel <- params$boundary_malignant_cnv_labels %||% params$malignant_cnv_labels
forbidden_labels <- c("Normal", "Filtered")
if (is.null(MalLabel)) {
  role_ok <- if ("infercnv_role" %in% colnames(TumorST@meta.data)) {
    TumorST@meta.data$infercnv_role == "Observation"
  } else {
    rep(TRUE, nrow(TumorST@meta.data))
  }
  valid <- is.finite(TumorST@meta.data$cnv_score) &
    !TumorST@meta.data$CNVLabel %in% forbidden_labels & role_ok
  label_scores <- vapply(
    split(TumorST@meta.data$cnv_score[valid], TumorST@meta.data$CNVLabel[valid]),
    stats::median,
    FUN.VALUE = numeric(1),
    na.rm = TRUE
  )
  label_scores <- label_scores[is.finite(label_scores)]
  if (length(label_scores) < 2) stop("Fewer than two valid Observation CNV labels for malignant seeding", call. = FALSE)
  MalLabel <- names(sort(label_scores, decreasing = TRUE))[1:2]
} else {
  MalLabel <- setdiff(as.character(MalLabel), forbidden_labels)
  if (length(MalLabel) < 1) stop("boundary_malignant_cnv_labels contains only forbidden Reference/Filtered labels", call. = FALSE)
}

boundary_params <- data.frame(
  parameter = c(
    "sample_name",
    "boundary_run_id",
    "boundary_malignant_cnv_labels",
    "boundary_mal_cluster_fraction",
    "boundary_umap_mal_ratio",
    "boundary_expand_mal_radius",
    "boundary_max_rounds"
  ),
  value = c(
    sample_name,
    if (nzchar(boundary_run_id)) boundary_run_id else "legacy",
    paste(MalLabel, collapse = ","),
    as.character(boundary_mal_cluster_fraction),
    as.character(boundary_umap_mal_ratio),
    as.character(boundary_expand_mal_radius),
    as.character(boundary_max_rounds)
  )
)
readr::write_tsv(boundary_params, file.path(out_dir, "params.tsv"))
readr::write_tsv(boundary_params, file.path(intermediate_out_dir, "params.tsv"))

MalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$CNVLabel %in% MalLabel, ])
NormalCluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[, c("seurat_clusters", "NormalScore")], TumorST@meta.data$seurat_clusters), function(x) mean(x$NormalScore))), decreasing = TRUE)[1]]
NormalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == NormalCluster, ])

CNV_seurat_df <- as.data.frame.array(table(TumorST@meta.data$CNVLabel, TumorST@meta.data$seurat_clusters))[MalLabel, , drop = FALSE]
ClusterID <- c()
for (cluster in levels(TumorST@meta.data$seurat_clusters)) {
  if (sum(CNV_seurat_df[, cluster]) > table(TumorST@meta.data$seurat_clusters)[cluster] * boundary_mal_cluster_fraction) {
    ClusterID <- c(ClusterID, cluster)
  }
}

# 计算每个恶性 cluster 的 UMAP 中心，以及正常参考中心。
CiMal <- tibble::tibble(cluster = ClusterID) |>
  dplyr::mutate(sub_MalCellID = purrr::map(cluster, function(x) intersect(MalCellID, rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == x, ])))) |>
  dplyr::mutate(sub_CiMal = purrr::map(sub_MalCellID, function(x) apply(UMAPembeddings[x, , drop = FALSE], 2, mean)))
CiNormal <- apply(UMAPembeddings[NormalCellID, , drop = FALSE], 2, mean)

MalCellIDsi <- purrr::map2(CiMal$sub_MalCellID, CiMal$sub_CiMal, function(x, y) {
  lapply(x, function(id) {
    pos <- UMAPembeddings[id, ]
    rt <- sqrt(sum((pos - y)^2))
    rn <- sqrt(sum((pos - CiNormal)^2))
    if (rt < boundary_umap_mal_ratio * rn) id
  }) |> unlist()
}) |> unlist()

MalCellIDL <- lapply(MalCellIDsi, function(name) ifelse(length(df_j[[name]][df_j[[name]] %in% MalCellIDsi]) == 0, name, NA)) |> unlist() |> stats::na.omit()
BdyCellID <- NULL
nbrs_of_MalL <- nbrs(df_j = df_j, MalCellIDAdd = MalCellIDL, CellIDRaw = c(MalCellIDsi, NormalCellID, BdyCellID))
ClusterL <- do.call(rbind, lapply(names(nbrs_of_MalL), function(celll) {
  sub <- TumorST@meta.data[celll, ]$seurat_clusters
  do.call(rbind, lapply(nbrs_of_MalL[[celll]], function(idl) {
    pos <- UMAPembeddings[idl, ]
    rt <- sqrt(sum((pos - unlist(CiMal[CiMal$cluster == sub, ]$sub_CiMal))^2))
    rn <- sqrt(sum((pos - CiNormal)^2))
    data.frame(CellID = idl, Location = ifelse(rt < boundary_umap_mal_ratio * rn, "Mal", "Bdy"))
  }))
}))

ClusterL_tab <- as.data.frame.array(table(ClusterL$CellID, ClusterL$Location))
ClusterL_tab$Location <- if ("Mal" %in% colnames(ClusterL_tab)) ifelse(ClusterL_tab$Mal > 0, "Mal", "Bdy") else rep("Bdy", nrow(ClusterL_tab))
ClusterL <- data.frame(CellID = rownames(ClusterL_tab), Location = as.character(ClusterL_tab$Location))

MalCellID <- c(MalCellIDsi, ClusterL[ClusterL$Location == "Mal", ]$CellID)
BdyCellID <- as.character(ClusterL[ClusterL$Location == "Bdy", ]$CellID)
MalCellIDN <- MalCellID
n <- 1
TumorSTn <- TumorST
Clustern <- rbind(
  data.frame(CellID = NormalCellID, Location = rep("Normal", length(NormalCellID))),
  data.frame(CellID = MalCellID, Location = rep("Mal", length(MalCellID))),
  data.frame(CellID = BdyCellID, Location = rep("Bdy", length(BdyCellID)))
)

# 逐圈扩展：从 Mal 边缘找未标注邻居，根据 UMAP 距离判断是继续恶性扩展还是边界。
repeat {
  if (length(MalCellIDN) < 3) break
  nbrs_of_Mal <- nbrs(df_j = df_j, MalCellIDAdd = MalCellIDN, CellIDRaw = c(MalCellID, NormalCellID, BdyCellID))
  if (length(unique(unlist(nbrs_of_Mal))) < 3) break
  TumorSTn <- subset(TumorST, cells = c(unique(unlist(nbrs_of_Mal)), MalCellID, NormalCellID, BdyCellID))
  TumorSTn@meta.data$Label <- Clustern$Location[match(rownames(TumorSTn@meta.data), Clustern$CellID)]
  TumorSTn@meta.data$Label <- factor(TumorSTn@meta.data$Label, levels = if (n == 1) c("Normal", "Bdy", "Mal") else c("Normal", "Bdy", "Mal", paste0("Mal", 1:(n - 1))))
  ClusterAdd <- ClusterUpdate(
    x = n,
    position = position,
    df_j = df_j,
    UMAPembeddings = UMAPembeddings,
    MalCellIDN = MalCellIDN,
    BdyCellID = BdyCellID,
    NormalCellID = NormalCellID,
    MalCellID = MalCellID,
    expand_mal_radius = boundary_expand_mal_radius
  )
  Clustern <- rbind(Clustern, ClusterAdd)
  TumorSTn@meta.data$LabelNew <- Clustern$Location[match(rownames(TumorSTn@meta.data), as.character(Clustern$CellID))]
  TumorSTn@meta.data$LabelNew <- factor(TumorSTn@meta.data$LabelNew, levels = c("Normal", "Bdy", "Mal", paste0("Mal", 1:n)))

  cols_n <- c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:(n + 1)])
  names(cols_n) <- levels(TumorSTn@meta.data$LabelNew)
  pdf(file.path(out_dir, paste0(sample_name, "_out_", n, ".pdf")), width = 7, height = 7)
  round_df <- TumorSTn@meta.data %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::mutate(LabelNew = factor(LabelNew, levels = names(cols_n)))
  print(make_spatial_boundary_plot(round_df, "LabelNew", cols_n))
  dev.off()

  MalCellID <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% c("Mal", paste0("Mal", 1:n)), ])
  MalCellIDN <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% paste0("Mal", n), ])
  BdyCellID <- rownames(TumorSTn@meta.data[TumorSTn$LabelNew == "Bdy", ])
  n <- n + 1
  if (n > boundary_max_rounds) break
}

# 最终折叠为三类：Mal 为恶性核心/扩展层，Bdy 为边界及邻近正常侧，剩余为 nMal。
position_all <- position
dists_all <- compute_interspot_distances(position = position_all, scale.factor = 1.05)
df_j_all <- find_neighbors(position = position_all, radius = dists_all$radius, method = "manhattan")
Mal_barcode <- rownames(TumorSTn@meta.data)[grep("Mal", TumorSTn@meta.data$LabelNew)]
Bdy_barcode <- rownames(TumorSTn@meta.data)[grep("Bdy", TumorSTn@meta.data$LabelNew)]
Normal_Bdy_barcode <- unique(unlist(nbrs(df_j = df_j_all, MalCellIDAdd = Mal_barcode, CellIDRaw = c(Bdy_barcode, Mal_barcode))))
nMal_barcode <- rownames(TumorST@meta.data)[!rownames(TumorST@meta.data) %in% c(Mal_barcode, Bdy_barcode, Normal_Bdy_barcode)]
Barcode_Ann <- data.frame(barcode = c(Mal_barcode, Bdy_barcode, Normal_Bdy_barcode, nMal_barcode), Location = c(rep("Mal", length(Mal_barcode)), rep("Bdy", length(c(Bdy_barcode, Normal_Bdy_barcode))), rep("nMal", length(nMal_barcode))))
TumorST@meta.data$Location <- Barcode_Ann$Location[match(rownames(TumorST@meta.data), Barcode_Ann$barcode)]
TumorST@meta.data$Location <- factor(TumorST@meta.data$Location, levels = c("Mal", "Bdy", "nMal"))
TumorST@misc$boundary_params <- boundary_params

pdf(file.path(out_dir, paste0(sample_name, "_BoundaryDefine.pdf")), width = 7, height = 7)
boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
final_df <- TumorST@meta.data %>%
  tibble::rownames_to_column("cell_ID") %>%
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))
print(make_spatial_boundary_plot(final_df, "Location", boundary_cols))
dev.off()

location_counts <- as.data.frame(table(TumorST@meta.data$Location, useNA = "ifany"))
colnames(location_counts) <- c("Location", "n")
readr::write_tsv(location_counts, file.path(out_dir, "location_counts.tsv"))
readr::write_tsv(location_counts, file.path(intermediate_out_dir, "location_counts.tsv"))
readr::write_rds(TumorSTn, file.path(intermediate_out_dir, "05_TumorST_boundary_subset.rds.gz"), compress = "gz")
readr::write_rds(TumorST, file.path(intermediate_out_dir, "05_TumorST_boundary_defined.rds.gz"), compress = "gz")
