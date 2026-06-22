# 00：R 版流水线全局配置。
# 本脚本不产生生物学结果；它负责固定随机种子、定位 scripts_format_R 根目录、
# 创建 input/intermediate/output/resources/lib 等目录，并集中管理所有主脚本共享的参数。
# input/<样本名>/ 存放输入；intermediate/<样本名>/ 通常存放下一步真正读取的
# RDS/TSV 中间结果；output/<样本名>/ 主要存放人工检查图表和表格。
options(stringsAsFactors = FALSE)
set.seed(666)

# 空值合并操作符：x 为 NULL 时使用默认值 y，11_lsgi_gradient.R 中也会用到。
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# 兼容 Rscript script.R、source("script.R") 和从脚本目录直接运行三种入口方式。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file_arg) > 0) sub("^--file=", "", script_file_arg[[1]]) else NULL
script_dir <- if (!is.null(script_file)) {
  normalizePath(dirname(script_file), winslash = "/", mustWork = FALSE)
} else if (!is.null(sys.frame(1)$ofile)) {
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

sample_name <- Sys.getenv("COTTRAZM_SAMPLE_NAME", unset = "CRC1")

# 所有路径均相对 scripts_format_R。
input_root <- file.path(script_dir, "input")
intermediate_root <- file.path(script_dir, "intermediate")
output_root <- file.path(script_dir, "output")
sample_input_dir <- file.path(input_root, sample_name)
sample_intermediate_dir <- file.path(intermediate_root, sample_name)
sample_output_dir <- file.path(output_root, sample_name)

paths <- list(
  input = sample_input_dir,
  spaceranger = file.path(sample_input_dir, "spaceranger_outs"),
  single_cell = file.path(sample_input_dir, "single_cell"),
  intermediate = sample_intermediate_dir,
  output = sample_output_dir,
  resources = file.path(script_dir, "resources"),
  lib = file.path(script_dir, "lib")
)

invisible(lapply(c(input_root, intermediate_root, output_root, unlist(paths)), dir.create, recursive = TRUE, showWarnings = FALSE))

# 主流程参数。修改聚类分辨率、inferCNV 线程数、反卷积细胞类型名、重构区域等，优先改这里。
params <- list(
  cluster_resolution = 1.5,
  infercnv_assay = "Spatial",
  infercnv_threads = 30,
  cnv_k = 8,
  python_conda_env = "cottrazm-py",
  malignant_cnv_labels = NULL,
  decon_malignant_cluster = "Malignant epithelial cells",
  decon_tissue_cluster = "Epithelial cells",
  decon_stromal_cluster = "Fibroblast cells",
  recon_locations = "Bdy",
  diff_assay = "Spatial",
  diff_logfc_cutoff = 0.25,
  diff_fdr_cutoff = 0.05,
  volcano_p_cutoff_log10 = 2,
  volcano_label_n = 10,
  pie_scale = 0.4,
  scatterpie_alpha = 0.8,
  pie_border_color = "grey"
)

# 全流程可能用到的包集合；每个脚本会按需传入子集给 load_required_packages()。
required_packages <- c(
  "Seurat", "magrittr", "dplyr", "ggplot2", "tibble", "purrr", "readr",
  "openxlsx", "cowplot", "ggpubr", "assertthat", "reticulate", "Matrix",
  "infercnv", "ape", "dendextend", "quadprog", "Rfast", "scatterpie",
  "jpeg", "png", "grid", "stringr", "ggrepel", "clusterProfiler",
  "org.Hs.eg.db"
)

load_required_packages <- function(pkgs = required_packages) {
  # 先检查命名空间是否可用，再 suppressPackageStartupMessages(library())，让报错更早、更清楚。
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(pkgs, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
}

# 离散类别图使用的调色板，例如 Seurat cluster、CNVLabel。
cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)
