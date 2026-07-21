# 00：R 版流水线全局配置。
# 本脚本不产生生物学结果；它负责固定随机种子、定位 scripts_format_R 根目录、
# 创建 input/intermediate/output/resources/lib 等目录，并集中管理所有主脚本共享的参数。
# input/<样本名>/ 存放输入；intermediate/<样本名>/ 通常存放下一步真正读取的
# RDS/TSV 中间结果；output/<样本名>/ 存放图表、表格和可直接用于下游分析的结构化结果矩阵。
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
if (!basename(script_dir) == "scripts_format_R") {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

sample_name <- Sys.getenv("COTTRAZM_SAMPLE_NAME", unset = "CRC1")
conda_root <- Sys.getenv("COTTRAZM_CONDA_ROOT", unset = "/lulabdata3/huangkeyun/zhangys/tools/miniforge3")
python_conda_env <- Sys.getenv("COTTRAZM_CONDA_ENV", unset = "BoundaryGrad")
python_bin <- Sys.getenv(
  "COTTRAZM_PYTHON",
  unset = file.path(conda_root, "envs", python_conda_env, "bin", "python")
)

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

env_chr <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) default else value
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) {
    stop(name, " must be numeric; got ", sQuote(value), call. = FALSE)
  }
  parsed
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

env_chr_vec <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  if (tolower(value) %in% c("null", "none", "auto")) {
    return(NULL)
  }
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parsed[nzchar(parsed)]
}

# 主流程参数。修改聚类分辨率、inferCNV 线程数、反卷积细胞类型名、重构区域等，优先改这里。
params <- list(
  cluster_resolution = 1.5,
  infercnv_assay = "Spatial",
  infercnv_threads = 30,
  infercnv_partition_method = "random_trees",
  infercnv_analysis_mode = "subclusters",
  infercnv_min_counts = 100,
  infercnv_min_features = 100,
  infercnv_reference_fraction = 0.06,
  infercnv_reference_min_spots = 150,
  infercnv_reference_max_spots = 400,
  infercnv_reference_immune_markers = c(
    "PTPRC", "CD2", "CD3D", "CD3E", "CD3G",
    "CD5", "CD7", "CD79A", "MS4A1", "CD19"
  ),
  infercnv_reference_epithelial_markers = c(
    "EPCAM", "KRT5", "KRT7", "KRT8", "KRT14",
    "KRT15", "KRT17", "KRT18", "KRT19"
  ),
  cnv_k = 8,
  python_conda_env = python_conda_env,
  python_bin = python_bin,
  # inferCNV/CNV 聚类结果中被视为恶性肿瘤区域的标签。
  # 这里的标签必须与上游 CNVLabel/CNV cluster 结果中的字符标签完全一致。
  # 默认设为 NULL，表示不在 00_config.R 中写死标签；05_define_boundary.R 会根据
  # 每个 CNVLabel 的 cnv_score 中位数自动选出最高的 2 个有效 Observation 标签。
  # 注意：CNVLabel 的数字编号本身通常只是聚类编号，不一定代表恶性程度高低；
  # 因此自动选择依据是 cnv_score，而不是简单取数值最大的 label 名。
  # 如果已经人工复核过 CNV 空间分布、marker 表达和组织形态，也可以手动改为
  # c("6", "8") 这类明确标签，覆盖自动选择。
  malignant_cnv_labels = NULL,
  boundary_run_id = env_chr("COTTRAZM_BOUNDARY_RUN_ID", NULL),
  # 边界识别步骤专用的恶性 CNV 标签。
  # 默认同样设为 NULL，使 05_define_boundary.R 自动从 cnv_score 最高的 2 个
  # CNVLabel 生成恶性种子。也可以通过环境变量 COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS
  # 覆盖，例如 "6,8"；若设为 "null"、"none" 或 "auto"，会返回 NULL 并保持自动选择。
  boundary_malignant_cnv_labels = env_chr_vec("COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS", NULL),
  boundary_mal_cluster_fraction = env_num("COTTRAZM_BOUNDARY_MAL_CLUSTER_FRACTION", 0.3),
  boundary_umap_mal_ratio = env_num("COTTRAZM_BOUNDARY_UMAP_MAL_RATIO", 0.5),
  boundary_expand_mal_radius = env_num("COTTRAZM_BOUNDARY_EXPAND_MAL_RADIUS", 1.0),
  boundary_max_rounds = env_int("COTTRAZM_BOUNDARY_MAX_ROUNDS", 6),
  # 反卷积矩阵或单细胞注释中代表恶性上皮细胞的细胞类型名称。
  # 该名称必须与 07_spatial_deconvolution.R 读取到的细胞类型列名/注释名完全一致；
  # 如果输入数据使用 "Malignant"、"Tumor epithelial" 等其他命名，需要在这里同步修改。
  decon_malignant_cluster = "Malignant epithelial cells",
  # 反卷积中用于表示非恶性/总体上皮组织成分的细胞类型名称。
  # 下游会用它和恶性上皮、基质成分一起构建空间组成结果；名称不匹配会导致对应成分缺失。
  decon_tissue_cluster = "Epithelial cells",
  # 反卷积中代表基质细胞的细胞类型名称，默认使用成纤维细胞。
  # 如果单细胞参考中基质细胞被命名为 "Fibroblasts"、"CAF" 或其他标签，需要改成实际名称。
  decon_stromal_cluster = "Fibroblast cells",
  # 空间重构步骤选取的 Location 区域。
  # 默认只重构边界区域 Bdy；该值必须存在于 TumorST@meta.data$Location 中。
  # 若需要同时分析多个区域，可改为字符向量，例如 c("Bdy", "Tumor")。
  recon_locations = "Bdy",
  diff_assay = "Spatial",
  diff_logfc_cutoff = 0.25,
  diff_fdr_cutoff = 0.05,
  volcano_p_cutoff_log10 = 2,
  volcano_label_n = 10,
  pie_scale = 0.4,
  scatterpie_alpha = 0.8,
  pie_border_color = "grey",
  # LSGI 箭头外观参数。
  lsgi_arrow_length_scale = 1.4,
  lsgi_arrow_linewidth = 1.0,
  lsgi_arrow_head_cm = 0.20,
  lsgi_arrow_head_angle = 30,
  # 用局部线性回归 R2 映射箭头头部张角，表达箭头方向估计的可信度。
  # R2 <= 0.3 映射到 22.5 度；R2 = 0.5 保持默认 30 度；R2 >= 0.7 映射到 45 度。
  lsgi_arrow_head_angle_by_r2 = TRUE,
  lsgi_arrow_head_angle_min = 22.5,
  lsgi_arrow_head_angle_mid = 30,
  lsgi_arrow_head_angle_max = 45,
  lsgi_arrow_head_angle_r2_min = 0.3,
  lsgi_arrow_head_angle_r2_mid = 0.5,
  lsgi_arrow_head_angle_r2_max = 0.7,
  lsgi_arrow_head_angle_step = 1,
  lsgi_arrow_closed = TRUE,
  # LSGI 箭头长度归一化方式：
  # "global" 表示所有细胞组分的箭头一起归一化，长度可跨组分比较；
  # "by_component" 表示每个细胞组分内部单独归一化，只比较同一组分内的梯度强弱。
  lsgi_arrow_length_normalization = "global"
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
