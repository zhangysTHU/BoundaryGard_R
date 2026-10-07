# 00：R 版流水线全局配置。
# 本脚本不产生生物学结果；它负责固定随机种子、定位 scripts_format_R 根目录、
# 创建 input/intermediate/output/resources 等目录，并集中管理所有主脚本共享的参数。
# input/<样本名>/ 存放输入；intermediate/<样本名>/ 通常存放下一步真正读取的
# RDS/TSV 中间结果；output/<样本名>/ 存放图表、表格和可直接用于下游分析的结构化结果矩阵。
options(stringsAsFactors = FALSE)
set.seed(666)

# 空值合并操作符：x 为 NULL 时使用默认值 y，04_lsgi_gradient.R 中也会用到。
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

# 默认路径仍相对 scripts_format_R；批量或 benchmark 运行必须通过独立 run root
# 或三个显式 root 覆盖，避免不同 scenario 共享可写目录。
resolve_configured_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  path <- if (identical(value, "")) default else path.expand(value)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

configured_run_root <- Sys.getenv("COTTRAZM_RUN_DIR", unset = "")
run_root <- if (nzchar(configured_run_root)) {
  normalizePath(path.expand(configured_run_root), winslash = "/", mustWork = FALSE)
} else {
  NA_character_
}
default_input_root <- if (is.na(run_root)) file.path(script_dir, "input") else file.path(run_root, "input")
default_intermediate_root <- if (is.na(run_root)) file.path(script_dir, "intermediate") else file.path(run_root, "intermediate")
default_output_root <- if (is.na(run_root)) file.path(script_dir, "output") else file.path(run_root, "output")
input_root <- resolve_configured_path("COTTRAZM_INPUT_ROOT", default_input_root)
intermediate_root <- resolve_configured_path("COTTRAZM_INTERMEDIATE_ROOT", default_intermediate_root)
output_root <- resolve_configured_path("COTTRAZM_OUTPUT_ROOT", default_output_root)
resources_root <- resolve_configured_path("COTTRAZM_RESOURCES_ROOT", file.path(script_dir, "resources"))
if (!grepl("^[A-Za-z0-9._-]+$", sample_name)) {
  stop("COTTRAZM_SAMPLE_NAME contains unsafe path characters: ", sQuote(sample_name), call. = FALSE)
}
sample_input_dir <- file.path(input_root, sample_name)
sample_intermediate_dir <- file.path(intermediate_root, sample_name)
sample_output_dir <- file.path(output_root, sample_name)

paths <- list(
  run = run_root,
  input = sample_input_dir,
  spaceranger = file.path(sample_input_dir, "spaceranger_outs"),
  single_cell = file.path(sample_input_dir, "single_cell"),
  intermediate = sample_intermediate_dir,
  output = sample_output_dir,
  resources = resources_root
)

configured_dirs <- unique(c(
  input_root, intermediate_root, output_root, resources_root,
  unlist(paths[!vapply(paths, function(x) anyNA(x), logical(1))], use.names = FALSE)
))
invisible(lapply(configured_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

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

env_bool <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  parsed <- tolower(trimws(value))
  if (parsed %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (parsed %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  stop(name, " must be boolean; got ", sQuote(value), call. = FALSE)
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
  plot_mode = {
    value <- tolower(Sys.getenv("COTTRAZM_PLOT_MODE", unset = "full"))
    if (!value %in% c("none", "summary", "full")) {
      stop("COTTRAZM_PLOT_MODE must be none, summary, or full; got ", sQuote(value), call. = FALSE)
    }
    value
  },
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
  # 默认设为 NULL，表示不在 00_config.R 中写死标签；02_boundary_definition.R 会根据
  # 每个 CNVLabel 的 cnv_score 中位数自动选出最高的 2 个有效 Observation 标签。
  # 注意：CNVLabel 的数字编号本身通常只是聚类编号，不一定代表恶性程度高低；
  # 因此自动选择依据是 cnv_score，而不是简单取数值最大的 label 名。
  # 如果已经人工复核过 CNV 空间分布、marker 表达和组织形态，也可以手动改为
  # c("6", "8") 这类明确标签，覆盖自动选择。
  malignant_cnv_labels = NULL,
  boundary_run_id = env_chr("COTTRAZM_BOUNDARY_RUN_ID", NULL),
  # 边界识别步骤专用的恶性 CNV 标签。
  # 默认同样设为 NULL，使 02_boundary_definition.R 自动从 cnv_score 最高的 2 个
  # CNVLabel 生成恶性种子。也可以通过环境变量 COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS
  # 覆盖，例如 "6,8"；若设为 "null"、"none" 或 "auto"，会返回 NULL 并保持自动选择。
  boundary_malignant_cnv_labels = env_chr_vec("COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS", NULL),
  boundary_mal_cluster_fraction = env_num("COTTRAZM_BOUNDARY_MAL_CLUSTER_FRACTION", 0.3),
  boundary_umap_mal_ratio = env_num("COTTRAZM_BOUNDARY_UMAP_MAL_RATIO", 0.5),
  boundary_expand_mal_radius = env_num("COTTRAZM_BOUNDARY_EXPAND_MAL_RADIUS", 1.0),
  boundary_max_rounds = env_int("COTTRAZM_BOUNDARY_MAX_ROUNDS", 6),
  # 反卷积矩阵或单细胞注释中代表恶性上皮细胞的细胞类型名称。
  # 该名称必须与 03_spatial_deconvolution.R 读取到的细胞类型列名/注释名完全一致；
  # 如果输入数据使用 "Malignant"、"Tumor epithelial" 等其他命名，需要在这里同步修改。
  decon_malignant_cluster = "Malignant epithelial cells",
  # 反卷积中用于表示非恶性/总体上皮组织成分的细胞类型名称。
  # 下游会用它和恶性上皮、基质成分一起构建空间组成结果；名称不匹配会导致对应成分缺失。
  decon_tissue_cluster = "Epithelial cells",
  # 反卷积中代表基质细胞的细胞类型名称，默认使用成纤维细胞。
  # 如果单细胞参考中基质细胞被命名为 "Fibroblasts"、"CAF" 或其他标签，需要改成实际名称。
  decon_stromal_cluster = "Fibroblast cells",
  # 可选 08 空间重构步骤选取的 Location 区域。
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
  # Downstream spatial-map export and geometry contract.  The 8 x 7 inch
  # canvas preserves a fixed full-image map area plus a right-hand legend
  # strip.  Biological spot size is always read from Space Ranger rather than
  # encoded as a physical geom_point size.
  spatial_plot_width = 8,
  spatial_plot_height = 7,
  spatial_plot_map_width = 6.45,
  spatial_plot_legend_width = 1.55,
  spatial_plot_spot_segments = 24,
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
  # 05_boundary_related_lsgi_arrows.R 默认使用的边界相关箭头筛选策略。
  # 当前默认 local_broad：LSGI local regression 邻域中 Bdy spot 数量 >= 5 且比例 >= 0.10。
  # 可选值：
  # - local_broad, local_relaxed, local_primary, local_strict
  # - partition_any, partition_relaxed, partition_primary, partition_strict
  # - primary_consensus, primary_union
  # - all 表示遍历全部内置策略；也可用逗号分隔多个策略，例如 "local_broad,primary_union"。
  boundary_lsgi_arrow_strategy = env_chr("COTTRAZM_BOUNDARY_LSGI_ARROW_STRATEGY", "local_broad"),
  # 06_tumor_arrow_guided_boundary_profile.R 默认参数。
  # 使用 05 中已经筛选出的 tumor/cancer epithelial arrows 作为局部跨边界 profile 的方向轴。
  tumor_profile_strategy = env_chr("COTTRAZM_TUMOR_PROFILE_STRATEGY", "local_broad"),
  tumor_profile_tumor_components = env_chr_vec(
    "COTTRAZM_TUMOR_PROFILE_TUMOR_COMPONENTS",
    c("Cancer.Epithelial", "Tumor")
  ),
  tumor_profile_methods = env_chr_vec(
    "COTTRAZM_TUMOR_PROFILE_METHODS",
    c("cell_component", "marker_module", "single_gene")
  ),
  tumor_profile_positive_control_features = env_chr_vec(
    "COTTRAZM_TUMOR_PROFILE_POSITIVE_CONTROL_FEATURES",
    c("Cancer.Epithelial")
  ),
  tumor_profile_features = env_chr_vec(
    "COTTRAZM_TUMOR_PROFILE_FEATURES",
    c(
      "REDOX_SCORE", "BM_FRS", "MYELOID_FOLR2_SUPPORT",
      "CAF_REDOX_SUPPORT", "ENDOTHELIAL_REDOX_SUPPORT",
      "SLC7A11", "GPX4", "ACSL4", "GCLC", "GCLM", "GSS", "TXNRD1", "NQO1"
    )
  ),
  tumor_profile_tube_half_width_mode = env_chr("COTTRAZM_TUMOR_PROFILE_TUBE_HALF_WIDTH_MODE", "half_pitch"),
  tumor_profile_zero_width_modes = env_chr_vec(
    "COTTRAZM_TUMOR_PROFILE_ZERO_WIDTH_MODES",
    c("spot_radius", "half_pitch")
  ),
  tumor_profile_bin_width_mode = env_chr("COTTRAZM_TUMOR_PROFILE_BIN_WIDTH_MODE", "half_pitch"),
  tumor_profile_bdy_cluster_gap_mode = env_chr("COTTRAZM_TUMOR_PROFILE_BDY_CLUSTER_GAP_MODE", "one_pitch"),
  tumor_profile_min_profile_spots = env_int("COTTRAZM_TUMOR_PROFILE_MIN_PROFILE_SPOTS", 3),
  tumor_profile_require_normal_spot = env_bool("COTTRAZM_TUMOR_PROFILE_REQUIRE_NORMAL_SPOT", TRUE),
  tumor_profile_normal_location = env_chr("COTTRAZM_TUMOR_PROFILE_NORMAL_LOCATION", "nMal"),
  # LSGI 箭头长度归一化方式：
  # "global" 表示所有细胞组分的箭头一起归一化，长度可跨组分比较；
  # "by_component" 表示每个细胞组分内部单独归一化，只比较同一组分内的梯度强弱。
  lsgi_arrow_length_normalization = "global",
  # LSGI embedding/component 来源。cell_component 使用 07_DeconData；
  # nmf 自动在 lsgi_nmf_ranks 中选择 k；其余三类从 resources/lsgi_embedding_catalog/*.tsv 读取。
  lsgi_component_methods = c("cell_component", "nmf", "marker_module", "pathway", "single_gene"),
  lsgi_embedding_catalog_dir = file.path(paths$resources, "lsgi_embedding_catalog"),
  lsgi_expression_assay = "Spatial",
  lsgi_expression_layer = "counts",
  lsgi_expression_scale_factor = 10000,
  # gene set/module 至少需要命中这么多个正向基因才进入 LSGI；single_gene 固定按 1 个基因处理。
  lsgi_min_feature_genes = 2,
  # NMF 自动 k 选择参数；04_lsgi_gradient.R 会在这些 k 上做 singlet cross-validation。
  lsgi_nmf_ranks = 6:10,
  lsgi_nmf_cv_replicates = 3,
  lsgi_nmf_cv_tol = 1e-3,
  lsgi_nmf_cv_maxit = 100,
  lsgi_nmf_final_tol = 1e-4,
  lsgi_nmf_final_maxit = 300,
  lsgi_nmf_rank_error_tolerance = 0.01,
  lsgi_nmf_test_density = 0.05,
  lsgi_nmf_l1 = 0.01,
  lsgi_nmf_l2 = 0,
  lsgi_nmf_threads = 1,
  lsgi_nmf_precision = "double",
  lsgi_nmf_assay = "Spatial",
  lsgi_nmf_layer = "counts",
  lsgi_nmf_top_genes = 2000,
  lsgi_nmf_min_gene_spots = 10,
  lsgi_nmf_scale_factor = 10000,
  lsgi_nmf_seed = 666,
  # 以下 ID 必须存在于 resources/lsgi_embedding_catalog/marker_modules.tsv。
  # 设为 character(0) 或 NULL 时，04_lsgi_gradient.R 会计算该 TSV 中所有条目。
  lsgi_marker_module_ids = c(
    "CXCL9_MACROPHAGE", "SPP1_MACROPHAGE", "CD8_CYTOTOXICITY",
    "TUMOR_IFN_M5", "TUMOR_HLA_M6", "CAF_ECM_BARRIER",
    "CONTRACTILE_BARRIER", "PEMT_INVASION", "LAMININ_332_INTERFACE",
    "HYPOXIA", "REDOX_SCORE", "BM_FRS", "MYELOID_FOLR2_SUPPORT",
    "CAF_REDOX_SUPPORT", "ENDOTHELIAL_REDOX_SUPPORT",
    "TLS_ORGANIZATION", "MATURE_ADIPOCYTE",
    "ADIPOCYTE_LIPOLYSIS_FIELD", "NORMAL_LIKE_ENDOTHELIUM",
    "TUMOR_ENDOTHELIUM", "PERIVASCULAR_LYMPHATIC", "RESIDUAL_RECURRENCE"
  ),
  # 以下 ID 必须存在于 resources/lsgi_embedding_catalog/pathways.tsv。
  lsgi_pathway_ids = c(
    "HALLMARK_IFN_GAMMA_RESPONSE", "REACTOME_MHC_I_ANTIGEN_PRESENTATION",
    "REACTOME_MHC_II_ANTIGEN_PRESENTATION", "REACTOME_CHEMOKINE_RECEPTORS",
    "HALLMARK_EMT", "HALLMARK_TGF_BETA_SIGNALING",
    "REACTOME_ECM_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION",
    "KEGG_FOCAL_ADHESION", "HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS",
    "HALLMARK_TNFA_NFKB", "HALLMARK_COMPLEMENT", "HALLMARK_UPR",
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_ADIPOGENESIS", "HALLMARK_CHOLESTEROL_HOMEOSTASIS"
  ),
  # 以下 ID 必须存在于 resources/lsgi_embedding_catalog/single_genes.tsv。
  lsgi_single_gene_ids = c(
    "CXCL9", "CXCL10", "IFNG", "CCL5", "GZMB", "IDO1", "HLA_DRA", "CD74",
    "SPP1", "GPNMB", "MMP9", "ASPN", "FAP", "MMP11", "COL1A1", "COL5A1",
    "POSTN", "ACTA2", "MMP2", "SPARC", "TGFB3", "S100A4", "SERPINE1",
    "LAMC2", "LAMA3", "ITGB4", "MMP14", "KRT17", "ANXA1", "LGALS3",
    "VIM", "SNAI1", "ZEB1", "CDH2", "NDRG1", "VEGFA", "EGLN3", "CA9",
    "SLC7A11", "GPX4", "ACSL4", "GCLC", "GCLM", "GSS", "TXNRD1", "NQO1",
    "ACKR1", "SELENOP", "APOD", "FGF7", "PLVAP", "KDR", "RGS5", "PDGFRB",
    "LYVE1", "CCL21", "PLIN1", "ADIPOQ", "PPARG", "FABP4", "HNF4A",
    "AQP7", "LIPE", "BNIP3", "CCL13", "IGF1", "FGF2", "S100A9",
    "S100A7", "SLPI", "CHI3L1", "SERPINA3", "AZGP1"
  )
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
