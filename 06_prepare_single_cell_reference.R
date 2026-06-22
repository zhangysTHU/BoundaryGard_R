# 06：准备 07 反卷积需要的单细胞参考。
# 推荐输入：
# - input/single_cell/sig_exp.rds.gz：gene x celltype signature 表达矩阵。
# - input/single_cell/clustermarkers_list.rds.gz：命名 list，元素名为 celltype，值为 marker genes。
# 备用输入：
# - input/single_cell/single_cell_seurat.rds + clustermarkers_list.rds.gz；
#   脚本会按 define_types.txt 指定的 metadata 列（默认 Majortypes）求平均表达。
# 输出：
# - intermediate/06_sig_exp.rds.gz：07/08 使用的 signature matrix。
# - intermediate/06_clustermarkers_list.rds.gz：07 marker enrichment 和 08 重构基因集合使用。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "readr"))

prepared_sig <- file.path(paths$single_cell, "sig_exp.rds.gz")
prepared_markers <- file.path(paths$single_cell, "clustermarkers_list.rds.gz")
out_sig <- file.path(paths$intermediate, "06_sig_exp.rds.gz")
out_markers <- file.path(paths$intermediate, "06_clustermarkers_list.rds.gz")

# 情况 1：已有整理好的 Cottrazm vignette 风格参考文件，直接复制到 intermediate。
if (file.exists(prepared_sig) && file.exists(prepared_markers)) {
  sig_exp <- readr::read_rds(prepared_sig)
  clustermarkers_list <- readr::read_rds(prepared_markers)
} else {
  # 情况 2：从 Seurat 单细胞对象中计算各细胞类型的平均表达 signature。
  sc_file <- file.path(paths$single_cell, "single_cell_seurat.rds")
  markers_file <- file.path(paths$single_cell, "clustermarkers_list.rds.gz")
  define_types_file <- file.path(paths$single_cell, "define_types.txt")
  if (!file.exists(sc_file) || !file.exists(markers_file)) {
    stop("Provide either input/single_cell/sig_exp.rds.gz and clustermarkers_list.rds.gz, or single_cell_seurat.rds plus clustermarkers_list.rds.gz.", call. = FALSE)
  }
  se.obj <- readr::read_rds(sc_file)
  clustermarkers_list <- readr::read_rds(markers_file)
  DefineTypes <- if (file.exists(define_types_file)) readLines(define_types_file, warn = FALSE)[1] else "Majortypes"
  sig_scran <- unique(unlist(clustermarkers_list))
  norm_exp <- 2^(se.obj@assays$RNA@data) - 1
  id <- se.obj@meta.data[, DefineTypes]
  ExprSubset <- norm_exp[sig_scran, ]
  sig_exp <- NULL
  for (cell_type in unique(id)) {
    sig_exp <- cbind(sig_exp, apply(ExprSubset, 1, function(y) mean(y[which(id == cell_type)])))
  }
  colnames(sig_exp) <- unique(id)
}

readr::write_rds(sig_exp, out_sig, compress = "gz")
readr::write_rds(clustermarkers_list, out_markers, compress = "gz")
