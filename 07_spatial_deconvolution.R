# 07：空间 spot 反卷积，估计每个 spot 中各参考细胞类型比例。
# 输入：
# - intermediate/05_TumorST_boundary_defined.rds.gz：需要 metadata$Location 和 seurat_clusters。
# - intermediate/06_sig_exp.rds.gz：gene x celltype signature。
# - intermediate/06_clustermarkers_list.rds.gz：celltype marker list。
# 输出：
# - intermediate/07_DeconData.rds.gz：data.frame，第一列 cell_ID，其余列为细胞类型比例；08/10/11 使用。
# - intermediate/07_decon_inputs.rds.gz：保存反卷积中间矩阵，便于调试。
# - intermediate/07_TumorST_for_decon.rds.gz：NormalizeData 后、带 Decon_topics 的对象。
# - output/07_spatial_deconvolution/DeconData.xlsx：反卷积结果的 Excel 版本。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "Matrix", "Rfast", "quadprog", "data.table", "magrittr", "dplyr", "tibble", "purrr", "openxlsx", "readr"))
source(file.path(paths$lib, "decon_helpers.R"))

out_dir <- file.path(paths$output, "07_spatial_deconvolution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TumorST <- readr::read_rds(file.path(paths$intermediate, "05_TumorST_boundary_defined.rds.gz"))
sig_exp <- readr::read_rds(file.path(paths$intermediate, "06_sig_exp.rds.gz"))
clustermarkers_list <- readr::read_rds(file.path(paths$intermediate, "06_clustermarkers_list.rds.gz"))

# 反卷积在 log-normalized 表达上做 marker enrichment，同时保留 nCount_Spatial 辅助低深度判断。
TumorST <- Seurat::NormalizeData(TumorST, assay = "Spatial")
TumorST@meta.data$Decon_topics <- paste(TumorST@meta.data$Location, TumorST@meta.data$seurat_clusters, sep = "_")
spatial_data <- Seurat::GetAssayData(TumorST, assay = "Spatial", layer = "data")
expr_values <- as.matrix(spatial_data)
nolog_expr <- 2^(expr_values) - 1
meta_data <- TumorST@meta.data[, c("nCount_Spatial", "Decon_topics", "Location")]

# meta_data 中每个细胞类型的列是该类型 top marker 在各 spot 中的平均表达，用于 topic 级候选筛选。
for (cluster in names(clustermarkers_list)) {
  cluster_markers <- clustermarkers_list[[cluster]][1:min(25, length(clustermarkers_list[[cluster]]))]
  cluster_score <- apply(spatial_data[rownames(spatial_data) %in% cluster_markers, , drop = FALSE], 2, mean)
  meta_data <- cbind(meta_data, cluster_score)
}
colnames(meta_data) <- c("nCount_Spatial", "Decon_topics", "Location", names(clustermarkers_list))

# 只使用空间表达和单细胞 signature 都存在的基因。
intersect_gene <- intersect(rownames(sig_exp), rownames(nolog_expr))
filter_sig <- sig_exp[intersect_gene, , drop = FALSE]
filter_expr <- nolog_expr[intersect_gene, , drop = FALSE]
filter_log_expr <- expr_values[intersect_gene, , drop = FALSE]
enrich_matrix <- get_enrich_matrix(filter_sig = filter_sig, clustermarkers_list = clustermarkers_list)
enrich_result <- enrich_analysis(filter_log_expr = filter_log_expr, enrich_matrix = enrich_matrix)

# 第一阶段：用 marker enrichment 和 DWLS 得到初始比例。
dwls_results <- spot_proportion_initial(
  enrich_matrix = enrich_matrix,
  enrich_result = enrich_result,
  filter_expr = filter_expr,
  filter_sig = filter_sig,
  clustermarkers_list = clustermarkers_list,
  meta_data = meta_data,
  malignant_cluster = params$decon_malignant_cluster,
  tissue_cluster = params$decon_tissue_cluster,
  stromal_cluster = params$decon_stromal_cluster
)

# 第二阶段：把初始比例二值化为候选细胞类型，再逐 spot 精修反卷积比例。
binary_matrix <- ifelse(dwls_results >= 0.01, 1, 0)
spot_proportion <- spot_deconvolution(
  expr = filter_expr,
  meta_data = meta_data,
  ct_exp = filter_sig,
  enrich_matrix = enrich_matrix,
  binary_matrix = binary_matrix
)

DeconData <- as.data.frame(t(spot_proportion))
DeconData <- tibble::rownames_to_column(DeconData, var = "cell_ID")

readr::write_rds(TumorST, file.path(paths$intermediate, "07_TumorST_for_decon.rds.gz"), compress = "gz")
readr::write_rds(list(filter_sig = filter_sig, filter_expr = filter_expr, filter_log_expr = filter_log_expr, enrich_matrix = enrich_matrix, enrich_result = enrich_result, meta_data = meta_data), file.path(paths$intermediate, "07_decon_inputs.rds.gz"), compress = "gz")
readr::write_rds(DeconData, file.path(paths$intermediate, "07_DeconData.rds.gz"), compress = "gz")
openxlsx::write.xlsx(DeconData, file.path(out_dir, "DeconData.xlsx"), overwrite = TRUE)
