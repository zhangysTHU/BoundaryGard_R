# 02：H&E 形态校正聚类。
# 输入：intermediate/01_TumorST_preprocessed.rds.gz。
# 输出：
# - intermediate/02_TumorST_clustered.rds.gz：Seurat 对象，新增 Morph assay、PCA/UMAP、seurat_clusters、NormalScore。
# - intermediate/InferCNV/CellAnnotation.txt：两列无表头，CellID 和 seurat_clusters，供 03 inferCNV 分组。
# - output/02_morphology_cluster/CRC1_tile/：stLearn 切出的 H&E tiles。
# - output/02_morphology_cluster/*Cluster.pdf 和 *NormalScore.pdf：聚类与正常细胞 marker 分数图。
# 下游：03 用 cluster/NormalScore 选择 inferCNV reference；05 用 UMAP 和 cluster 做边界。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "Matrix", "reticulate", "ggplot2", "ggpubr", "cowplot", "readr"))

out_dir <- file.path(paths$output, "02_morphology_cluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TumorST <- readr::read_rds(file.path(paths$intermediate, "01_TumorST_preprocessed.rds.gz"))

# 通过 reticulate 调 Python/stLearn 生成 SME 形态校正矩阵；结果缓存为 MatrixMarket，避免重复切图。
reticulate::use_condaenv(params$python_conda_env, required = TRUE)
reticulate::source_python(file.path(paths$resources, "python", "Rusedtile.py"))
sme_mtx_file <- file.path(out_dir, paste0(sample_name, "_raw_SME_normalizeA.mtx"))
Adjusted_expr_mtx <- if (file.exists(sme_mtx_file)) {
  Matrix::readMM(sme_mtx_file)
} else {
  ME_normalize(inDir = paste0(paths$spaceranger, .Platform$file.sep), outDir = paste0(out_dir, .Platform$file.sep), sample = sample_name)
}
# stLearn 输出是 spot x gene，这里转成 Seurat 需要的 gene x spot 并放入 Morph assay。
spatial_counts <- Seurat::GetAssayData(TumorST, assay = "Spatial", layer = "counts")

ok <- try({
  rownames(Adjusted_expr_mtx) <- colnames(spatial_counts)
  colnames(Adjusted_expr_mtx) <- rownames(spatial_counts)
}, silent = TRUE)
if (inherits(ok, "try-error")) {
  Adjusted_expr_mtx <- Matrix::readMM(sme_mtx_file)
  rownames(Adjusted_expr_mtx) <- colnames(spatial_counts)
  colnames(Adjusted_expr_mtx) <- rownames(spatial_counts)
}

Adjusted_expr_mtxF <- t(as.matrix(Adjusted_expr_mtx))
MorphMatirxSeurat <- Seurat::CreateSeuratObject(counts = as(Adjusted_expr_mtxF, "sparseMatrix"))
MorphMatirxSeurat <- subset(MorphMatirxSeurat, cells = rownames(TumorST@meta.data))
TumorST@assays$Morph <- MorphMatirxSeurat@assays$RNA

# 标准 Seurat 聚类流程：NormalizeData -> HVG -> ScaleData -> PCA -> neighbors -> UMAP -> clusters。
TumorST <- Seurat::NormalizeData(TumorST, assay = "Morph")
TumorST <- Seurat::FindVariableFeatures(TumorST, mean.function = Seurat::ExpMean, dispersion.function = Seurat::LogVMR, assay = "Morph")
TumorST <- Seurat::ScaleData(TumorST, assay = "Morph")
TumorST <- Seurat::RunPCA(TumorST, npcs = 50, verbose = FALSE, assay = "Morph")
TumorST <- Seurat::FindNeighbors(TumorST, reduction = "pca", dims = 1:50, assay = "Morph")
TumorST <- Seurat::RunUMAP(TumorST, dims = 1:50, assay = "Morph")
TumorST <- Seurat::FindClusters(TumorST, resolution = params$cluster_resolution, algorithm = 1, graph.name = "Morph_snn")
TumorST@meta.data$seurat_clusters <- TumorST@meta.data[, paste0("Morph_snn_res.", params$cluster_resolution)]

pdf(file.path(out_dir, paste0(sample_name, "_Spatial_SeuratCluster.pdf")), width = 7, height = 7)
print(Seurat::SpatialDimPlot(TumorST, group.by = "seurat_clusters", cols = cluster_cols, pt.size.factor = 1, alpha = 0.8) +
  ggplot2::scale_fill_manual(values = cluster_cols) +
  ggplot2::labs(title = paste0("Resolution = ", params$cluster_resolution)))
dev.off()

pdf(file.path(out_dir, paste0(sample_name, "_UMAP_SeuratCluster.pdf")), width = 7, height = 7)
print(Seurat::DimPlot(TumorST, group.by = "seurat_clusters", cols = cluster_cols) +
  ggplot2::labs(title = paste0("Resolution = ", params$cluster_resolution)) +
  ggplot2::scale_fill_manual(values = cluster_cols))
dev.off()

# 用免疫/B 细胞 marker 近似正常细胞富集程度，NormalScore 最高的 cluster 后面作为 CNV reference。
normal_features <- c("PTPRC", "CD2", "CD3D", "CD3E", "CD3G", "CD5", "CD7", "CD79A", "MS4A1", "CD19")
morph_data <- Seurat::GetAssayData(TumorST, assay = "Morph", layer = "data")
present_normal_features <- rownames(morph_data) %in% normal_features
TumorST@meta.data$NormalScore <- if (any(present_normal_features)) {
  apply(morph_data[present_normal_features, , drop = FALSE], 2, mean)
} else {
  0
}

pdf(file.path(out_dir, paste0(sample_name, "_NormalScore.pdf")), width = 6, height = 4)
normal_score_df <- TumorST@meta.data[, c("seurat_clusters", "NormalScore")]
print(Seurat::VlnPlot(TumorST, features = "NormalScore", pt.size = 0, group.by = "seurat_clusters", cols = cluster_cols) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_hline(yintercept = max(unlist(lapply(split(normal_score_df, normal_score_df$seurat_clusters), function(x) median(x$NormalScore)))), linetype = "dashed") +
  ggpubr::stat_compare_means() +
  Seurat::NoLegend())
dev.off()

cellAnnotation <- data.frame(CellID = rownames(TumorST@meta.data), DefineTypes = TumorST@meta.data[, "seurat_clusters"])
dir.create(file.path(paths$intermediate, "InferCNV"), recursive = TRUE, showWarnings = FALSE)
write.table(cellAnnotation, file.path(paths$intermediate, "InferCNV", "CellAnnotation.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
readr::write_rds(TumorST, file.path(paths$intermediate, "02_TumorST_clustered.rds.gz"), compress = "gz")
