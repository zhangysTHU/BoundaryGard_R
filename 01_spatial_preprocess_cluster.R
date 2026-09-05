# 01: Spatial preprocessing, QC, morphology-adjusted clustering.

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "Matrix", "reticulate", "ggplot2", "ggpubr", "cowplot", "openxlsx", "readr"))

# 01A: preprocess ST.
out_dir <- file.path(paths$output, "01_preprocess")
dir.create(file.path(out_dir, "QC"), recursive = TRUE, showWarnings = FALSE)

# Space Ranger 可能提供 matrix 目录，也可能提供 filtered_feature_bc_matrix.h5；两者择一读取。
matrix_dir <- file.path(paths$spaceranger, "filtered_feature_bc_matrix")
h5_file <- file.path(paths$spaceranger, "filtered_feature_bc_matrix.h5")
spatial_dir <- file.path(paths$spaceranger, "spatial")

if (dir.exists(matrix_dir)) {
  xdata <- Seurat::Read10X(data.dir = matrix_dir)
} else if (file.exists(h5_file)) {
  xdata <- Seurat::Read10X_h5(filename = h5_file)
} else {
  stop("No Space Ranger expression matrix found under input/spaceranger_outs.", call. = FALSE)
}

# 构建 Seurat 对象，并把 H&E 图像对象挂到 image slot，供 SpatialDimPlot/SpatialFeaturePlot 使用。
TumorST <- Seurat::CreateSeuratObject(counts = xdata, project = sample_name, min.cells = 0, assay = "Spatial")
Ximage <- Seurat::Read10X_Image(image.dir = spatial_dir)
Seurat::DefaultAssay(Ximage) <- "Spatial"
Ximage <- Ximage[colnames(TumorST)]
TumorST[["image"]] <- Ximage

# 线粒体比例是空间/单细胞 QC 常用指标，后续不直接过滤，但用于人工检查样本质量。
TumorST[["Mito.percent"]] <- Seurat::PercentageFeatureSet(TumorST, pattern = "^MT-")

pdf(file.path(out_dir, "QC", "Vlnplot.pdf"), width = 6, height = 4)
p <- Seurat::VlnPlot(TumorST, features = c("nFeature_Spatial", "nCount_Spatial", "Mito.percent"), pt.size = 0, combine = FALSE)
p <- lapply(p, function(x) x + Seurat::NoLegend() + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 0)))
print(cowplot::plot_grid(plotlist = p, ncol = 3))
dev.off()

pdf(file.path(out_dir, "QC", "featureplot.pdf"), width = 7, height = 7)
p <- Seurat::SpatialFeaturePlot(TumorST, features = c("nFeature_Spatial", "nCount_Spatial", "Mito.percent"), combine = FALSE)
p <- lapply(p, function(x) x + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 0)))
print(cowplot::plot_grid(plotlist = p, ncol = 3))
dev.off()

QCData <- TumorST@meta.data[, c("nCount_Spatial", "nFeature_Spatial", "Mito.percent")]
openxlsx::write.xlsx(QCData, file.path(out_dir, "QC", "QCData.xlsx"), overwrite = TRUE)
readr::write_rds(TumorST, file.path(paths$intermediate, "01_TumorST_preprocessed.rds.gz"), compress = "gz")

# 01B: morphology-adjusted clustering.
out_dir <- file.path(paths$output, "02_morphology_cluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TumorST <- readr::read_rds(file.path(paths$intermediate, "01_TumorST_preprocessed.rds.gz"))

# 通过 reticulate 调 Python/stLearn 生成 SME 形态校正矩阵；优先绑定统一环境中的 python，
# 这样即使当前 shell 没有先 conda activate，R 版流程也能稳定找到 stlearn。
if (file.exists(params$python_bin)) {
  reticulate::use_python(params$python_bin, required = TRUE)
} else {
  reticulate::use_condaenv(params$python_conda_env, required = TRUE)
}
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
