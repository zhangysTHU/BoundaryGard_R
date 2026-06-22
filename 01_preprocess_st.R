# 01：读取 10x Visium / Space Ranger 输出并计算基础 QC。
# 输入：input/spaceranger_outs/，应包含 filtered_feature_bc_matrix(.h5 或目录) 和 spatial/ 图像坐标。
# 输出：
# - intermediate/01_TumorST_preprocessed.rds.gz：Seurat 对象，assay=Spatial，含 H&E image 和 QC metadata。
# - output/01_preprocess/QC/QCData.xlsx：每个 spot 的 nCount_Spatial、nFeature_Spatial、Mito.percent。
# - output/01_preprocess/QC/Vlnplot.pdf：三个 QC 指标的小提琴图。
# - output/01_preprocess/QC/featureplot.pdf：QC 指标在组织空间位置上的分布图。
# 下游：02 读取该 Seurat 对象做形态校正聚类。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "ggplot2", "cowplot", "openxlsx", "readr"))

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
