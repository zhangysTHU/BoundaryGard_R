# 04：把 03 的 CNVLabel/cnv_score 合并回 Seurat 对象，并生成 CNV 可视化。
# 输入：
# - intermediate/02_TumorST_clustered.rds.gz：聚类后的 Seurat 对象。
# - intermediate/03_cnv_calls.tsv：优先读取，三列 cell_ID、CNVLabel、cnv_score。
# - 若 03_cnv_calls.tsv 不存在，则尝试从 output/03_infercnv/output_Spatial 的 inferCNV HMM 文件重建。
# 输出：
# - intermediate/04_TumorST_cnv_scored.rds.gz：metadata 增加 CNVLabel、cnv_score，供 05 使用。
# - output/04_cnv_score/*.pdf：CNV 空间图、UMAP 图和 CNV score 分布图。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "ape", "dendextend", "ggplot2", "ggpubr", "readr"))

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
