# 05：根据 CNV、UMAP 和空间邻接关系划分 Mal / Bdy / nMal。
# 输入：intermediate/04_TumorST_cnv_scored.rds.gz，需要 CNVLabel、cnv_score、seurat_clusters、NormalScore、UMAP。
# 输出：
# - intermediate/05_TumorST_boundary_defined.rds.gz：完整对象，metadata$Location 为 Mal/Bdy/nMal，供 07/09/10/11 使用。
# - intermediate/05_TumorST_boundary_subset.rds.gz：边界扩展过程中涉及的子集对象，供复核。
# - output/05_boundary/CRC1_BoundaryDefine.pdf 以及逐轮 CRC1_out_n.pdf：边界划分图。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "magrittr", "dplyr", "purrr", "tibble", "ggplot2", "assertthat", "readr"))
source(file.path(paths$lib, "boundary_helpers.R"))

TumorST <- readr::read_rds(file.path(paths$intermediate, "04_TumorST_cnv_scored.rds.gz"))
out_dir <- file.path(paths$output, "05_boundary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# UMAP 用来衡量表达/形态相似性；空间邻接表用来限制边界只能沿相邻 spot 扩展。
UMAPembeddings <- as.data.frame(TumorST@reductions$umap@cell.embeddings)
colnames(UMAPembeddings) <- c("x", "y")
position <- load_spaceranger_positions(paths$spaceranger, rownames(TumorST@meta.data))
position$spot.ids <- seq_len(nrow(position))
dists <- compute_interspot_distances(position = position, scale.factor = 1.05)
df_j <- find_neighbors(position = position, radius = dists$radius, method = "manhattan")

# 若未手动指定恶性 CNV 标签，就从有效 Observation 标签中取 CNV score 中位数最高的两个。
# Reference/Filtered 绝不能成为恶性种子；全 NA 的低质量 spot 也不参与排序。
MalLabel <- params$malignant_cnv_labels
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
  if (length(MalLabel) < 1) stop("malignant_cnv_labels contains only forbidden Reference/Filtered labels", call. = FALSE)
}

MalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$CNVLabel %in% MalLabel, ])
NormalCluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[, c("seurat_clusters", "NormalScore")], TumorST@meta.data$seurat_clusters), function(x) mean(x$NormalScore))), decreasing = TRUE)[1]]
NormalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == NormalCluster, ])

CNV_seurat_df <- as.data.frame.array(table(TumorST@meta.data$CNVLabel, TumorST@meta.data$seurat_clusters))[MalLabel, , drop = FALSE]
ClusterID <- c()
for (cluster in levels(TumorST@meta.data$seurat_clusters)) {
  if (sum(CNV_seurat_df[, cluster]) > table(TumorST@meta.data$seurat_clusters)[cluster] * 0.5) {
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
    if (rt < 1 / 3 * rn) id
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
    data.frame(CellID = idl, Location = ifelse(rt < 1 / 3 * rn, "Mal", "Bdy"))
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
  ClusterAdd <- ClusterUpdate(x = n, position = position, df_j = df_j, UMAPembeddings = UMAPembeddings, MalCellIDN = MalCellIDN, BdyCellID = BdyCellID, NormalCellID = NormalCellID, MalCellID = MalCellID)
  Clustern <- rbind(Clustern, ClusterAdd)
  TumorSTn@meta.data$LabelNew <- Clustern$Location[match(rownames(TumorSTn@meta.data), as.character(Clustern$CellID))]
  TumorSTn@meta.data$LabelNew <- factor(TumorSTn@meta.data$LabelNew, levels = c("Normal", "Bdy", "Mal", paste0("Mal", 1:n)))

  cols_n <- c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:(n + 1)])
  pdf(file.path(out_dir, paste0(sample_name, "_out_", n, ".pdf")), width = 7, height = 7)
  print(Seurat::SpatialDimPlot(TumorSTn, group.by = "LabelNew", cols = cols_n) + ggplot2::scale_fill_manual(values = cols_n))
  dev.off()

  MalCellID <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% c("Mal", paste0("Mal", 1:n)), ])
  MalCellIDN <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% paste0("Mal", n), ])
  BdyCellID <- rownames(TumorSTn@meta.data[TumorSTn$LabelNew == "Bdy", ])
  n <- n + 1
  if (n > 6) break
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

pdf(file.path(out_dir, paste0(sample_name, "_BoundaryDefine.pdf")), width = 7, height = 7)
boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
print(Seurat::SpatialDimPlot(TumorST, group.by = "Location", cols = boundary_cols) +
  ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE))
dev.off()

readr::write_rds(TumorSTn, file.path(paths$intermediate, "05_TumorST_boundary_subset.rds.gz"), compress = "gz")
readr::write_rds(TumorST, file.path(paths$intermediate, "05_TumorST_boundary_defined.rds.gz"), compress = "gz")
