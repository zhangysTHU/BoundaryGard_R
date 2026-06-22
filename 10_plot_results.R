# 10：汇总绘图。
# 输入：
# - intermediate/05_TumorST_boundary_defined.rds.gz：Location 和空间图像坐标。
# - intermediate/07_DeconData.rds.gz：反卷积比例，用于 bar/pie 图。
# - intermediate/09_DiffGenes.rds.gz：差异表，用于火山图。
# 输出：
# - output/10_plots/DeconBarplot.pdf：各 Location 的细胞组成百分比。
# - output/10_plots/DeconPieplot.pdf：每个 spot 的细胞组成 pie，叠加低分辨率 H&E。
# - output/10_plots/DiffVolcano_<Location>.pdf：各 Location 的差异基因火山图。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "ggplot2", "ggrepel", "scatterpie", "cowplot", "jpeg", "png", "grid", "stringr", "dplyr", "tibble", "readr"))

out_dir <- file.path(paths$output, "10_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TumorST <- readr::read_rds(file.path(paths$intermediate, "05_TumorST_boundary_defined.rds.gz"))
DeconData <- readr::read_rds(file.path(paths$intermediate, "07_DeconData.rds.gz"))
DiffGenes <- readr::read_rds(file.path(paths$intermediate, "09_DiffGenes.rds.gz"))
plot_col <- colnames(DeconData)[2:ncol(DeconData)]

# 1) 各空间区域中细胞类型比例的堆叠柱状图。
metadata <- DeconData
metadata[, "Location"] <- TumorST@meta.data$Location[match(metadata$cell_ID, rownames(TumorST@meta.data))]
Location <- do.call(rbind, lapply(unique(metadata$Location), function(x) {
  metadata_split <- as.data.frame(metadata[metadata$Location == x, ])
  colSums(metadata_split[, plot_col], na.rm = TRUE) |>
    data.frame() |>
    tibble::rownames_to_column() |>
    stats::setNames(c("Types", "Sum")) |>
    dplyr::mutate(Per = 100 * Sum / sum(Sum))
}))
Location$Group <- rep(unique(metadata$Location), each = length(plot_col))
Location$Types <- factor(Location$Types, levels = plot_col)
barplot <- ggplot2::ggplot(Location, ggplot2::aes(x = Group, y = Per, fill = Types)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(out_dir, "DeconBarplot.pdf"), barplot, width = 7, height = 5)

# 2) 在低分辨率 H&E 图上叠加每个 spot 的反卷积 pie chart。
img_path <- file.path(paths$spaceranger, "spatial", "tissue_lowres_image.png")
if (file.exists(img_path)) {
  DeconData_sub <- DeconData[DeconData$cell_ID %in% rownames(TumorST@meta.data), ]
  slice <- names(TumorST@images)[1]
  image_coordinates <- tryCatch(
    data.frame(TumorST@images[[slice]]@coordinates),
    error = function(e) {
      coords <- Seurat::GetTissueCoordinates(TumorST)
      rownames(coords) <- coords$cell
      data.frame(
        row = coords$y,
        col = coords$x,
        imagerow = coords$y,
        imagecol = coords$x,
        row.names = coords$cell
      )
    }
  )
  spatial_coord <- image_coordinates |>
    tibble::rownames_to_column("cell_ID") |>
    dplyr::mutate(
      imagerow_scaled = imagerow * TumorST@images[[slice]]@scale.factors$lowres,
      imagecol_scaled = imagecol * TumorST@images[[slice]]@scale.factors$lowres
    ) |>
    dplyr::inner_join(DeconData_sub, by = "cell_ID")
  img <- png::readPNG(img_path)
  img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  pieplot <- ggplot2::ggplot() +
    ggplot2::annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
    scatterpie::geom_scatterpie(data = spatial_coord, ggplot2::aes(x = imagecol_scaled, y = imagerow_scaled), cols = plot_col, color = params$pie_border_color, alpha = params$scatterpie_alpha, pie_scale = params$pie_scale, lwd = 0.1) +
    ggplot2::scale_y_reverse() +
    ggplot2::ylim(nrow(img), 0) +
    ggplot2::xlim(0, ncol(img)) +
    ggplot2::theme_void() +
    ggplot2::coord_fixed(ratio = 1, expand = TRUE, clip = "on")
  ggplot2::ggsave(file.path(out_dir, "DeconPieplot.pdf"), pieplot, width = 7, height = 7)
}

# 3) 对每个 Location 的差异表画火山图。
for (loc in intersect(names(DiffGenes), c("Mal", "Bdy", "nMal"))) {
  LocationDiff <- DiffGenes[[loc]]
  LocationDiff <- LocationDiff[LocationDiff$Symbol %in% grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", LocationDiff$Symbol, invert = TRUE, value = TRUE), ]
  dataset <- LocationDiff |>
    tibble::rownames_to_column() |>
    stats::setNames(c("gene", colnames(LocationDiff)))
  dataset$color <- ifelse(-log10(LocationDiff$FDR) < params$volcano_p_cutoff_log10 | abs(LocationDiff$Diff) <= params$diff_logfc_cutoff, "grey", ifelse(LocationDiff$Diff > 0, loc, "Other"))
  datasetN <- dataset |>
    dplyr::group_by(color) |>
    dplyr::slice_max(order_by = abs(Diff), n = params$volcano_label_n, with_ties = FALSE) |>
    dplyr::ungroup()
  datasetN <- datasetN[datasetN$color != "grey", ]
  p <- ggplot2::ggplot(dataset, ggplot2::aes(x = Diff, y = -log10(pvalue), color = color)) +
    ggplot2::geom_point(ggplot2::aes(fill = color), size = 1) +
    ggplot2::scale_color_manual(values = c("red", "grey", "blue")) +
    ggrepel::geom_text_repel(data = datasetN, ggplot2::aes(label = Symbol), max.overlaps = 100, size = 3, segment.color = "black", show.legend = FALSE, color = "black", fontface = "bold") +
    ggplot2::geom_vline(xintercept = c(-params$diff_logfc_cutoff, params$diff_logfc_cutoff), lty = 4, col = "black", lwd = 0.8) +
    ggplot2::geom_hline(yintercept = params$volcano_p_cutoff_log10, lty = 4, col = "black", lwd = 0.8) +
    ggplot2::labs(x = "log2(FoldChange)", y = "-log10(pvalue)") +
    ggplot2::theme(panel.background = ggplot2::element_rect(color = "black", fill = NA), legend.position = "bottom", legend.title = ggplot2::element_blank())
  ggplot2::ggsave(file.path(out_dir, paste0("DiffVolcano_", loc, ".pdf")), p, width = 6, height = 5)
}
