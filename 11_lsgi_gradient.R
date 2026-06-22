# 11：LSGI cell-component gradient 分析。
# 输入：
# - intermediate/05_TumorST_boundary_defined.rds.gz：提供 Location、空间图像和 spot 坐标。
# - intermediate/07_DeconData.rds.gz：cell_ID + 细胞类型比例，作为 LSGI embeddings。
# - ../LSGI-master/R/LSGI.R：外部 LSGI 方法源码。
# 输出：
# - intermediate/11_lsgi_cell_component_result.rds.gz：LSGI 预处理结果，可通过 --reuse-lsgi=TRUE 复用。
# - output/11_lsgi_gradient/grid_info.csv：LSGI meta-grid 信息。
# - output/11_lsgi_gradient/cell_component_gradient_arrows.csv：通过 R2 阈值筛选的细胞组分梯度箭头。
# - output/11_lsgi_gradient/cell_component_gradient_distance.csv 和 heatmap.pdf：组分梯度距离/相似性。
# - output/11_lsgi_gradient/*LSGIGradient.pdf：边界图上叠加细胞组分梯度箭头。
# - output/11_lsgi_gradient/run_summary.txt：本次参数和输出摘要。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "readr", "dplyr", "tibble", "ggplot2", "png", "grid", "viridis", "ComplexHeatmap", "reshape2"))

# 支持 --key=value 形式命令行参数，用于临时覆盖 LSGI 网格、R2 阈值、箭头样式等。
parse_cli_options <- function(args) {
  opts <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 1) {
      opts[[kv[1]]] <- TRUE
    } else {
      opts[[kv[1]]] <- paste(kv[-1], collapse = "=")
    }
  }
  opts
}

get_opt <- function(opts, dashed, underscored = gsub("-", "_", dashed), default = NULL) {
  if (!is.null(opts[[dashed]])) return(opts[[dashed]])
  if (!is.null(opts[[underscored]])) return(opts[[underscored]])
  default
}

as_bool <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

cli_opts <- parse_cli_options(commandArgs(trailingOnly = TRUE))

out_dir <- file.path(paths$output, "11_lsgi_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# LSGI 作为外部方法源码放在项目根目录的 LSGI-master 下。
lsgi_root <- normalizePath(file.path(script_dir, "..", "LSGI-master"), winslash = "/", mustWork = TRUE)
lsgi_script <- file.path(lsgi_root, "R", "LSGI.R")
if (!file.exists(lsgi_script)) {
  stop("Cannot find LSGI source script: ", lsgi_script, call. = FALSE)
}

# LSGI 内部需要 balanced_clustering；若 anticlust 不存在，就用 kmeans fallback 保证可运行。
if (requireNamespace("anticlust", quietly = TRUE)) {
  suppressPackageStartupMessages(library(anticlust))
  grid_clustering_backend <- "anticlust::balanced_clustering"
} else {
  # LSGI calls balanced_clustering() inside get.grid.coords(). This fallback keeps
  # the adapter runnable when the optional anticlust dependency is unavailable.
  balanced_clustering <- function(x, K) {
    x <- as.data.frame(x)
    K <- max(1L, min(as.integer(K), nrow(x)))
    if (K == 1L) {
      return(rep(1L, nrow(x)))
    }
    scaled_x <- scale(x)
    stats::kmeans(scaled_x, centers = K, iter.max = 100, nstart = 5)$cluster
  }
  grid_clustering_backend <- "stats::kmeans fallback"
}
source(lsgi_script)

# 读取边界结果和反卷积比例，并把 Seurat 图像坐标整理为低分辨率图上的 X/Y。
TumorST <- readr::read_rds(file.path(paths$intermediate, "05_TumorST_boundary_defined.rds.gz"))
DeconData <- readr::read_rds(file.path(paths$intermediate, "07_DeconData.rds.gz"))

slice <- names(TumorST@images)[1]
if (is.na(slice) || !nzchar(slice)) {
  stop("TumorST does not contain a spatial image slot.", call. = FALSE)
}

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

scale_factor <- TumorST@images[[slice]]@scale.factors$lowres %||% 1
spatial_df <- image_coordinates |>
  tibble::rownames_to_column("cell_ID") |>
  dplyr::mutate(
    X = imagecol * scale_factor,
    Y = imagerow * scale_factor,
    Location = TumorST@meta.data$Location[match(cell_ID, rownames(TumorST@meta.data))]
  )

decon_cols <- setdiff(colnames(DeconData), "cell_ID")
DeconData <- DeconData[DeconData$cell_ID %in% spatial_df$cell_ID, c("cell_ID", decon_cols), drop = FALSE]
DeconData <- DeconData[match(intersect(spatial_df$cell_ID, DeconData$cell_ID), DeconData$cell_ID), , drop = FALSE]
spatial_df <- spatial_df[match(DeconData$cell_ID, spatial_df$cell_ID), , drop = FALSE]

if (nrow(spatial_df) < 10 || nrow(DeconData) != nrow(spatial_df)) {
  stop("Too few matched spots between boundary object and deconvolution matrix.", call. = FALSE)
}

# LSGI 输入：spatial_coords 是 spot 的二维空间位置，embeddings 是每个 spot 的细胞组分比例。
spatial_coords <- spatial_df[, c("X", "Y"), drop = FALSE]
rownames(spatial_coords) <- spatial_df$cell_ID

embeddings <- as.matrix(DeconData[, decon_cols, drop = FALSE])
storage.mode(embeddings) <- "numeric"
rownames(embeddings) <- DeconData$cell_ID
embeddings <- embeddings[, colSums(abs(embeddings), na.rm = TRUE) > 0, drop = FALSE]
if (ncol(embeddings) < 1) {
  stop("No non-zero cell-component columns found in 07_DeconData.rds.gz.", call. = FALSE)
}

# LSGI 参数可来自命令行，也可来自 00_config.R 的 params$lsgi_*，否则使用默认值。
n_grids_scale <- as.numeric(get_opt(cli_opts, "n-grids-scale", default = params$lsgi_n_grids_scale %||% 10))
n_cells_per_meta <- as.numeric(get_opt(cli_opts, "n-cells-per-meta", default = params$lsgi_n_cells_per_meta %||% min(50, nrow(spatial_coords))))
r_squared_thresh <- as.numeric(get_opt(cli_opts, "r-squared-thresh", default = params$lsgi_r_squared_thresh %||% 0.3))
minimum_fctr <- as.numeric(get_opt(cli_opts, "minimum-fctr", default = params$lsgi_minimum_fctr %||% 3))
arrow_length_scale <- as.numeric(get_opt(cli_opts, "arrow-length-scale", default = params$lsgi_arrow_length_scale %||% 1.4))
arrow_linewidth <- as.numeric(get_opt(cli_opts, "arrow-linewidth", default = params$lsgi_arrow_linewidth %||% 1.0))
arrow_head_cm <- as.numeric(get_opt(cli_opts, "arrow-head-cm", default = params$lsgi_arrow_head_cm %||% 0.20))
arrow_closed <- as_bool(get_opt(cli_opts, "arrow-closed", default = params$lsgi_arrow_closed %||% TRUE))
reuse_lsgi <- as_bool(get_opt(cli_opts, "reuse-lsgi", default = FALSE))
arrow_type <- if (arrow_closed) "closed" else "open"

# local.traj.preprocessing 是 LSGI 的耗时预处理，可缓存复用。
lsgi_result_path <- file.path(paths$intermediate, "11_lsgi_cell_component_result.rds.gz")
if (reuse_lsgi && file.exists(lsgi_result_path)) {
  lsgi_res <- readr::read_rds(lsgi_result_path)
} else {
  lsgi_res <- local.traj.preprocessing(
    spatial_coords = spatial_coords,
    embeddings = embeddings,
    n.grids.scale = n_grids_scale,
    n.cells.per.meta = n_cells_per_meta
  )
  readr::write_rds(lsgi_res, lsgi_result_path, compress = "gz")
}
utils::write.csv(lsgi_res$grid.info, file.path(out_dir, "grid_info.csv"), row.names = FALSE)

# get.ind.rsqrs 估计每个局部网格中各细胞组分梯度方向，并用 R2 衡量线性趋势强度。
lin_res <- get.ind.rsqrs(lsgi_res)
lin_res <- stats::na.omit(lin_res)
arrow_df <- lin_res[lin_res$rsquared > r_squared_thresh, , drop = FALSE]
if (nrow(arrow_df) > 0) {
  arrow_df <- arrow_df |>
    dplyr::group_by(fctr) |>
    dplyr::filter(dplyr::n() >= minimum_fctr) |>
    dplyr::ungroup() |>
    as.data.frame()
}
utils::write.csv(arrow_df, file.path(out_dir, "cell_component_gradient_arrows.csv"), row.names = FALSE)

# 组分梯度之间的距离/相似性矩阵，用热图展示不同细胞组分空间变化是否同向。
dist_mat <- tryCatch(
  avg.dist.calc(lsgi_res, r_squared_thresh = r_squared_thresh, minimum.fctr = minimum_fctr),
  error = function(e) {
    warning("LSGI distance calculation skipped: ", conditionMessage(e), call. = FALSE)
    NULL
  }
)
if (!is.null(dist_mat) && nrow(dist_mat) > 0) {
  utils::write.csv(dist_mat, file.path(out_dir, "cell_component_gradient_distance.csv"), row.names = FALSE)
  grDevices::pdf(file.path(out_dir, "cell_component_gradient_distance_heatmap.pdf"), width = 7, height = 6)
  ComplexHeatmap::draw(plt.dist.heat(dist_mat))
  grDevices::dev.off()
}

boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
point_df <- spatial_df |>
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))

arrow_plot_df <- arrow_df
if (nrow(arrow_plot_df) > 0) {
  arrow_plot_df <- arrow_plot_df |>
    dplyr::mutate(
      X_end = X + vx.u * arrow_length_scale,
      Y_end = Y + vy.u * arrow_length_scale
    )
}

# 把筛选后的 LSGI 梯度箭头叠加到已有 ggplot 边界图上。
add_gradient_arrows <- function(p, arrow_data = arrow_plot_df) {
  if (nrow(arrow_data) == 0) {
    return(p + ggplot2::labs(subtitle = paste0("No component gradients passed R2 > ", r_squared_thresh)))
  }
  p +
    ggplot2::geom_segment(
      data = arrow_data,
      ggplot2::aes(
        x = X,
        y = Y,
        xend = X_end,
        yend = Y_end,
        color = fctr
      ),
      inherit.aes = FALSE,
      linewidth = arrow_linewidth,
      lineend = "round",
      arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), type = arrow_type)
    ) +
    ggplot2::labs(color = "LSGI component")
}

# 1) 不带 H&E 背景的边界 + 梯度箭头图。
boundary_base <- ggplot2::ggplot(point_df, ggplot2::aes(x = X, y = Y, fill = Location)) +
  ggplot2::geom_point(shape = 21, size = 1.8, stroke = 0.1, color = "grey25", alpha = 0.9) +
  ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE) +
  ggplot2::scale_y_reverse() +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::ggtitle(paste0(sample_name, " boundary with LSGI cell-component gradients"))

boundary_gradient <- add_gradient_arrows(boundary_base)
ggplot2::ggsave(
  file.path(out_dir, paste0(sample_name, "_BoundaryDefine_LSGIGradient.pdf")),
  boundary_gradient,
  width = 8,
  height = 7
)

# 2) H&E 背景上的边界 + 梯度箭头图。
img_path <- file.path(paths$spaceranger, "spatial", "tissue_lowres_image.png")
if (file.exists(img_path)) {
  img <- png::readPNG(img_path)
  img_grob <- grid::rasterGrob(
    img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc")
  )
  point_df_he <- point_df |>
    dplyr::mutate(Y = -Y)
  arrow_plot_df_he <- arrow_plot_df
  if (nrow(arrow_plot_df_he) > 0) {
    arrow_plot_df_he <- arrow_plot_df_he |>
      dplyr::mutate(
        Y = -Y,
        Y_end = -Y_end
      )
  }
  he_base <- ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = -nrow(img),
      ymax = 0
    ) +
    ggplot2::geom_point(
      data = point_df_he,
      ggplot2::aes(x = X, y = Y, fill = Location),
      shape = 21,
      size = 1.8,
      stroke = 0.1,
      color = "grey20",
      alpha = 0.82
    ) +
    ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, ncol(img)),
      ylim = c(-nrow(img), 0),
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::ggtitle(paste0(sample_name, " HE-boundary with LSGI cell-component gradients"))

  he_gradient <- add_gradient_arrows(he_base, arrow_plot_df_he)
  ggplot2::ggsave(
    file.path(out_dir, paste0(sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf")),
    he_gradient,
    width = 8,
    height = 7
  )
}

# 3) LSGI 自带的纯梯度图，便于和 adapter 绘图对照。
grDevices::pdf(file.path(out_dir, "cell_component_gradients_plain_lsgi.pdf"), width = 8, height = 7)
print(plt.factors.gradient.ind(
  info = lsgi_res,
  r_squared_thresh = r_squared_thresh,
  minimum.fctr = minimum_fctr,
  arrow.length.scale = arrow_length_scale
) + ggplot2::ggtitle("LSGI cell-component gradients"))
grDevices::dev.off()

# 记录本次参数、匹配 spot 数和主要输出，方便复现实验设置。
sink(file.path(out_dir, "run_summary.txt"))
cat("Cottrazm + LSGI cell-component gradient analysis\n")
cat("================================================\n\n")
cat("Sample:", sample_name, "\n")
cat("Matched spots:", nrow(spatial_coords), "\n")
cat("Cell components:", paste(colnames(embeddings), collapse = ", "), "\n")
cat("Grid clustering backend:", grid_clustering_backend, "\n")
cat("n.grids.scale:", n_grids_scale, "\n")
cat("n.cells.per.meta:", n_cells_per_meta, "\n")
cat("R-squared threshold:", r_squared_thresh, "\n")
cat("Minimum arrows per component:", minimum_fctr, "\n")
cat("Selected gradient arrows:", nrow(arrow_df), "\n\n")
cat("Arrow length scale:", arrow_length_scale, "\n")
cat("Arrow linewidth:", arrow_linewidth, "\n")
cat("Arrow head cm:", arrow_head_cm, "\n\n")
cat("Arrow closed:", arrow_closed, "\n")
cat("Reused LSGI result:", reuse_lsgi && file.exists(lsgi_result_path), "\n\n")
cat("Outputs:\n")
cat("- intermediate/11_lsgi_cell_component_result.rds.gz\n")
cat("- output/11_lsgi_gradient/grid_info.csv\n")
cat("- output/11_lsgi_gradient/cell_component_gradient_arrows.csv\n")
cat("- output/11_lsgi_gradient/", sample_name, "_BoundaryDefine_LSGIGradient.pdf\n", sep = "")
cat("- output/11_lsgi_gradient/", sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf\n", sep = "")
sink()

message("Done. LSGI gradient outputs written to ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), ".")
