# 11b：LSGI cell-component gradient 分析的增强副本。
# 在原始 11_lsgi_gradient.R 基础上，额外输出：
# - grid_local_spots.pdf / grid_local_spots_overlay.pdf：每个 grid 实际用于局部回归的 spot 分页图。
# - grid_partition.pdf / grid_partition_overlay.pdf：每个 spot 唯一归属最近 grid 的分区图。
# - spot_grid_membership.csv：按 spot 汇总其参与的 local grids 和唯一最近 grid。
# - grid_spot_summary.csv：按 grid 汇总其 local spots 和 partition spots。
# - grid_local_spot_membership.csv / grid_partition_membership.csv：grid-spot 明细表。
# - output/11_lsgi_gradient/arrow_tables/cell_component_arrows_by_grid.csv：按 grid 排列的完整 arrow 明细。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "readr", "dplyr", "tibble", "ggplot2", "png", "grid", "viridis", "ComplexHeatmap", "reshape2", "magrittr"))

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

lsgi_root_candidates <- c(
  file.path(script_dir, "..", "LSGI-master"),
  file.path(script_dir, "..", "Cottrazm-main", "LSGI-master"),
  file.path(script_dir, "..", "Cottrazm-main-archived", "LSGI-master")
)
lsgi_root <- NULL
for (candidate in lsgi_root_candidates) {
  candidate_norm <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
  if (dir.exists(candidate_norm) && file.exists(file.path(candidate_norm, "R", "LSGI.R"))) {
    lsgi_root <- candidate_norm
    break
  }
}
if (is.null(lsgi_root)) {
  stop(
    "Cannot find LSGI source directory. Tried: ",
    paste(normalizePath(lsgi_root_candidates, winslash = "/", mustWork = FALSE), collapse = ", "),
    call. = FALSE
  )
}
lsgi_script <- file.path(lsgi_root, "R", "LSGI.R")

if (requireNamespace("anticlust", quietly = TRUE)) {
  suppressPackageStartupMessages(library(anticlust))
  grid_clustering_backend <- paste0("anticlust::balanced_clustering ", as.character(utils::packageVersion("anticlust")))
} else {
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
      row = coords$x,
      col = coords$y,
      imagerow = coords$x,
      imagecol = coords$y,
      row.names = coords$cell
    )
  }
)

scale_factor <- TumorST@images[[slice]]@scale.factors$lowres %||% 1
spatial_df <- image_coordinates %>%
  tibble::rownames_to_column("cell_ID") %>%
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

spatial_coords <- spatial_df[, c("X", "Y"), drop = FALSE]
rownames(spatial_coords) <- spatial_df$cell_ID

embeddings <- as.matrix(DeconData[, decon_cols, drop = FALSE])
storage.mode(embeddings) <- "numeric"
rownames(embeddings) <- DeconData$cell_ID
embeddings <- embeddings[, colSums(abs(embeddings), na.rm = TRUE) > 0, drop = FALSE]
if (ncol(embeddings) < 1) {
  stop("No non-zero cell-component columns found in 07_DeconData.rds.gz.", call. = FALSE)
}

n_grids_scale <- as.numeric(get_opt(cli_opts, "n-grids-scale", default = params$lsgi_n_grids_scale %||% 10))
n_cells_per_meta <- as.numeric(get_opt(cli_opts, "n-cells-per-meta", default = params$lsgi_n_cells_per_meta %||% min(50, nrow(spatial_coords))))
r_squared_thresh <- as.numeric(get_opt(cli_opts, "r-squared-thresh", default = params$lsgi_r_squared_thresh %||% 0.3))
minimum_fctr <- as.numeric(get_opt(cli_opts, "minimum-fctr", default = params$lsgi_minimum_fctr %||% 3))
arrow_length_scale <- as.numeric(get_opt(cli_opts, "arrow-length-scale", default = params$lsgi_arrow_length_scale %||% 1.4))
arrow_linewidth <- as.numeric(get_opt(cli_opts, "arrow-linewidth", default = params$lsgi_arrow_linewidth %||% 1.0))
arrow_head_cm <- as.numeric(get_opt(cli_opts, "arrow-head-cm", default = params$lsgi_arrow_head_cm %||% 0.20))
arrow_head_angle <- as.numeric(get_opt(cli_opts, "arrow-head-angle", default = params$lsgi_arrow_head_angle %||% 30))
arrow_head_angle_by_r2 <- as_bool(get_opt(cli_opts, "arrow-head-angle-by-r2", default = params$lsgi_arrow_head_angle_by_r2 %||% TRUE))
arrow_head_angle_min <- as.numeric(get_opt(cli_opts, "arrow-head-angle-min", default = params$lsgi_arrow_head_angle_min %||% 22.5))
arrow_head_angle_mid <- as.numeric(get_opt(cli_opts, "arrow-head-angle-mid", default = params$lsgi_arrow_head_angle_mid %||% arrow_head_angle))
arrow_head_angle_max <- as.numeric(get_opt(cli_opts, "arrow-head-angle-max", default = params$lsgi_arrow_head_angle_max %||% 45))
arrow_head_angle_r2_min <- as.numeric(get_opt(cli_opts, "arrow-head-angle-r2-min", default = params$lsgi_arrow_head_angle_r2_min %||% 0.3))
arrow_head_angle_r2_mid <- as.numeric(get_opt(cli_opts, "arrow-head-angle-r2-mid", default = params$lsgi_arrow_head_angle_r2_mid %||% 0.5))
arrow_head_angle_r2_max <- as.numeric(get_opt(cli_opts, "arrow-head-angle-r2-max", default = params$lsgi_arrow_head_angle_r2_max %||% 0.7))
arrow_head_angle_step <- as.numeric(get_opt(cli_opts, "arrow-head-angle-step", default = params$lsgi_arrow_head_angle_step %||% 1))
arrow_closed <- as_bool(get_opt(cli_opts, "arrow-closed", default = params$lsgi_arrow_closed %||% TRUE))
arrow_length_normalization <- as.character(get_opt(
  cli_opts,
  "arrow-length-normalization",
  default = params$lsgi_arrow_length_normalization %||% "global"
))
reuse_lsgi <- as_bool(get_opt(cli_opts, "reuse-lsgi", default = FALSE))
arrow_type <- if (arrow_closed) "closed" else "open"
if (!arrow_length_normalization %in% c("global", "by_component")) {
  stop("--arrow-length-normalization must be either 'global' or 'by_component'.", call. = FALSE)
}
if (!is.finite(arrow_head_angle) || arrow_head_angle <= 0 || arrow_head_angle >= 180) {
  stop("--arrow-head-angle must be a finite value between 0 and 180 degrees.", call. = FALSE)
}
if (arrow_head_angle_by_r2) {
  if (
    !all(is.finite(c(arrow_head_angle_min, arrow_head_angle_mid, arrow_head_angle_max))) ||
      arrow_head_angle_min <= 0 ||
      arrow_head_angle_min > arrow_head_angle_mid ||
      arrow_head_angle_mid > arrow_head_angle_max ||
      arrow_head_angle_max >= 180
  ) {
    stop("R2-mapped arrow head angles must satisfy 0 < min <= mid <= max < 180.", call. = FALSE)
  }
  if (
    !all(is.finite(c(arrow_head_angle_r2_min, arrow_head_angle_r2_mid, arrow_head_angle_r2_max))) ||
      arrow_head_angle_r2_min >= arrow_head_angle_r2_mid ||
      arrow_head_angle_r2_mid >= arrow_head_angle_r2_max
  ) {
    stop("R2 anchors for arrow head angle must satisfy min < mid < max.", call. = FALSE)
  }
  if (!is.finite(arrow_head_angle_step) || arrow_head_angle_step <= 0) {
    stop("--arrow-head-angle-step must be a finite positive value.", call. = FALSE)
  }
}

lsgi_result_path <- file.path(paths$intermediate, "11_lsgi_cell_component_result.rds.gz")
did_reuse_lsgi <- FALSE
if (reuse_lsgi && file.exists(lsgi_result_path)) {
  lsgi_res <- readr::read_rds(lsgi_result_path)
  did_reuse_lsgi <- TRUE
} else {
  lsgi_res <- local.traj.preprocessing(
    spatial_coords = spatial_coords,
    embeddings = embeddings,
    n.grids.scale = n_grids_scale,
    n.cells.per.meta = n_cells_per_meta
  )
  readr::write_rds(lsgi_res, lsgi_result_path, compress = "gz")
}

grid_ids <- paste0("grid_", seq_len(nrow(lsgi_res$grid.info)))
dist_to_grid <- as.matrix(lsgi_res$dist.to.grid)
rownames(dist_to_grid) <- rownames(spatial_coords)
colnames(dist_to_grid) <- grid_ids

spot_base_df <- spatial_df %>%
  dplyr::left_join(DeconData, by = "cell_ID")

grid_info_df <- lsgi_res$grid.info %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    grid = grid_ids,
    grid_index = seq_len(dplyr::n())
  ) %>%
  dplyr::select(grid, grid_index, X, Y, vx, vy, R_squared, Assignment, qsum, sf, vx.u, vy.u)
utils::write.csv(grid_info_df, file.path(out_dir, "grid_info.csv"), row.names = FALSE)

grid_palette <- setNames(
  grDevices::rainbow(length(grid_ids), s = 0.7, v = 0.9, end = if (length(grid_ids) > 1) 0.95 else 0.01),
  grid_ids
)

grid_centers_df <- grid_info_df %>%
  dplyr::transmute(
    grid,
    grid_center_X = X,
    grid_center_Y = Y
  )

local_membership_list <- vector("list", length(grid_ids))
for (i in seq_along(grid_ids)) {
  cells <- lsgi_res$local.linear.info$cell[[i]]
  local_membership_list[[i]] <- tibble::tibble(
    grid = grid_ids[i],
    cell_ID = cells,
    local_rank = seq_along(cells),
    dist_to_grid = as.numeric(dist_to_grid[cells, i])
  )
}

local_membership_df <- dplyr::bind_rows(local_membership_list) %>%
  dplyr::left_join(spot_base_df[, c("cell_ID", "X", "Y", "Location"), drop = FALSE], by = "cell_ID") %>%
  dplyr::left_join(grid_centers_df, by = "grid") %>%
  dplyr::arrange(grid, local_rank)
utils::write.csv(local_membership_df, file.path(out_dir, "grid_local_spot_membership.csv"), row.names = FALSE)

nearest_idx <- max.col(-dist_to_grid, ties.method = "first")
partition_df <- tibble::tibble(
  cell_ID = rownames(dist_to_grid),
  nearest_grid = grid_ids[nearest_idx],
  nearest_grid_distance = dist_to_grid[cbind(seq_len(nrow(dist_to_grid)), nearest_idx)]
) %>%
  dplyr::left_join(spot_base_df[, c("cell_ID", "X", "Y", "Location"), drop = FALSE], by = "cell_ID") %>%
  dplyr::left_join(grid_centers_df, by = c("nearest_grid" = "grid")) %>%
  dplyr::arrange(nearest_grid, cell_ID)
utils::write.csv(partition_df, file.path(out_dir, "grid_partition_membership.csv"), row.names = FALSE)

spot_membership_summary <- local_membership_df %>%
  dplyr::group_by(cell_ID) %>%
  dplyr::summarise(
    local_grid_count = dplyr::n(),
    local_grids = paste(grid, collapse = ";"),
    local_grid_ranks = paste(local_rank, collapse = ";"),
    .groups = "drop"
  )

spot_table_df <- spot_base_df %>%
  dplyr::left_join(spot_membership_summary, by = "cell_ID") %>%
  dplyr::left_join(partition_df[, c("cell_ID", "nearest_grid", "nearest_grid_distance", "grid_center_X", "grid_center_Y"), drop = FALSE], by = "cell_ID") %>%
  dplyr::mutate(
    local_grid_count = ifelse(is.na(local_grid_count), 0L, local_grid_count),
    local_grids = dplyr::coalesce(local_grids, ""),
    local_grid_ranks = dplyr::coalesce(local_grid_ranks, "")
  )
utils::write.csv(spot_table_df, file.path(out_dir, "spot_grid_membership.csv"), row.names = FALSE)

grid_local_summary <- local_membership_df %>%
  dplyr::group_by(grid) %>%
  dplyr::summarise(
    n_local_spots = dplyr::n(),
    local_spots = paste(cell_ID, collapse = ";"),
    .groups = "drop"
  )

grid_partition_summary <- partition_df %>%
  dplyr::group_by(nearest_grid) %>%
  dplyr::summarise(
    n_partition_spots = dplyr::n(),
    partition_spots = paste(cell_ID, collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::rename(grid = nearest_grid)

grid_table_df <- grid_info_df %>%
  dplyr::left_join(grid_local_summary, by = "grid") %>%
  dplyr::left_join(grid_partition_summary, by = "grid") %>%
  dplyr::mutate(
    n_local_spots = ifelse(is.na(n_local_spots), 0L, n_local_spots),
    local_spots = dplyr::coalesce(local_spots, ""),
    n_partition_spots = ifelse(is.na(n_partition_spots), 0L, n_partition_spots),
    partition_spots = dplyr::coalesce(partition_spots, "")
  )
utils::write.csv(grid_table_df, file.path(out_dir, "grid_spot_summary.csv"), row.names = FALSE)

lin_res <- get.ind.rsqrs(lsgi_res)
lin_res <- stats::na.omit(lin_res)

rescale_arrow_length <- function(df) {
  df <- df %>%
    dplyr::mutate(
      gradient_strength = sqrt(vx^2 + vy^2),
      original_vx.u = vx.u,
      original_vy.u = vy.u,
      original_scaled_length = sqrt(vx.u^2 + vy.u^2)
    )

  scale_one_group <- function(group_df) {
    finite_strength <- is.finite(group_df$gradient_strength)
    if (!any(finite_strength)) {
      group_df$gradient_strength_norm <- NA_real_
      group_df$arrow_length_multiplier <- NA_real_
      return(group_df)
    }

    min_strength <- min(group_df$gradient_strength[finite_strength], na.rm = TRUE)
    max_strength <- max(group_df$gradient_strength[finite_strength], na.rm = TRUE)
    if (!is.finite(min_strength) || !is.finite(max_strength) || max_strength <= min_strength) {
      strength_norm <- ifelse(finite_strength, 0.5, NA_real_)
    } else {
      strength_norm <- (group_df$gradient_strength - min_strength) / (max_strength - min_strength)
      strength_norm[!finite_strength] <- NA_real_
    }

    group_df$gradient_strength_norm <- strength_norm
    group_df$arrow_length_multiplier <- 0.5 + 1.5 * strength_norm
    group_df
  }

  df <- if (identical(arrow_length_normalization, "global")) {
    scale_one_group(df)
  } else {
    df %>%
      dplyr::group_by(fctr) %>%
      dplyr::group_modify(~ scale_one_group(.x)) %>%
      dplyr::ungroup()
  }

  df %>%
    dplyr::mutate(
      vx.u = original_vx.u * arrow_length_multiplier,
      vy.u = original_vy.u * arrow_length_multiplier,
      scaled_length = sqrt(vx.u^2 + vy.u^2),
      plotted_length = scaled_length * arrow_length_scale
    )
}

lin_res <- rescale_arrow_length(lin_res)

map_arrow_head_angle <- function(rsquared) {
  if (!arrow_head_angle_by_r2) {
    return(rep(arrow_head_angle, length(rsquared)))
  }

  r2_clipped <- pmin(pmax(rsquared, arrow_head_angle_r2_min), arrow_head_angle_r2_max)
  angle <- ifelse(
    r2_clipped <= arrow_head_angle_r2_mid,
    arrow_head_angle_min +
      (r2_clipped - arrow_head_angle_r2_min) /
        (arrow_head_angle_r2_mid - arrow_head_angle_r2_min) *
        (arrow_head_angle_mid - arrow_head_angle_min),
    arrow_head_angle_mid +
      (r2_clipped - arrow_head_angle_r2_mid) /
        (arrow_head_angle_r2_max - arrow_head_angle_r2_mid) *
        (arrow_head_angle_max - arrow_head_angle_mid)
  )
  angle <- round(angle / arrow_head_angle_step) * arrow_head_angle_step
  pmin(pmax(angle, arrow_head_angle_min), arrow_head_angle_max)
}

lin_res <- lin_res %>%
  dplyr::mutate(arrow_head_angle_mapped = map_arrow_head_angle(rsquared))

pass_count_df <- lin_res %>%
  dplyr::filter(rsquared > r_squared_thresh) %>%
  dplyr::count(fctr, name = "n_grids_passing_r_squared")

all_arrow_df <- lin_res %>%
  dplyr::left_join(pass_count_df, by = "fctr") %>%
  dplyr::mutate(
    n_grids_passing_r_squared = ifelse(is.na(n_grids_passing_r_squared), 0L, n_grids_passing_r_squared),
    component = fctr,
    raw_length = gradient_strength,
    passes_r_squared = rsquared > r_squared_thresh,
    passes_minimum_fctr = n_grids_passing_r_squared >= minimum_fctr,
    included_in_filtered_output = passes_r_squared & passes_minimum_fctr
  ) %>%
  dplyr::arrange(grid, component)

grid_arrow_df <- all_arrow_df %>%
  dplyr::transmute(
    grid,
    Assignment = as.character(component),
    gradient_strength,
    gradient_strength_norm,
    arrow_length_multiplier,
    arrow_head_angle_mapped,
    vx.u_rescaled = vx.u,
    vy.u_rescaled = vy.u,
    scaled_length_rescaled = scaled_length,
    plotted_length
  )
grid_info_df <- grid_info_df %>%
  dplyr::mutate(
    Assignment = as.character(Assignment),
    original_vx.u = vx.u,
    original_vy.u = vy.u,
    original_scaled_length = sqrt(vx.u^2 + vy.u^2)
  ) %>%
  dplyr::left_join(grid_arrow_df, by = c("grid", "Assignment"), suffix = c("", "_all_component")) %>%
  dplyr::mutate(
    vx.u = dplyr::coalesce(vx.u_rescaled, vx.u),
    vy.u = dplyr::coalesce(vy.u_rescaled, vy.u),
    scaled_length = dplyr::coalesce(scaled_length_rescaled, sqrt(vx.u^2 + vy.u^2))
  ) %>%
  dplyr::select(-vx.u_rescaled, -vy.u_rescaled, -scaled_length_rescaled)
utils::write.csv(grid_info_df, file.path(out_dir, "grid_info.csv"), row.names = FALSE)

arrow_dir <- file.path(out_dir, "arrow_tables")
dir.create(arrow_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  all_arrow_df[, c(
    "grid", "component", "X", "Y", "vx", "vy", "rsquared", "raw_length",
    "gradient_strength", "gradient_strength_norm", "arrow_length_multiplier",
    "arrow_head_angle_mapped",
    "sf", "original_vx.u", "original_vy.u", "original_scaled_length",
    "vx.u", "vy.u", "scaled_length", "plotted_length",
    "passes_r_squared", "passes_minimum_fctr", "included_in_filtered_output"
  )],
  file.path(arrow_dir, "cell_component_arrows_by_grid.csv"),
  row.names = FALSE
)

arrow_df <- all_arrow_df[all_arrow_df$included_in_filtered_output, c(
  "vx", "vy", "rsquared", "component", "grid", "X", "Y", "qsum", "sf",
  "gradient_strength", "gradient_strength_norm", "arrow_length_multiplier",
  "arrow_head_angle_mapped",
  "original_vx.u", "original_vy.u", "original_scaled_length",
  "vx.u", "vy.u", "scaled_length", "plotted_length"
), drop = FALSE]
colnames(arrow_df)[colnames(arrow_df) == "component"] <- "fctr"
utils::write.csv(arrow_df, file.path(out_dir, "cell_component_gradient_arrows.csv"), row.names = FALSE)

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
point_df <- spot_base_df %>%
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))

arrow_plot_df <- arrow_df
if (nrow(arrow_plot_df) > 0) {
  arrow_plot_df <- arrow_plot_df %>%
    dplyr::mutate(
      X_end = X + vx.u * arrow_length_scale,
      Y_end = Y + vy.u * arrow_length_scale
    )
}

add_gradient_arrows <- function(p, arrow_data = arrow_plot_df) {
  if (nrow(arrow_data) == 0) {
    return(p + ggplot2::labs(subtitle = paste0("No component gradients passed R2 > ", r_squared_thresh)))
  }
  if (!"arrow_head_angle_mapped" %in% colnames(arrow_data)) {
    arrow_data$arrow_head_angle_mapped <- arrow_head_angle
  }
  angle_values <- sort(unique(arrow_data$arrow_head_angle_mapped[is.finite(arrow_data$arrow_head_angle_mapped)]))
  if (length(angle_values) == 0) {
    angle_values <- arrow_head_angle
  }

  out <- p
  for (angle_value in angle_values) {
    layer_data <- arrow_data[is.finite(arrow_data$arrow_head_angle_mapped) & arrow_data$arrow_head_angle_mapped == angle_value, , drop = FALSE]
    if (nrow(layer_data) == 0) next
    out <- out +
      ggplot2::geom_segment(
        data = layer_data,
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
        arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), angle = angle_value, type = arrow_type)
      )
  }
  out + ggplot2::labs(color = "LSGI component")
}

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

img_path <- file.path(paths$spaceranger, "spatial", "tissue_lowres_image.png")
has_image <- file.exists(img_path)
if (has_image) {
  img <- png::readPNG(img_path)
  img_grob <- grid::rasterGrob(
    img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc")
  )
  point_df_he <- point_df %>%
    dplyr::mutate(Y = -Y)
  arrow_plot_df_he <- arrow_plot_df
  if (nrow(arrow_plot_df_he) > 0) {
    arrow_plot_df_he <- arrow_plot_df_he %>%
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

plain_gradient_base <- ggplot2::ggplot(point_df, ggplot2::aes(x = X, y = Y)) +
  ggplot2::geom_point(size = 1.5, shape = 20, stroke = 0, color = "lightgrey") +
  ggplot2::scale_y_reverse() +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::ggtitle("LSGI cell-component gradients")
ggplot2::ggsave(
  file.path(out_dir, "cell_component_gradients_plain_lsgi.pdf"),
  add_gradient_arrows(plain_gradient_base),
  width = 8,
  height = 7
)

partition_plot <- ggplot2::ggplot(partition_df, ggplot2::aes(x = X, y = Y, color = nearest_grid)) +
  ggplot2::geom_point(size = 1.6, alpha = 0.9) +
  ggplot2::geom_point(
    data = grid_info_df,
    ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    shape = 4,
    size = 1.8,
    stroke = 0.8,
    color = "black"
  ) +
  ggplot2::scale_color_manual(values = grid_palette, guide = "none") +
  ggplot2::scale_y_reverse() +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::ggtitle(paste0(sample_name, " grid partition (nearest-grid assignment)"))
ggplot2::ggsave(file.path(out_dir, "grid_partition.pdf"), partition_plot, width = 8, height = 7)

if (has_image) {
  partition_df_he <- partition_df %>%
    dplyr::mutate(Y = -Y)
  grid_info_he <- grid_info_df %>%
    dplyr::mutate(Y = -Y)
  partition_overlay_plot <- ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = -nrow(img),
      ymax = 0
    ) +
    ggplot2::geom_point(
      data = partition_df_he,
      ggplot2::aes(x = X, y = Y, color = nearest_grid),
      size = 1.6,
      alpha = 0.88
    ) +
    ggplot2::geom_point(
      data = grid_info_he,
      ggplot2::aes(x = X, y = Y),
      inherit.aes = FALSE,
      shape = 4,
      size = 1.8,
      stroke = 0.8,
      color = "black"
    ) +
    ggplot2::scale_color_manual(values = grid_palette, guide = "none") +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, ncol(img)),
      ylim = c(-nrow(img), 0),
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste0(sample_name, " grid partition overlay"))
  ggplot2::ggsave(file.path(out_dir, "grid_partition_overlay.pdf"), partition_overlay_plot, width = 8, height = 7)
} else {
  ggplot2::ggsave(file.path(out_dir, "grid_partition_overlay.pdf"), partition_plot, width = 8, height = 7)
}

make_local_plot <- function(grid_id, overlay = FALSE) {
  local_df <- local_membership_df[local_membership_df$grid == grid_id, , drop = FALSE]
  center_df <- grid_info_df[grid_info_df$grid == grid_id, , drop = FALSE]
  grid_color <- grid_palette[[grid_id]]
  title_text <- paste0(sample_name, " ", grid_id, " local spots (n=", nrow(local_df), ", component=", as.character(center_df$Assignment[1]), ")")

  if (!overlay || !has_image) {
    return(
      ggplot2::ggplot() +
        ggplot2::geom_point(
          data = point_df,
          ggplot2::aes(x = X, y = Y),
          color = "grey82",
          size = 1.3,
          alpha = 0.65
        ) +
        ggplot2::geom_point(
          data = local_df,
          ggplot2::aes(x = X, y = Y),
          color = grid_color,
          size = 1.9,
          alpha = 0.95
        ) +
        ggplot2::geom_point(
          data = center_df,
          ggplot2::aes(x = X, y = Y),
          inherit.aes = FALSE,
          shape = 4,
          size = 2.3,
          stroke = 0.9,
          color = "black"
        ) +
        ggplot2::scale_y_reverse() +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::ggtitle(if (overlay && !has_image) paste0(title_text, " overlay") else title_text)
    )
  }

  local_df_he <- local_df %>%
    dplyr::mutate(Y = -Y)
  center_df_he <- center_df %>%
    dplyr::mutate(Y = -Y)
  point_df_he <- point_df %>%
    dplyr::mutate(Y = -Y)

  ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = -nrow(img),
      ymax = 0
    ) +
    ggplot2::geom_point(
      data = point_df_he,
      ggplot2::aes(x = X, y = Y),
      color = "grey82",
      size = 1.2,
      alpha = 0.55
    ) +
    ggplot2::geom_point(
      data = local_df_he,
      ggplot2::aes(x = X, y = Y),
      color = grid_color,
      size = 1.9,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      data = center_df_he,
      ggplot2::aes(x = X, y = Y),
      inherit.aes = FALSE,
      shape = 4,
      size = 2.3,
      stroke = 0.9,
      color = "black"
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, ncol(img)),
      ylim = c(-nrow(img), 0),
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste0(title_text, " overlay"))
}

grDevices::pdf(file.path(out_dir, "grid_local_spots.pdf"), width = 8, height = 7)
for (grid_id in grid_ids) {
  print(make_local_plot(grid_id, overlay = FALSE))
}
grDevices::dev.off()

grDevices::pdf(file.path(out_dir, "grid_local_spots_overlay.pdf"), width = 8, height = 7)
for (grid_id in grid_ids) {
  print(make_local_plot(grid_id, overlay = TRUE))
}
grDevices::dev.off()

sink(file.path(out_dir, "run_summary_11b.txt"))
cat("Cottrazm + LSGI detailed grid-membership analysis\n")
cat("=================================================\n\n")
cat("Sample:", sample_name, "\n")
cat("Matched spots:", nrow(spatial_coords), "\n")
cat("Cell components:", paste(colnames(embeddings), collapse = ", "), "\n")
cat("Grid clustering backend:", grid_clustering_backend, "\n")
cat("n.grids.scale:", n_grids_scale, "\n")
cat("n.cells.per.meta:", n_cells_per_meta, "\n")
cat("R-squared threshold:", r_squared_thresh, "\n")
cat("Minimum arrows per component:", minimum_fctr, "\n")
cat("Selected gradient arrows:", nrow(arrow_df), "\n")
cat("Generated grids:", nrow(grid_info_df), "\n")
cat("H&E overlay available:", has_image, "\n\n")
cat("Arrow length scale:", arrow_length_scale, "\n")
cat("Arrow length normalization:", arrow_length_normalization, "\n")
cat("Arrow length multiplier range: 0.5 to 2 times the previous default plotted length\n")
cat("Arrow linewidth:", arrow_linewidth, "\n")
cat("Arrow head cm:", arrow_head_cm, "\n")
cat("Arrow head angle:", arrow_head_angle, "\n")
cat("Arrow head angle by R2:", arrow_head_angle_by_r2, "\n")
cat("Arrow head angle range:", arrow_head_angle_min, "to", arrow_head_angle_max, "\n")
cat("Arrow head angle R2 anchors:", arrow_head_angle_r2_min, arrow_head_angle_r2_mid, arrow_head_angle_r2_max, "\n")
cat("Arrow head angle step:", arrow_head_angle_step, "\n")
cat("Arrow closed:", arrow_closed, "\n")
cat("Reused LSGI result:", did_reuse_lsgi, "\n\n")
cat("Outputs:\n")
cat("- intermediate/11_lsgi_cell_component_result.rds.gz\n")
cat("- output/11_lsgi_gradient/grid_info.csv\n")
cat("- output/11_lsgi_gradient/spot_grid_membership.csv\n")
cat("- output/11_lsgi_gradient/grid_spot_summary.csv\n")
cat("- output/11_lsgi_gradient/grid_local_spot_membership.csv\n")
cat("- output/11_lsgi_gradient/grid_partition_membership.csv\n")
cat("- output/11_lsgi_gradient/grid_partition.pdf\n")
cat("- output/11_lsgi_gradient/grid_partition_overlay.pdf\n")
cat("- output/11_lsgi_gradient/grid_local_spots.pdf\n")
cat("- output/11_lsgi_gradient/grid_local_spots_overlay.pdf\n")
cat("- output/11_lsgi_gradient/arrow_tables/cell_component_arrows_by_grid.csv\n")
sink()

message("Done. Detailed LSGI grid outputs written to ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), ".")
