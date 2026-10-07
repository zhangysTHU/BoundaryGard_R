# 06: Tumor-arrow-guided boundary profile.
# This downstream analysis uses boundary-related tumor/cancer epithelial arrows
# from 05 as local axes, aligns each axis to one nearby Bdy spot when possible,
# and summarizes selected embedding features along the local axis.
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "R", "boundarygrad_core.R"))
source(file.path(script_dir, "R", "spatial_plot_core.R"))
load_required_packages(c("readr", "dplyr", "tibble", "ggplot2", "jsonlite", "png", "grid", "patchwork"))

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

parse_chr_vec <- function(x, default = character()) {
  if (is.null(x)) return(default)
  values <- if (length(x) > 1) {
    as.character(x)
  } else {
    strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  }
  values <- trimws(values)
  values[nzchar(values)]
}

sanitize_filename <- function(x, prefix = "item") {
  out <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  empty <- !nzchar(out) | is.na(out)
  out[empty] <- paste0(prefix, seq_len(sum(empty)))
  make.unique(out, sep = "_")
}

sanitize_component_names <- function(x, prefix = "component") {
  out <- make.names(as.character(x), unique = TRUE)
  empty <- !nzchar(out) | is.na(out)
  out[empty] <- paste0(prefix, seq_len(sum(empty)))
  out
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path, call. = FALSE)
  }
  path
}

read_csv_df <- function(path, ...) {
  as.data.frame(readr::read_csv(path, show_col_types = FALSE, ...), stringsAsFactors = FALSE)
}

read_profile_embedding_matrix <- function(method_id, lsgi_dir, output_dir) {
  if (!identical(method_id, "cell_component")) {
    embedding_path <- file.path(lsgi_dir, method_id, "embedding_matrix.csv")
    if (!file.exists(embedding_path)) return(NULL)
    embedding_df <- read_csv_df(embedding_path)
    if (!"cell_ID" %in% colnames(embedding_df)) {
      stop("Embedding matrix is missing cell_ID column: ", embedding_path, call. = FALSE)
    }
    embedding_df$cell_ID <- as.character(embedding_df$cell_ID)
    return(list(path = embedding_path, data = embedding_df, source = "11_lsgi_gradient"))
  }

  embedding_path <- file.path(lsgi_dir, "cell_component", "embedding_matrix.csv")
  if (file.exists(embedding_path)) {
    embedding_df <- read_csv_df(embedding_path)
    if (!"cell_ID" %in% colnames(embedding_df)) {
      stop("Embedding matrix is missing cell_ID column: ", embedding_path, call. = FALSE)
    }
    embedding_df$cell_ID <- as.character(embedding_df$cell_ID)
    return(list(path = embedding_path, data = embedding_df, source = "11_lsgi_gradient"))
  }

  decon_path <- file.path(output_dir, "07_spatial_deconvolution", "spot_matrix_pre_lsgi.csv")
  if (!file.exists(decon_path)) return(NULL)
  decon_df <- read_csv_df(decon_path)
  if (!"cell_ID" %in% colnames(decon_df)) {
    stop("Deconvolution spot matrix is missing cell_ID column: ", decon_path, call. = FALSE)
  }
  component_cols <- setdiff(colnames(decon_df), c("cell_ID", "X", "Y", "Location"))
  if (length(component_cols) == 0) {
    stop("Deconvolution spot matrix has no cell component columns: ", decon_path, call. = FALSE)
  }
  embedding_df <- decon_df[, c("cell_ID", component_cols), drop = FALSE]
  colnames(embedding_df) <- c("cell_ID", sanitize_component_names(component_cols, prefix = "cell_component_"))
  embedding_df$cell_ID <- as.character(embedding_df$cell_ID)
  list(path = decon_path, data = embedding_df, source = "07_spatial_deconvolution")
}

empty_df <- function(cols) {
  out <- as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE)
  out
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  if (all(abs(x) < .Machine$double.eps^0.5)) return(1)
  out <- tryCatch(
    suppressWarnings(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  as.numeric(out)
}

estimate_spot_dimensions <- function(spatial_dir) {
  scale_path <- require_file(file.path(spatial_dir, "scalefactors_json.json"), "Space Ranger scalefactors JSON")
  scale_factors <- jsonlite::fromJSON(scale_path)
  if (is.null(scale_factors$spot_diameter_fullres) || is.null(scale_factors$tissue_lowres_scalef)) {
    stop("scalefactors_json.json is missing spot_diameter_fullres or tissue_lowres_scalef.", call. = FALSE)
  }
  lowres_scalef <- as.numeric(scale_factors$tissue_lowres_scalef)
  spot_diameter_px <- as.numeric(scale_factors$spot_diameter_fullres) * lowres_scalef

  list_path <- file.path(spatial_dir, "tissue_positions_list.csv")
  header_path <- file.path(spatial_dir, "tissue_positions.csv")
  if (file.exists(list_path)) {
    tissue_pos <- readr::read_csv(
      list_path,
      col_names = c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"),
      show_col_types = FALSE
    )
  } else if (file.exists(header_path)) {
    tissue_pos <- readr::read_csv(header_path, show_col_types = FALSE)
  } else {
    stop("Cannot find tissue_positions_list.csv or tissue_positions.csv under: ", spatial_dir, call. = FALSE)
  }
  tissue_pos <- as.data.frame(tissue_pos, stringsAsFactors = FALSE)
  required_cols <- c("in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  missing_cols <- setdiff(required_cols, colnames(tissue_pos))
  if (length(missing_cols) > 0) {
    stop("Tissue positions file is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  tissue_pos <- tissue_pos[as.numeric(tissue_pos$in_tissue) == 1, , drop = FALSE]
  if (nrow(tissue_pos) < 2) {
    stop("Too few in-tissue spots to estimate spot pitch.", call. = FALSE)
  }

  tissue_pos$array_row <- as.numeric(tissue_pos$array_row)
  tissue_pos$array_col <- as.numeric(tissue_pos$array_col)
  tissue_pos$pxl_row_in_fullres <- as.numeric(tissue_pos$pxl_row_in_fullres)
  tissue_pos$pxl_col_in_fullres <- as.numeric(tissue_pos$pxl_col_in_fullres)

  keys <- paste(tissue_pos$array_row, tissue_pos$array_col, sep = ",")
  row_index <- seq_len(nrow(tissue_pos))
  names(row_index) <- keys
  offsets <- matrix(c(0, 2, 1, 1, 1, -1), ncol = 2, byrow = TRUE)
  adjacent_dist <- numeric()
  for (i in seq_len(nrow(offsets))) {
    target_keys <- paste(tissue_pos$array_row + offsets[i, 1], tissue_pos$array_col + offsets[i, 2], sep = ",")
    ok <- target_keys %in% names(row_index)
    if (!any(ok)) next
    src_idx <- which(ok)
    dst_idx <- unname(row_index[target_keys[ok]])
    dist_fullres <- sqrt(
      (tissue_pos$pxl_col_in_fullres[src_idx] - tissue_pos$pxl_col_in_fullres[dst_idx])^2 +
        (tissue_pos$pxl_row_in_fullres[src_idx] - tissue_pos$pxl_row_in_fullres[dst_idx])^2
    )
    adjacent_dist <- c(adjacent_dist, dist_fullres * lowres_scalef)
  }

  if (length(adjacent_dist) == 0) {
    xy <- as.matrix(tissue_pos[, c("pxl_col_in_fullres", "pxl_row_in_fullres"), drop = FALSE]) * lowres_scalef
    if (nrow(xy) > 1000) {
      set.seed(666)
      xy <- xy[sample(seq_len(nrow(xy)), 1000), , drop = FALSE]
    }
    dist_mat <- as.matrix(stats::dist(xy))
    diag(dist_mat) <- NA_real_
    adjacent_dist <- apply(dist_mat, 1, min, na.rm = TRUE)
  }

  adjacent_dist <- adjacent_dist[is.finite(adjacent_dist) & adjacent_dist > 0]
  if (length(adjacent_dist) == 0) {
    stop("Failed to estimate spot pitch from tissue positions.", call. = FALSE)
  }

  list(
    scalefactors_path = scale_path,
    tissue_positions_path = if (file.exists(list_path)) list_path else header_path,
    tissue_lowres_scalef = lowres_scalef,
    spot_diameter_px = spot_diameter_px,
    spot_radius_px = spot_diameter_px / 2,
    spot_pitch_px = stats::median(adjacent_dist, na.rm = TRUE),
    n_adjacent_pairs = length(adjacent_dist)
  )
}

resolve_px_mode <- function(mode, dims, label) {
  mode_chr <- trimws(as.character(mode))
  numeric_value <- suppressWarnings(as.numeric(mode_chr))
  if (length(mode_chr) == 1 && is.finite(numeric_value) && numeric_value > 0) {
    return(numeric_value)
  }
  key <- gsub("-", "_", tolower(mode_chr))
  out <- switch(
    key,
    spot_radius = dims$spot_radius_px,
    spot_diameter = dims$spot_diameter_px,
    half_pitch = dims$spot_pitch_px / 2,
    one_pitch = dims$spot_pitch_px,
    spot_pitch = dims$spot_pitch_px,
    two_pitch = dims$spot_pitch_px * 2,
    NA_real_
  )
  if (!is.finite(out) || out <= 0) {
    stop("Unsupported ", label, " mode: ", mode_chr, call. = FALSE)
  }
  out
}

assign_axis_clusters <- function(x, gap_px) {
  if (length(x) == 0) return(integer())
  ord <- order(x)
  clusters_sorted <- cumsum(c(TRUE, diff(x[ord]) > gap_px))
  out <- integer(length(x))
  out[ord] <- clusters_sorted
  out
}

summarise_profile_set <- function(patch_bin_df, patches_df, set_id, patch_filter) {
  patch_ids <- patches_df$patch_id[patch_filter(patches_df)]
  if (length(patch_ids) == 0) {
    return(data.frame())
  }
  df <- patch_bin_df[patch_bin_df$patch_id %in% patch_ids, , drop = FALSE]
  if (nrow(df) == 0) {
    return(data.frame())
  }
  df %>%
    dplyr::group_by(component_method, feature_key, feature, distance_bin_index, bin_start_px, bin_end_px, bin_center_px) %>%
    dplyr::summarise(
      analysis_set = set_id,
      n_patches = dplyr::n_distinct(patch_id),
      n_patch_bins = dplyr::n(),
      total_spots = sum(n_spots, na.rm = TRUE),
      mean_score = mean(bin_mean, na.rm = TRUE),
      median_score = stats::median(bin_mean, na.rm = TRUE),
      sd_score = safe_sd(bin_mean),
      se_score = safe_se(bin_mean),
      ci_low = ifelse(is.finite(se_score), mean_score - 1.96 * se_score, NA_real_),
      ci_high = ifelse(is.finite(se_score), mean_score + 1.96 * se_score, NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::select(analysis_set, dplyr::everything()) %>%
    dplyr::arrange(component_method, feature, bin_center_px)
}

make_band_tests <- function(patch_bin_df, bin_width_px) {
  if (nrow(patch_bin_df) == 0) {
    return(data.frame(
      component_method = character(),
      feature_key = character(),
      feature = character(),
      test_id = character(),
      n_patches = integer(),
      estimate_mean = numeric(),
      estimate_median = numeric(),
      p_value = numeric(),
      p_adj = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  combos <- unique(patch_bin_df[, c("component_method", "feature_key", "feature"), drop = FALSE])
  rows <- list()
  row_i <- 0L
  for (i in seq_len(nrow(combos))) {
    method_id <- combos$component_method[i]
    feature_key <- combos$feature_key[i]
    feature_id <- combos$feature[i]
    df <- patch_bin_df[
      patch_bin_df$component_method == method_id & patch_bin_df$feature_key == feature_key & patch_bin_df$feature == feature_id,
      ,
      drop = FALSE
    ]
    patch_ids <- unique(df$patch_id)

    boundary_minus_tumor <- numeric()
    boundary_minus_stroma <- numeric()
    tumor_slopes <- numeric()
    stroma_slopes <- numeric()

    for (patch_id in patch_ids) {
      patch_df <- df[df$patch_id == patch_id, , drop = FALSE]
      boundary_vals <- patch_df$bin_mean[abs(patch_df$bin_center_px) <= bin_width_px / 2]
      tumor_vals <- patch_df$bin_mean[patch_df$bin_center_px > bin_width_px / 2]
      stroma_vals <- patch_df$bin_mean[patch_df$bin_center_px < -bin_width_px / 2]
      if (length(boundary_vals) > 0 && length(tumor_vals) > 0) {
        boundary_minus_tumor <- c(boundary_minus_tumor, mean(boundary_vals, na.rm = TRUE) - mean(tumor_vals, na.rm = TRUE))
      }
      if (length(boundary_vals) > 0 && length(stroma_vals) > 0) {
        boundary_minus_stroma <- c(boundary_minus_stroma, mean(boundary_vals, na.rm = TRUE) - mean(stroma_vals, na.rm = TRUE))
      }

      tumor_df <- patch_df[patch_df$bin_center_px > bin_width_px / 2, , drop = FALSE]
      if (nrow(tumor_df) >= 2 && length(unique(tumor_df$bin_center_px)) >= 2) {
        tumor_slopes <- c(tumor_slopes, stats::coef(stats::lm(bin_mean ~ bin_center_px, data = tumor_df))[["bin_center_px"]])
      }
      stroma_df <- patch_df[patch_df$bin_center_px < -bin_width_px / 2, , drop = FALSE]
      if (nrow(stroma_df) >= 2 && length(unique(stroma_df$bin_center_px)) >= 2) {
        stroma_slopes <- c(stroma_slopes, stats::coef(stats::lm(bin_mean ~ bin_center_px, data = stroma_df))[["bin_center_px"]])
      }
    }

    tests <- list(
      boundary_near_vs_tumor_side = boundary_minus_tumor,
      boundary_near_vs_stroma_side = boundary_minus_stroma,
      tumor_side_slope = tumor_slopes,
      stroma_side_slope = stroma_slopes
    )
    for (test_id in names(tests)) {
      values <- tests[[test_id]]
      values <- values[is.finite(values)]
      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        component_method = method_id,
        feature_key = feature_key,
        feature = feature_id,
        test_id = test_id,
        n_patches = length(values),
        estimate_mean = ifelse(length(values) > 0, mean(values), NA_real_),
        estimate_median = ifelse(length(values) > 0, stats::median(values), NA_real_),
        p_value = safe_wilcox_p(values),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- dplyr::bind_rows(rows)
  out$p_adj <- NA_real_
  finite_p <- is.finite(out$p_value)
  out$p_adj[finite_p] <- stats::p.adjust(out$p_value[finite_p], method = "BH")
  out
}

plot_distance_profile <- function(summary_df, out_path, title_text) {
  if (nrow(summary_df) == 0) return(invisible(FALSE))
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = bin_center_px, y = mean_score)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey40") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      fill = "#9ecae1",
      alpha = 0.35,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(color = "#08519c", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(size = n_patches), color = "#08519c", alpha = 0.8, na.rm = TRUE) +
    ggplot2::scale_size_continuous(range = c(1.2, 3.5)) +
    ggplot2::labs(
      title = title_text,
      x = "Aligned distance along tumor arrow axis (lowres px)",
      y = "Mean embedding score",
      size = "Patches"
    ) +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(out_path, p, width = 6.5, height = 4.5)
  invisible(TRUE)
}

profile_class_order <- c("cross_boundary", "one_sided", "no_zero", "no_direction", "no_local_spots")
profile_class_labels <- c(
  cross_boundary = "cross-boundary\nzero found",
  one_sided = "zero found\none side only",
  no_zero = "no zero\nfound",
  no_direction = "no direction",
  no_local_spots = "no local spots"
)

profile_class_label <- function(x) {
  x_chr <- as.character(x)
  out <- unname(profile_class_labels[x_chr])
  missing <- is.na(out)
  out[missing] <- x_chr[missing]
  out
}

profile_class_label_factor <- function(x) {
  factor(profile_class_label(x), levels = unname(profile_class_labels[profile_class_order]))
}

make_patch_rectangles <- function(patches_df) {
  empty <- data.frame(
    patch_id = character(),
    sample = character(),
    strategy_id = character(),
    grid = character(),
    tumor_component = character(),
    profile_class = character(),
    passes_min_spots_filter = logical(),
    passes_normal_filter = logical(),
    profile_exclusion_reason = character(),
    primary_in_analysis = logical(),
    corner_order = integer(),
    X = numeric(),
    Y = numeric(),
    stringsAsFactors = FALSE
  )
  if (nrow(patches_df) == 0) return(empty)
  keep <- is.finite(patches_df$arrow_unit_vx) &
    is.finite(patches_df$arrow_unit_vy) &
    is.finite(patches_df$zero_axis_X) &
    is.finite(patches_df$zero_axis_Y) &
    is.finite(patches_df$negative_limit_px) &
    is.finite(patches_df$positive_limit_px) &
    is.finite(patches_df$tube_half_width_px) &
    patches_df$negative_limit_px < patches_df$positive_limit_px &
    patches_df$tube_half_width_px > 0
  if (!any(keep)) return(empty)

  rows <- vector("list", sum(keep))
  row_i <- 0L
  for (idx in which(keep)) {
    patch <- patches_df[idx, , drop = FALSE]
    ux <- as.numeric(patch$arrow_unit_vx[[1]])
    uy <- as.numeric(patch$arrow_unit_vy[[1]])
    tx <- -uy
    ty <- ux
    w <- as.numeric(patch$tube_half_width_px[[1]])
    x0 <- as.numeric(patch$zero_axis_X[[1]])
    y0 <- as.numeric(patch$zero_axis_Y[[1]])
    neg <- as.numeric(patch$negative_limit_px[[1]])
    pos <- as.numeric(patch$positive_limit_px[[1]])

    start_x <- x0 + neg * ux
    start_y <- y0 + neg * uy
    end_x <- x0 + pos * ux
    end_y <- y0 + pos * uy
    x <- c(start_x + w * tx, end_x + w * tx, end_x - w * tx, start_x - w * tx, start_x + w * tx)
    y <- c(start_y + w * ty, end_y + w * ty, end_y - w * ty, start_y - w * ty, start_y + w * ty)

    row_i <- row_i + 1L
    rows[[row_i]] <- data.frame(
      patch_id = rep(patch$patch_id[[1]], 5),
      sample = rep(patch$sample[[1]], 5),
      strategy_id = rep(patch$strategy_id[[1]], 5),
      grid = rep(patch$grid[[1]], 5),
      tumor_component = rep(patch$tumor_component[[1]], 5),
      profile_class = rep(patch$profile_class[[1]], 5),
      passes_min_spots_filter = rep(patch$passes_min_spots_filter[[1]], 5),
      passes_normal_filter = rep(patch$passes_normal_filter[[1]], 5),
      profile_exclusion_reason = rep(patch$profile_exclusion_reason[[1]], 5),
      primary_in_analysis = rep(patch$primary_in_analysis[[1]], 5),
      corner_order = seq_len(5),
      X = x,
      Y = y,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

make_profile_plot_arrow_df <- function(df, component_levels, arrow_length_scale) {
  empty <- data.frame()
  if (nrow(df) == 0) return(empty)
  if (!"arrow_head_angle_mapped" %in% colnames(df)) {
    df$arrow_head_angle_mapped <- arrow_head_angle
  }
  out <- sp_prepare_arrows(
    df,
    arrow_length_scale = arrow_length_scale,
    x_col = "arrow_X",
    y_col = "arrow_Y",
    vx_col = "arrow_vx_u",
    vy_col = "arrow_vy_u",
    flip_y = TRUE
  )
  out$component <- factor(as.character(out$tumor_component), levels = component_levels)
  out$arrow_head_angle_mapped <- ifelse(
    is.finite(out$arrow_head_angle_mapped),
    out$arrow_head_angle_mapped,
    arrow_head_angle
  )
  out
}

add_profile_arrows <- function(p, arrow_df, component_palette = NULL, legend_title = "Tumor arrow", fixed_color = NULL, color_drop = FALSE) {
  out <- p
  if (nrow(arrow_df) == 0) return(out)
  angles <- sort(unique(arrow_df$arrow_head_angle_mapped[is.finite(arrow_df$arrow_head_angle_mapped)]))
  if (length(angles) == 0) angles <- arrow_head_angle
  for (angle_value in angles) {
    layer_data <- arrow_df[is.finite(arrow_df$arrow_head_angle_mapped) & arrow_df$arrow_head_angle_mapped == angle_value, , drop = FALSE]
    if (nrow(layer_data) == 0) next
    if (is.null(fixed_color)) {
      out <- out +
        ggplot2::geom_segment(
          data = layer_data,
          ggplot2::aes(x = X, y = Y, xend = X_end, yend = Y_end, color = component),
          inherit.aes = FALSE,
          linewidth = arrow_linewidth,
          lineend = "round",
          na.rm = TRUE,
          arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), angle = angle_value, type = arrow_type)
        )
    } else {
      out <- out +
        ggplot2::geom_segment(
          data = layer_data,
          ggplot2::aes(x = X, y = Y, xend = X_end, yend = Y_end),
          inherit.aes = FALSE,
          color = fixed_color,
          alpha = 0.75,
          linewidth = arrow_linewidth,
          lineend = "round",
          na.rm = TRUE,
          arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), angle = angle_value, type = arrow_type)
        )
    }
  }
  if (is.null(fixed_color)) {
    out <- out +
      ggplot2::scale_color_manual(values = component_palette, drop = color_drop) +
      ggplot2::labs(color = legend_title)
  }
  out
}

make_boundary_spot_base <- function(point_df, title_text) {
  point_df$Location <- factor(as.character(point_df$Location), levels = names(boundary_cols))
  spot_polygons <- sp_make_spot_polygons(
    point_df,
    plot_context,
    segments = params$spatial_plot_spot_segments %||% 24L
  )
  sp_boundary_base(
    spot_polygons,
    plot_context,
    boundary_cols,
    title = title_text,
    show_image = plot_context$has_image
  )
}

cli_opts <- parse_cli_options(commandArgs(trailingOnly = TRUE))
plot_mode <- bg_validate_plot_mode(get_opt(cli_opts, "plot-mode", default = params$plot_mode %||% "full"))

lsgi_dir <- file.path(paths$output, "11_lsgi_gradient")
boundary_arrow_dir <- file.path(paths$output, "12_boundary_related_lsgi_arrows")
strategy_id <- as.character(get_opt(cli_opts, "strategy", default = params$tumor_profile_strategy %||% "local_broad"))
tumor_components <- parse_chr_vec(
  get_opt(cli_opts, "tumor-components", default = params$tumor_profile_tumor_components %||% c("Cancer.Epithelial", "Tumor")),
  default = c("Cancer.Epithelial", "Tumor")
)
profile_methods <- parse_chr_vec(
  get_opt(cli_opts, "methods", default = params$tumor_profile_methods %||% c("cell_component", "marker_module", "single_gene")),
  default = c("cell_component", "marker_module", "single_gene")
)
profile_features <- parse_chr_vec(
  get_opt(cli_opts, "features", default = params$tumor_profile_features %||% character()),
  default = character()
)
positive_control_features <- parse_chr_vec(
  get_opt(cli_opts, "positive-control-features", default = params$tumor_profile_positive_control_features %||% c("Cancer.Epithelial")),
  default = c("Cancer.Epithelial")
)
tube_half_width_mode <- as.character(get_opt(cli_opts, "tube-half-width-mode", default = params$tumor_profile_tube_half_width_mode %||% "half_pitch"))
zero_width_modes <- parse_chr_vec(
  get_opt(cli_opts, "zero-width-modes", default = params$tumor_profile_zero_width_modes %||% c("spot_radius", "half_pitch")),
  default = c("spot_radius", "half_pitch")
)
bin_width_mode <- as.character(get_opt(cli_opts, "bin-width-mode", default = params$tumor_profile_bin_width_mode %||% "half_pitch"))
bdy_cluster_gap_mode <- as.character(get_opt(cli_opts, "bdy-cluster-gap-mode", default = params$tumor_profile_bdy_cluster_gap_mode %||% "one_pitch"))
min_profile_spots <- as.integer(as.numeric(get_opt(cli_opts, "min-profile-spots", default = params$tumor_profile_min_profile_spots %||% 3)))
require_normal_spot <- as_bool(get_opt(cli_opts, "require-normal-spot", default = params$tumor_profile_require_normal_spot %||% TRUE))
normal_location <- as.character(get_opt(cli_opts, "normal-location", default = params$tumor_profile_normal_location %||% "nMal"))
plot_space <- identical(plot_mode, "full") &&
  as_bool(get_opt(cli_opts, "plot-space", default = TRUE))
clean_output <- as_bool(get_opt(cli_opts, "clean-output", default = TRUE))
out_dir <- as.character(get_opt(
  cli_opts,
  "out-dir",
  default = file.path(paths$output, "13_tumor_arrow_guided_boundary_profile")
))
arrow_linewidth <- as.numeric(params$lsgi_arrow_linewidth %||% 1.0)
arrow_head_cm <- as.numeric(params$lsgi_arrow_head_cm %||% 0.20)
arrow_length_scale <- as.numeric(get_opt(
  cli_opts,
  "arrow-length-scale",
  default = params$lsgi_arrow_length_scale %||% 1.4
))
arrow_head_angle <- as.numeric(params$lsgi_arrow_head_angle %||% 30)
arrow_closed <- as_bool(params$lsgi_arrow_closed %||% TRUE)
arrow_type <- if (arrow_closed) "closed" else "open"
boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")

if (length(profile_methods) == 1 && identical(tolower(profile_methods), "all")) {
  profile_methods <- c("cell_component", "marker_module", "single_gene", "pathway")
}
allowed_profile_methods <- c("cell_component", "marker_module", "single_gene", "pathway")
unknown_methods <- setdiff(profile_methods, allowed_profile_methods)
if (length(unknown_methods) > 0) {
  stop("Unsupported profile method(s): ", paste(unknown_methods, collapse = ", "), call. = FALSE)
}
if (is.na(min_profile_spots) || min_profile_spots < 1) {
  stop("--min-profile-spots must be an integer >= 1.", call. = FALSE)
}
if (!nzchar(normal_location)) {
  stop("--normal-location must be a non-empty Location label.", call. = FALSE)
}

bg_assert_safe_child(out_dir, paths$output, "profile output directory")
if (dir.exists(out_dir) && clean_output) {
  unlink(out_dir, recursive = TRUE, force = TRUE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots")
profile_plot_dir <- file.path(plot_dir, "distance_profiles")
if (!identical(plot_mode, "none")) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(profile_plot_dir, recursive = TRUE, showWarnings = FALSE)
}

spatial_dir <- file.path(paths$spaceranger, "spatial")
dims <- estimate_spot_dimensions(spatial_dir)
plot_context <- sp_read_context(paths$spaceranger, load_image = plot_space)
if (!isTRUE(all.equal(dims$spot_diameter_px, plot_context$spot_diameter_lowres, tolerance = 1e-10))) {
  stop("Inconsistent low-resolution spot diameter between profile and plotting geometry.", call. = FALSE)
}
tube_half_width_px <- resolve_px_mode(tube_half_width_mode, dims, "tube half-width")
zero_width_px <- vapply(zero_width_modes, resolve_px_mode, numeric(1), dims = dims, label = "zero-width")
bin_width_px <- resolve_px_mode(bin_width_mode, dims, "bin-width")
bdy_cluster_gap_px <- resolve_px_mode(bdy_cluster_gap_mode, dims, "Bdy cluster gap")

arrow_path <- require_file(
  file.path(boundary_arrow_dir, "cell_component", "arrow_tables", paste0(strategy_id, "_arrows.csv")),
  paste0("boundary-related cell_component arrow table for strategy ", strategy_id)
)
local_membership_path <- require_file(
  file.path(lsgi_dir, "grid_local_spot_membership.csv"),
  "LSGI grid local spot membership table"
)
spot_membership_path <- file.path(lsgi_dir, "spot_grid_membership.csv")

arrows <- read_csv_df(arrow_path)
local_df <- read_csv_df(local_membership_path)
required_arrow_cols <- c("grid", "component", "X", "Y", "vx", "vy", "vx.u", "vy.u")
missing_arrow_cols <- setdiff(required_arrow_cols, colnames(arrows))
if (length(missing_arrow_cols) > 0) {
  stop("Arrow table is missing column(s): ", paste(missing_arrow_cols, collapse = ", "), call. = FALSE)
}
required_local_cols <- c("grid", "cell_ID", "local_rank", "dist_to_grid", "X", "Y", "Location")
missing_local_cols <- setdiff(required_local_cols, colnames(local_df))
if (length(missing_local_cols) > 0) {
  stop("Local membership table is missing column(s): ", paste(missing_local_cols, collapse = ", "), call. = FALSE)
}

arrows$grid <- as.character(arrows$grid)
arrows$component <- as.character(arrows$component)
local_df$grid <- as.character(local_df$grid)
local_df$cell_ID <- as.character(local_df$cell_ID)
local_df$Location <- as.character(local_df$Location)

tumor_arrows <- arrows[arrows$component %in% tumor_components, , drop = FALSE]
if (nrow(tumor_arrows) == 0) {
  stop(
    "No tumor arrows found for component(s): ", paste(tumor_components, collapse = ", "),
    ". Available components: ", paste(sort(unique(arrows$component)), collapse = ", "),
    call. = FALSE
  )
}
tumor_arrows$arrow_row_id <- seq_len(nrow(tumor_arrows))
tumor_arrows$patch_id <- paste(
  strategy_id,
  tumor_arrows$grid,
  sanitize_filename(tumor_arrows$component),
  tumor_arrows$arrow_row_id,
  sep = "__"
)

patch_rows <- vector("list", nrow(tumor_arrows))
spot_rows <- vector("list", nrow(tumor_arrows))

for (i in seq_len(nrow(tumor_arrows))) {
  arrow <- tumor_arrows[i, , drop = FALSE]
  patch_id <- arrow$patch_id[[1]]
  grid_id <- arrow$grid[[1]]
  local_spots <- local_df[local_df$grid == grid_id, , drop = FALSE]
  arrow_x <- as.numeric(arrow$X[[1]])
  arrow_y <- as.numeric(arrow$Y[[1]])

  if ("vx.u" %in% colnames(arrow) && "vy.u" %in% colnames(arrow)) {
    direction <- c(as.numeric(arrow[["vx.u"]][[1]]), as.numeric(arrow[["vy.u"]][[1]]))
    direction_source <- "vx.u/vy.u"
  } else {
    direction <- c(NA_real_, NA_real_)
    direction_source <- "missing"
  }
  direction_norm <- sqrt(sum(direction^2))
  if (!is.finite(direction_norm) || direction_norm <= 0) {
    direction <- c(as.numeric(arrow$vx[[1]]), as.numeric(arrow$vy[[1]]))
    direction_source <- "vx/vy"
    direction_norm <- sqrt(sum(direction^2))
  }

  base_patch <- data.frame(
    patch_id = patch_id,
    sample = sample_name,
    strategy_id = strategy_id,
    grid = grid_id,
    tumor_component = arrow$component[[1]],
    arrow_X = arrow_x,
    arrow_Y = arrow_y,
    arrow_vx = as.numeric(arrow$vx[[1]]),
    arrow_vy = as.numeric(arrow$vy[[1]]),
    arrow_vx_u = as.numeric(arrow[["vx.u"]][[1]]),
    arrow_vy_u = as.numeric(arrow[["vy.u"]][[1]]),
    arrow_direction_source = direction_source,
    arrow_unit_vx = NA_real_,
    arrow_unit_vy = NA_real_,
    rsquared = if ("rsquared" %in% colnames(arrow)) as.numeric(arrow$rsquared[[1]]) else NA_real_,
    gradient_strength = if ("gradient_strength" %in% colnames(arrow)) as.numeric(arrow$gradient_strength[[1]]) else NA_real_,
    gradient_strength_norm = if ("gradient_strength_norm" %in% colnames(arrow)) as.numeric(arrow$gradient_strength_norm[[1]]) else NA_real_,
    plotted_length = if ("plotted_length" %in% colnames(arrow)) as.numeric(arrow$plotted_length[[1]]) else NA_real_,
    arrow_head_angle_mapped = if ("arrow_head_angle_mapped" %in% colnames(arrow)) as.numeric(arrow$arrow_head_angle_mapped[[1]]) else NA_real_,
    local_n_Mal_arrow = if ("local_n_Mal" %in% colnames(arrow)) as.integer(arrow$local_n_Mal[[1]]) else NA_integer_,
    local_n_Bdy_arrow = if ("local_n_Bdy" %in% colnames(arrow)) as.integer(arrow$local_n_Bdy[[1]]) else NA_integer_,
    local_n_nMal_arrow = if ("local_n_nMal" %in% colnames(arrow)) as.integer(arrow$local_n_nMal[[1]]) else NA_integer_,
    local_frac_Bdy_arrow = if ("local_frac_Bdy" %in% colnames(arrow)) as.numeric(arrow$local_frac_Bdy[[1]]) else NA_real_,
    spot_pitch_px = dims$spot_pitch_px,
    spot_radius_px = dims$spot_radius_px,
    tube_half_width_px = tube_half_width_px,
    bin_width_px = bin_width_px,
    local_radius_px = ifelse(nrow(local_spots) > 0, max(as.numeric(local_spots$dist_to_grid), na.rm = TRUE), NA_real_),
    zero_found = FALSE,
    zero_width_mode = NA_character_,
    zero_width_px = NA_real_,
    zero_cell_ID = NA_character_,
    zero_spot_X = NA_real_,
    zero_spot_Y = NA_real_,
    zero_axis_X = NA_real_,
    zero_axis_Y = NA_real_,
    zero_axis_s_from_arrow_px = NA_real_,
    zero_lateral_distance_px = NA_real_,
    negative_limit_px = NA_real_,
    positive_limit_px = NA_real_,
    n_bdy_axis_candidates = 0L,
    n_bdy_tube_clusters = 0L,
    truncated_negative_by_bdy_cluster = FALSE,
    truncated_positive_by_bdy_cluster = FALSE,
    n_profile_spots = 0L,
    n_negative_side_spots = 0L,
    n_positive_side_spots = 0L,
    n_boundary_near_spots = 0L,
    n_profile_Mal = 0L,
    n_profile_Bdy = 0L,
    n_profile_nMal = 0L,
    n_profile_normal_spots = 0L,
    median_profile_x_Mal = NA_real_,
    median_profile_x_Bdy = NA_real_,
    median_profile_x_nMal = NA_real_,
    profile_class = "no_direction",
    require_normal_spot = require_normal_spot,
    normal_location = normal_location,
    passes_min_spots_filter = FALSE,
    passes_normal_filter = !require_normal_spot,
    profile_exclusion_reason = "no_direction",
    primary_in_analysis = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(local_spots) == 0) {
    base_patch$profile_class <- "no_local_spots"
    base_patch$profile_exclusion_reason <- "no_local_spots"
    patch_rows[[i]] <- base_patch
    spot_rows[[i]] <- data.frame()
    next
  }
  if (!is.finite(direction_norm) || direction_norm <= 0) {
    patch_rows[[i]] <- base_patch
    spot_rows[[i]] <- data.frame()
    next
  }

  projection <- bg_project_axis(
    local_spots$X, local_spots$Y, arrow_x, arrow_y,
    direction[[1]], direction[[2]]
  )
  unit_v <- projection$unit_v
  unit_t <- projection$unit_t
  axis_s <- projection$axis_s
  lateral_signed <- projection$lateral_signed
  lateral_abs <- abs(lateral_signed)
  local_radius_px <- base_patch$local_radius_px[[1]]
  finite_extension <- is.finite(axis_s) & abs(axis_s) <= local_radius_px

  zero <- bg_select_zero_spot(
    local_spots$Location, lateral_abs, axis_s, finite_extension,
    local_spots$local_rank, zero_width_modes, zero_width_px
  )
  zero_found <- zero$found
  zero_idx <- zero$index
  zero_mode <- zero$mode
  zero_px <- zero$width
  zero_candidate_count <- zero$candidate_count
  centered_axis <- bg_center_axis_on_zero(axis_s, zero_idx)
  zero_s <- centered_axis$zero_s
  profile_x <- centered_axis$profile_x
  negative_limit <- -local_radius_px
  positive_limit <- local_radius_px
  bdy_tube_idx <- which(
    local_spots$Location == "Bdy" &
      lateral_abs <= tube_half_width_px &
      is.finite(profile_x) &
      abs(profile_x) <= local_radius_px
  )
  bdy_cluster_count <- 0L
  truncated_negative <- FALSE
  truncated_positive <- FALSE
  if (zero_found && length(bdy_tube_idx) > 0) {
    bdy_clusters <- assign_axis_clusters(profile_x[bdy_tube_idx], bdy_cluster_gap_px)
    bdy_cluster_count <- length(unique(bdy_clusters))
    zero_in_bdy <- match(zero_idx, bdy_tube_idx)
    if (!is.na(zero_in_bdy)) {
      zero_cluster <- bdy_clusters[[zero_in_bdy]]
      cluster_df <- data.frame(
        cluster = sort(unique(bdy_clusters)),
        min_x = vapply(sort(unique(bdy_clusters)), function(cl) min(profile_x[bdy_tube_idx][bdy_clusters == cl], na.rm = TRUE), numeric(1)),
        max_x = vapply(sort(unique(bdy_clusters)), function(cl) max(profile_x[bdy_tube_idx][bdy_clusters == cl], na.rm = TRUE), numeric(1)),
        stringsAsFactors = FALSE
      )
      zero_cluster_df <- cluster_df[cluster_df$cluster == zero_cluster, , drop = FALSE]
      other_clusters <- cluster_df[cluster_df$cluster != zero_cluster, , drop = FALSE]
      positive_clusters <- other_clusters[
        other_clusters$min_x > zero_cluster_df$max_x &
          other_clusters$min_x - zero_cluster_df$max_x >= bdy_cluster_gap_px,
        ,
        drop = FALSE
      ]
      negative_clusters <- other_clusters[
        other_clusters$max_x < zero_cluster_df$min_x &
          zero_cluster_df$min_x - other_clusters$max_x >= bdy_cluster_gap_px,
        ,
        drop = FALSE
      ]
      if (nrow(positive_clusters) > 0) {
        positive_limit <- min(positive_limit, min(positive_clusters$min_x, na.rm = TRUE) - dims$spot_radius_px)
        truncated_positive <- TRUE
      }
      if (nrow(negative_clusters) > 0) {
        negative_limit <- max(negative_limit, max(negative_clusters$max_x, na.rm = TRUE) + dims$spot_radius_px)
        truncated_negative <- TRUE
      }
      if (!is.finite(negative_limit) || !is.finite(positive_limit) || negative_limit >= positive_limit) {
        negative_limit <- -local_radius_px
        positive_limit <- local_radius_px
        truncated_negative <- FALSE
        truncated_positive <- FALSE
      }
    }
  }

  in_tube <- lateral_abs <= tube_half_width_px & is.finite(profile_x) & abs(profile_x) <= local_radius_px
  in_profile <- in_tube & profile_x >= negative_limit & profile_x <= positive_limit
  side_eps <- dims$spot_radius_px
  side <- ifelse(
    profile_x > side_eps,
    "tumor_direction",
    ifelse(profile_x < -side_eps, "opposite_direction", "boundary_near")
  )
  # Center bins at zero so the selected Bdy spot is not forced into the positive-side bin.
  distance_bin_index <- floor(profile_x / bin_width_px + 0.5)
  bin_center_px <- distance_bin_index * bin_width_px
  bin_start_px <- bin_center_px - bin_width_px / 2
  bin_end_px <- bin_center_px + bin_width_px / 2

  profile_locations <- local_spots$Location[in_profile]
  profile_x_kept <- profile_x[in_profile]
  n_profile_total <- sum(in_profile, na.rm = TRUE)
  n_neg <- sum(in_profile & profile_x < -side_eps, na.rm = TRUE)
  n_pos <- sum(in_profile & profile_x > side_eps, na.rm = TRUE)
  n_near <- sum(in_profile & abs(profile_x) <= side_eps, na.rm = TRUE)
  n_normal <- sum(profile_locations == normal_location, na.rm = TRUE)
  profile_class <- if (!zero_found) {
    "no_zero"
  } else if (n_neg > 0 && n_pos > 0) {
    "cross_boundary"
  } else {
    "one_sided"
  }
  passes_min_spots <- n_profile_total >= min_profile_spots
  passes_normal <- !require_normal_spot || n_normal > 0
  primary_in_analysis <- zero_found && profile_class == "cross_boundary" && passes_min_spots && passes_normal
  exclusion_reason <- character()
  if (!zero_found) {
    exclusion_reason <- c(exclusion_reason, "no_zero")
  }
  if (zero_found && profile_class != "cross_boundary") {
    exclusion_reason <- c(exclusion_reason, profile_class)
  }
  if (!passes_min_spots) {
    exclusion_reason <- c(exclusion_reason, "too_few_profile_spots")
  }
  if (!passes_normal) {
    exclusion_reason <- c(exclusion_reason, paste0("no_", normal_location, "_spot"))
  }
  if (primary_in_analysis) {
    exclusion_reason <- "primary"
  } else if (length(exclusion_reason) == 0) {
    exclusion_reason <- "excluded"
  } else {
    exclusion_reason <- paste(unique(exclusion_reason), collapse = ";")
  }

  base_patch$arrow_unit_vx <- unit_v[[1]]
  base_patch$arrow_unit_vy <- unit_v[[2]]
  base_patch$zero_found <- zero_found
  base_patch$zero_width_mode <- zero_mode
  base_patch$zero_width_px <- zero_px
  base_patch$zero_cell_ID <- if (zero_found) local_spots$cell_ID[[zero_idx]] else NA_character_
  base_patch$zero_spot_X <- if (zero_found) as.numeric(local_spots$X[[zero_idx]]) else NA_real_
  base_patch$zero_spot_Y <- if (zero_found) as.numeric(local_spots$Y[[zero_idx]]) else NA_real_
  base_patch$zero_axis_X <- arrow_x + unit_v[[1]] * zero_s
  base_patch$zero_axis_Y <- arrow_y + unit_v[[2]] * zero_s
  base_patch$zero_axis_s_from_arrow_px <- if (zero_found) zero_s else NA_real_
  base_patch$zero_lateral_distance_px <- if (zero_found) lateral_abs[[zero_idx]] else NA_real_
  base_patch$negative_limit_px <- negative_limit
  base_patch$positive_limit_px <- positive_limit
  base_patch$n_bdy_axis_candidates <- zero_candidate_count
  base_patch$n_bdy_tube_clusters <- bdy_cluster_count
  base_patch$truncated_negative_by_bdy_cluster <- truncated_negative
  base_patch$truncated_positive_by_bdy_cluster <- truncated_positive
  base_patch$n_profile_spots <- n_profile_total
  base_patch$n_negative_side_spots <- n_neg
  base_patch$n_positive_side_spots <- n_pos
  base_patch$n_boundary_near_spots <- n_near
  base_patch$n_profile_Mal <- sum(profile_locations == "Mal", na.rm = TRUE)
  base_patch$n_profile_Bdy <- sum(profile_locations == "Bdy", na.rm = TRUE)
  base_patch$n_profile_nMal <- sum(profile_locations == "nMal", na.rm = TRUE)
  base_patch$n_profile_normal_spots <- n_normal
  base_patch$median_profile_x_Mal <- ifelse(any(profile_locations == "Mal"), stats::median(profile_x_kept[profile_locations == "Mal"], na.rm = TRUE), NA_real_)
  base_patch$median_profile_x_Bdy <- ifelse(any(profile_locations == "Bdy"), stats::median(profile_x_kept[profile_locations == "Bdy"], na.rm = TRUE), NA_real_)
  base_patch$median_profile_x_nMal <- ifelse(any(profile_locations == "nMal"), stats::median(profile_x_kept[profile_locations == "nMal"], na.rm = TRUE), NA_real_)
  base_patch$profile_class <- profile_class
  base_patch$passes_min_spots_filter <- passes_min_spots
  base_patch$passes_normal_filter <- passes_normal
  base_patch$profile_exclusion_reason <- exclusion_reason
  base_patch$primary_in_analysis <- primary_in_analysis

  patch_rows[[i]] <- base_patch
  spot_rows[[i]] <- data.frame(
    patch_id = patch_id,
    sample = sample_name,
    strategy_id = strategy_id,
    grid = grid_id,
    tumor_component = arrow$component[[1]],
    cell_ID = local_spots$cell_ID,
    local_rank = as.integer(local_spots$local_rank),
    dist_to_grid_px = as.numeric(local_spots$dist_to_grid),
    X = as.numeric(local_spots$X),
    Y = as.numeric(local_spots$Y),
    Location = local_spots$Location,
    axis_s_from_arrow_px = axis_s,
    profile_x_px = profile_x,
    lateral_signed_px = lateral_signed,
    lateral_distance_px = lateral_abs,
    in_tube = in_tube,
    in_profile = in_profile,
    profile_side = side,
    distance_bin_index = distance_bin_index,
    bin_start_px = bin_start_px,
    bin_end_px = bin_end_px,
    bin_center_px = bin_center_px,
    zero_found = zero_found,
    zero_cell_ID = if (zero_found) local_spots$cell_ID[[zero_idx]] else NA_character_,
    profile_class = profile_class,
    passes_min_spots_filter = passes_min_spots,
    passes_normal_filter = passes_normal,
    profile_exclusion_reason = exclusion_reason,
    primary_in_analysis = base_patch$primary_in_analysis,
    stringsAsFactors = FALSE
  )
}

patches <- dplyr::bind_rows(patch_rows)
profile_spots <- dplyr::bind_rows(spot_rows)
patches$profile_class_label <- profile_class_label(patches$profile_class)
profile_spots$profile_class_label <- profile_class_label(profile_spots$profile_class)
patch_rectangles <- make_patch_rectangles(patches)
patch_rectangles$profile_class_label <- profile_class_label(patch_rectangles$profile_class)

readr::write_csv(patches, file.path(out_dir, "tumor_arrow_profile_patches.csv"))
readr::write_csv(profile_spots, file.path(out_dir, "tumor_arrow_profile_spots.csv"))
readr::write_csv(patch_rectangles, file.path(out_dir, "tumor_arrow_profile_rectangles.csv"))

feature_filter_all <- length(profile_features) == 1 && identical(tolower(profile_features), "all")
embedding_specs <- list()
feature_selection_rows <- list()
feature_selection_i <- 0L
found_features <- character()

for (method_id in profile_methods) {
  embedding_spec <- read_profile_embedding_matrix(method_id, lsgi_dir, paths$output)
  if (is.null(embedding_spec)) {
    warning("Skipping method ", method_id, " because no embedding source is available.", call. = FALSE)
    next
  }
  embedding_path <- embedding_spec$path
  embedding_df <- embedding_spec$data
  available_features <- setdiff(colnames(embedding_df), "cell_ID")
  selected_features <- if (feature_filter_all || length(profile_features) == 0) {
    available_features
  } else {
    intersect(profile_features, available_features)
  }
  if (identical(method_id, "cell_component") && !feature_filter_all && length(positive_control_features) > 0) {
    selected_features <- unique(c(intersect(positive_control_features, available_features), selected_features))
  }
  if (length(selected_features) == 0) {
    feature_selection_i <- feature_selection_i + 1L
    feature_selection_rows[[feature_selection_i]] <- data.frame(
      component_method = method_id,
      embedding_path = embedding_path,
      embedding_source = embedding_spec$source,
      feature_key = NA_character_,
      feature = NA_character_,
      selected = FALSE,
      reason = "no_requested_features_present",
      stringsAsFactors = FALSE
    )
    next
  }
  embedding_specs[[method_id]] <- list(
    path = embedding_path,
    source = embedding_spec$source,
    data = embedding_df,
    features = selected_features
  )
  found_features <- unique(c(found_features, selected_features))
  for (feature_id in selected_features) {
    feature_selection_i <- feature_selection_i + 1L
    feature_selection_rows[[feature_selection_i]] <- data.frame(
      component_method = method_id,
      embedding_path = embedding_path,
      embedding_source = embedding_spec$source,
      feature_key = paste(method_id, feature_id, sep = "::"),
      feature = feature_id,
      selected = TRUE,
      reason = ifelse(identical(method_id, "cell_component") && feature_id %in% positive_control_features, "positive_control", ""),
      stringsAsFactors = FALSE
    )
  }
}

feature_selection <- dplyr::bind_rows(feature_selection_rows)
if (nrow(feature_selection) == 0) {
  feature_selection <- data.frame(
    component_method = character(),
    embedding_path = character(),
    embedding_source = character(),
    feature_key = character(),
    feature = character(),
    selected = logical(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}
readr::write_csv(feature_selection, file.path(out_dir, "selected_features.csv"))

missing_requested <- if (feature_filter_all || length(profile_features) == 0) {
  character()
} else {
  setdiff(profile_features, found_features)
}
missing_requested_df <- data.frame(feature = missing_requested, stringsAsFactors = FALSE)
readr::write_csv(missing_requested_df, file.path(out_dir, "missing_requested_features.csv"))
selected_positive_controls <- feature_selection$feature[
  feature_selection$component_method == "cell_component" &
    feature_selection$selected &
    feature_selection$feature %in% positive_control_features
]
missing_positive_controls <- setdiff(positive_control_features, selected_positive_controls)
readr::write_csv(
  data.frame(feature = missing_positive_controls, stringsAsFactors = FALSE),
  file.path(out_dir, "missing_positive_control_features.csv")
)
if (length(embedding_specs) == 0) {
  stop("No embedding features are available for tumor-arrow-guided boundary profile.", call. = FALSE)
}

profile_spots_in <- profile_spots[profile_spots$in_profile, , drop = FALSE]
spot_value_rows <- list()
spot_value_i <- 0L
base_cols <- c(
  "patch_id", "sample", "strategy_id", "grid", "tumor_component", "cell_ID",
  "Location", "profile_x_px", "lateral_distance_px", "profile_side",
  "distance_bin_index", "bin_start_px", "bin_end_px", "bin_center_px",
  "zero_found", "profile_class", "profile_class_label", "passes_min_spots_filter",
  "passes_normal_filter", "profile_exclusion_reason", "primary_in_analysis"
)
profile_spots_base <- profile_spots_in[, base_cols, drop = FALSE]

for (method_id in names(embedding_specs)) {
  spec <- embedding_specs[[method_id]]
  embedding_df <- spec$data
  row_match <- match(profile_spots_base$cell_ID, embedding_df$cell_ID)
  for (feature_id in spec$features) {
    values <- as.numeric(embedding_df[[feature_id]][row_match])
    keep <- is.finite(values)
    if (!any(keep)) next
    spot_value_i <- spot_value_i + 1L
    spot_value_rows[[spot_value_i]] <- data.frame(
      profile_spots_base[keep, , drop = FALSE],
      component_method = method_id,
      feature_key = paste(method_id, feature_id, sep = "::"),
      feature = feature_id,
      embedding_value = values[keep],
      stringsAsFactors = FALSE
    )
  }
}

feature_spot_values <- dplyr::bind_rows(spot_value_rows)
if (nrow(feature_spot_values) == 0) {
  feature_spot_values <- data.frame(
    profile_spots_base[FALSE, , drop = FALSE],
    component_method = character(),
    feature_key = character(),
    feature = character(),
    embedding_value = numeric(),
    stringsAsFactors = FALSE
  )
}

patch_bin <- if (nrow(feature_spot_values) == 0) {
  data.frame()
} else {
  feature_spot_values %>%
    dplyr::group_by(
      component_method, feature, patch_id, strategy_id, grid, tumor_component,
      feature_key, profile_class, primary_in_analysis, zero_found,
      passes_min_spots_filter, passes_normal_filter, profile_exclusion_reason,
      distance_bin_index, bin_start_px, bin_end_px, bin_center_px
    ) %>%
    dplyr::summarise(
      n_spots = dplyr::n(),
      n_Mal = sum(Location == "Mal", na.rm = TRUE),
      n_Bdy = sum(Location == "Bdy", na.rm = TRUE),
      n_nMal = sum(Location == "nMal", na.rm = TRUE),
      bin_mean = mean(embedding_value, na.rm = TRUE),
      bin_median = stats::median(embedding_value, na.rm = TRUE),
      bin_sd = safe_sd(embedding_value),
      .groups = "drop"
    ) %>%
    dplyr::arrange(component_method, feature, patch_id, bin_center_px)
}
readr::write_csv(patch_bin, file.path(out_dir, "feature_profile_by_patch_bin.csv"))

profile_summary <- dplyr::bind_rows(
  summarise_profile_set(patch_bin, patches, "primary_cross_boundary", function(x) x$primary_in_analysis),
  summarise_profile_set(patch_bin, patches, "secondary_one_sided", function(x) x$profile_class == "one_sided" & x$n_profile_spots >= min_profile_spots),
  summarise_profile_set(patch_bin, patches, "secondary_no_zero", function(x) x$profile_class == "no_zero" & x$n_profile_spots >= min_profile_spots)
)
if (nrow(profile_summary) == 0) {
  profile_summary <- data.frame(
    analysis_set = character(),
    component_method = character(),
    feature_key = character(),
    feature = character(),
    distance_bin_index = integer(),
    bin_start_px = numeric(),
    bin_end_px = numeric(),
    bin_center_px = numeric(),
    n_patches = integer(),
    n_patch_bins = integer(),
    total_spots = integer(),
    mean_score = numeric(),
    median_score = numeric(),
    sd_score = numeric(),
    se_score = numeric(),
    ci_low = numeric(),
    ci_high = numeric(),
    stringsAsFactors = FALSE
  )
}
readr::write_csv(profile_summary, file.path(out_dir, "feature_profile_summary.csv"))

patch_bin_primary <- patch_bin[patch_bin$primary_in_analysis, , drop = FALSE]
boundary_tests <- make_band_tests(patch_bin_primary, bin_width_px)
readr::write_csv(boundary_tests, file.path(out_dir, "feature_boundary_tests.csv"))

primary_summary <- profile_summary[profile_summary$analysis_set == "primary_cross_boundary", , drop = FALSE]
if (!identical(plot_mode, "none") && nrow(primary_summary) > 0) {
  combos <- unique(primary_summary[, c("component_method", "feature"), drop = FALSE])
  for (i in seq_len(nrow(combos))) {
    method_id <- combos$component_method[i]
    feature_id <- combos$feature[i]
    plot_df <- primary_summary[
      primary_summary$component_method == method_id & primary_summary$feature == feature_id,
      ,
      drop = FALSE
    ]
    method_plot_dir <- file.path(profile_plot_dir, sanitize_filename(method_id))
    dir.create(method_plot_dir, recursive = TRUE, showWarnings = FALSE)
    plot_distance_profile(
      plot_df,
      file.path(method_plot_dir, paste0(sanitize_filename(feature_id), "_distance_profile.pdf")),
      paste(sample_name, method_id, feature_id, "tumor-arrow-guided profile")
    )
  }

  heatmap_df <- primary_summary
  heatmap_df$feature_key <- paste(heatmap_df$component_method, heatmap_df$feature, sep = "::")
  feature_order <- unique(heatmap_df$feature_key)
  bin_index_range <- seq(
    min(heatmap_df$distance_bin_index, na.rm = TRUE),
    max(heatmap_df$distance_bin_index, na.rm = TRUE)
  )
  heatmap_grid <- expand.grid(
    feature_key = feature_order,
    distance_bin_index = bin_index_range,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  feature_meta <- heatmap_df[match(feature_order, heatmap_df$feature_key), c("feature_key", "component_method", "feature"), drop = FALSE]
  heatmap_grid <- merge(heatmap_grid, feature_meta, by = "feature_key", all.x = TRUE, sort = FALSE)
  heatmap_grid$bin_center_px <- heatmap_grid$distance_bin_index * bin_width_px
  heatmap_grid$bin_start_px <- heatmap_grid$bin_center_px - bin_width_px / 2
  heatmap_grid$bin_end_px <- heatmap_grid$bin_center_px + bin_width_px / 2
  heatmap_values <- heatmap_df[, c("feature_key", "distance_bin_index", "mean_score", "n_patches"), drop = FALSE]
  heatmap_complete <- merge(
    heatmap_grid,
    heatmap_values,
    by = c("feature_key", "distance_bin_index"),
    all.x = TRUE,
    sort = FALSE
  )
  heatmap_complete$.feature_order <- match(heatmap_complete$feature_key, feature_order)
  heatmap_complete <- heatmap_complete[order(heatmap_complete$.feature_order, heatmap_complete$distance_bin_index), , drop = FALSE]
  heatmap_complete$z_score <- ave(
    heatmap_complete$mean_score,
    heatmap_complete$feature_key,
    FUN = function(x) {
      out <- rep(NA_real_, length(x))
      ok <- is.finite(x)
      if (!any(ok)) return(out)
      x_sd <- stats::sd(x[ok], na.rm = TRUE)
      if (!is.finite(x_sd) || x_sd == 0) {
        out[ok] <- 0
      } else {
        out[ok] <- (x[ok] - mean(x[ok], na.rm = TRUE)) / x_sd
      }
      out
    }
  )
  heatmap_complete$feature_key <- factor(heatmap_complete$feature_key, levels = rev(feature_order))
  p_heatmap <- ggplot2::ggplot(heatmap_complete, ggplot2::aes(x = bin_center_px, y = feature_key, fill = z_score)) +
    ggplot2::geom_tile(width = bin_width_px, height = 0.9, color = "white", linewidth = 0.1) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey30") +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, na.value = "grey92") +
    ggplot2::labs(
      title = paste(sample_name, "tumor-arrow-guided boundary profile heatmap"),
      x = "Aligned distance along tumor arrow axis (lowres px)",
      y = "Feature",
      fill = "Row z"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  ggplot2::ggsave(file.path(plot_dir, "feature_distance_profile_heatmap.pdf"), p_heatmap, width = 7.5, height = max(4, 0.22 * length(feature_order) + 2))
}

qc_counts <- as.data.frame(table(patches$profile_class), stringsAsFactors = FALSE)
colnames(qc_counts) <- c("profile_class", "n_patches")
qc_counts$profile_class <- as.character(qc_counts$profile_class)
qc_counts$profile_class_label <- profile_class_label_factor(qc_counts$profile_class)
p_qc <- ggplot2::ggplot(qc_counts, ggplot2::aes(x = profile_class_label, y = n_patches, fill = profile_class_label)) +
  ggplot2::geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = n_patches), vjust = -0.25, size = 3) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(
    title = paste(sample_name, "tumor-arrow profile patch QC"),
    x = "Profile class",
    y = "Patches"
  ) +
  ggplot2::theme_bw(base_size = 10)
if (!identical(plot_mode, "none")) {
  ggplot2::ggsave(file.path(plot_dir, "patch_qc_profile_class.pdf"), p_qc, width = 5.5, height = 4)
}

if (plot_space) {
  if (file.exists(spot_membership_path)) {
    all_spots <- read_csv_df(spot_membership_path)
    if (all(c("cell_ID", "X", "Y", "Location") %in% colnames(all_spots))) {
      all_spots <- unique(all_spots[, c("cell_ID", "X", "Y", "Location"), drop = FALSE])
    } else {
      all_spots <- unique(local_df[, c("cell_ID", "X", "Y", "Location"), drop = FALSE])
    }
  } else {
    all_spots <- unique(local_df[, c("cell_ID", "X", "Y", "Location"), drop = FALSE])
  }
  all_spots$X <- as.numeric(all_spots$X)
  all_spots$Y <- as.numeric(all_spots$Y)
  all_spots$Location <- factor(as.character(all_spots$Location), levels = names(boundary_cols))
  plot_patches <- patches[
    is.finite(patches$arrow_unit_vx) & is.finite(patches$arrow_unit_vy) &
      is.finite(patches$arrow_vx_u) & is.finite(patches$arrow_vy_u),
    ,
    drop = FALSE
  ]
  component_levels <- sort(unique(as.character(plot_patches$tumor_component)))
  component_palette <- stats::setNames(
    grDevices::hcl.colors(length(component_levels), palette = "Dark 3"),
    component_levels
  )
  plot_arrows <- make_profile_plot_arrow_df(plot_patches, component_levels, arrow_length_scale)
  zero_points <- plot_patches[plot_patches$zero_found, , drop = FALSE]
  zero_points_plot <- zero_points
  zero_points_plot$zero_axis_Y <- -zero_points_plot$zero_axis_Y
  profile_points <- unique(profile_spots[profile_spots$in_profile, c("cell_ID", "X", "Y", "Location"), drop = FALSE])
  profile_points$X <- as.numeric(profile_points$X)
  profile_points$Y <- as.numeric(profile_points$Y)
  profile_spot_polygons <- sp_make_spot_polygons(
    profile_points,
    plot_context,
    segments = params$spatial_plot_spot_segments %||% 24L
  )

  p_space <- make_boundary_spot_base(
    all_spots,
    paste(sample_name, "tumor-arrow-guided profile patches")
  ) +
    sp_spot_layer(
      profile_spot_polygons,
      fill = NA,
      colour = "grey5",
      alpha = 0.6,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = zero_points_plot,
      ggplot2::aes(x = zero_axis_X, y = zero_axis_Y),
      inherit.aes = FALSE,
      shape = 4,
      size = 1.7,
      stroke = 0.5,
      color = "#de2d26"
    )
  p_space <- add_profile_arrows(
    p_space,
    plot_arrows,
    component_palette,
    legend_title = "cell_component"
  )
  sp_save_plot(
    p_space,
    file.path(plot_dir, "tumor_arrow_profile_patches_space.pdf"),
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    map_width = params$spatial_plot_map_width %||% 6.45,
    legend_width = params$spatial_plot_legend_width %||% 1.55
  )

  plot_rectangles <- patch_rectangles
  plot_rectangles$Y <- -plot_rectangles$Y
  plot_rectangles$profile_class_label <- profile_class_label_factor(plot_rectangles$profile_class)
  rectangle_cols <- c(
    "cross-boundary\nzero found" = "#08519c",
    "zero found\none side only" = "#756bb1",
    "no zero\nfound" = "#636363",
    "no direction" = "#bdbdbd",
    "no local spots" = "#bdbdbd"
  )
  p_rectangles <- make_boundary_spot_base(
    all_spots,
    paste(sample_name, "tumor-arrow profile tube coverage")
  )
  if (nrow(plot_rectangles) > 0) {
    excluded_rectangles <- plot_rectangles[!plot_rectangles$primary_in_analysis, , drop = FALSE]
    primary_rectangles <- plot_rectangles[plot_rectangles$primary_in_analysis, , drop = FALSE]
    if (nrow(excluded_rectangles) > 0) {
      p_rectangles <- p_rectangles +
        ggplot2::geom_polygon(
          data = excluded_rectangles,
          ggplot2::aes(x = X, y = Y, group = patch_id),
          inherit.aes = FALSE,
          color = "grey45",
          fill = NA,
          linewidth = 0.3,
          linetype = "dashed",
          alpha = 0.35,
          na.rm = TRUE
        )
    }
    if (nrow(primary_rectangles) > 0) {
      p_rectangles <- p_rectangles +
      ggplot2::geom_polygon(
        data = primary_rectangles,
        ggplot2::aes(x = X, y = Y, group = patch_id, color = profile_class_label),
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 0.35,
        alpha = 0.75,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = rectangle_cols, drop = TRUE) +
      ggplot2::labs(color = "Profile class")
    }
  }
  p_rectangles <- p_rectangles +
    sp_spot_layer(
      profile_spot_polygons,
      fill = NA,
      colour = "grey5",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = zero_points_plot,
      ggplot2::aes(x = zero_axis_X, y = zero_axis_Y),
      inherit.aes = FALSE,
      shape = 4,
      size = 1.7,
      stroke = 0.5,
      color = "#de2d26",
      na.rm = TRUE
    )
  p_rectangles <- add_profile_arrows(
    p_rectangles,
    plot_arrows,
    fixed_color = "grey20"
  )
  sp_save_plot(
    p_rectangles,
    file.path(plot_dir, "tumor_arrow_profile_tube_rectangles_space.pdf"),
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    map_width = params$spatial_plot_map_width %||% 6.45,
    legend_width = params$spatial_plot_legend_width %||% 1.55
  )
}

summary_path <- file.path(out_dir, "run_summary_06.txt")
sink(summary_path)
cat("Tumor-arrow-guided boundary profile\n")
cat("Sample:", sample_name, "\n")
cat("Strategy:", strategy_id, "\n")
cat("Tumor components:", paste(tumor_components, collapse = ","), "\n")
cat("Profile methods:", paste(profile_methods, collapse = ","), "\n")
cat("Positive control features:", paste(positive_control_features, collapse = ","), "\n")
cat("Requested features:", ifelse(feature_filter_all, "all", paste(profile_features, collapse = ",")), "\n")
cat("Selected feature rows:", nrow(feature_selection[feature_selection$selected, , drop = FALSE]), "\n")
cat("Missing requested features:", length(missing_requested), "\n")
cat("Missing positive control features:", length(missing_positive_controls), "\n")
cat("Arrow table:", arrow_path, "\n")
cat("Local membership:", local_membership_path, "\n")
cat("Scalefactors:", dims$scalefactors_path, "\n")
cat("Tissue positions:", dims$tissue_positions_path, "\n")
cat("spot_diameter_px:", dims$spot_diameter_px, "\n")
cat("spot_radius_px:", dims$spot_radius_px, "\n")
cat("spot_pitch_px:", dims$spot_pitch_px, "\n")
cat("arrow_length_scale:", arrow_length_scale, "\n")
cat("spatial_canvas_inches:", params$spatial_plot_width %||% 8, "x", params$spatial_plot_height %||% 7, "\n")
cat("tube_half_width_mode:", tube_half_width_mode, "\n")
cat("tube_half_width_px:", tube_half_width_px, "\n")
cat("zero_width_modes:", paste(zero_width_modes, collapse = ","), "\n")
cat("zero_width_px:", paste(round(zero_width_px, 6), collapse = ","), "\n")
cat("bin_width_mode:", bin_width_mode, "\n")
cat("bin_width_px:", bin_width_px, "\n")
cat("bin_assignment:", "centered_on_zero", "\n")
cat("bdy_cluster_gap_mode:", bdy_cluster_gap_mode, "\n")
cat("bdy_cluster_gap_px:", bdy_cluster_gap_px, "\n")
cat("min_profile_spots:", min_profile_spots, "\n")
cat("require_normal_spot:", require_normal_spot, "\n")
cat("normal_location:", normal_location, "\n")
cat("Tumor arrows:", nrow(tumor_arrows), "\n")
cat("Patches:", nrow(patches), "\n")
cat("Patch rectangles:", length(unique(patch_rectangles$patch_id)), "\n")
cat("Zero found:", sum(patches$zero_found, na.rm = TRUE), "\n")
cat("Primary patches:", sum(patches$primary_in_analysis, na.rm = TRUE), "\n")
cat("Patches passing min spots:", sum(patches$passes_min_spots_filter, na.rm = TRUE), "\n")
cat("Patches passing normal filter:", sum(patches$passes_normal_filter, na.rm = TRUE), "\n")
cat("Profile class labels:\n")
print(data.frame(profile_class = profile_class_order, label = unname(profile_class_labels[profile_class_order])))
cat("Profile class counts:\n")
print(table(patches$profile_class, useNA = "ifany"))
cat("Profile exclusion reason counts:\n")
print(table(patches$profile_exclusion_reason, useNA = "ifany"))
cat("Profile spots in profile:", sum(profile_spots$in_profile, na.rm = TRUE), "\n")
cat("Median profile spots per patch:", stats::median(patches$n_profile_spots, na.rm = TRUE), "\n")
cat("\nOutputs:\n")
cat("- tumor_arrow_profile_patches.csv\n")
cat("- tumor_arrow_profile_spots.csv\n")
cat("- tumor_arrow_profile_rectangles.csv\n")
cat("- selected_features.csv\n")
cat("- missing_requested_features.csv\n")
cat("- missing_positive_control_features.csv\n")
cat("- feature_profile_by_patch_bin.csv\n")
cat("- feature_profile_summary.csv\n")
cat("- feature_boundary_tests.csv\n")
cat("- plots/\n")
sink()

cat("Tumor-arrow-guided boundary profile completed.\n")
cat("Output:", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Tumor arrows:", nrow(tumor_arrows), "\n")
cat("Zero found:", sum(patches$zero_found, na.rm = TRUE), "\n")
cat("Primary patches:", sum(patches$primary_in_analysis, na.rm = TRUE), "\n")
