# 05：筛选边界相关 LSGI arrows，并按 LSGI 多 module/method 输出格式绘图。
# 输入：
# - output/<sample>/11_lsgi_gradient/<module>/arrow_tables/arrows_by_grid.csv
# - output/<sample>/11_lsgi_gradient/grid_partition_membership.csv
# - output/<sample>/11_lsgi_gradient/grid_local_spot_membership.csv
# - output/<sample>/07_spatial_deconvolution/spot_matrix_pre_lsgi.csv：用于补充 spot 反卷积信息。
# 输出：
# - output/<sample>/12_boundary_related_lsgi_arrows/<module>/strategy_summary.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/<module>/arrow_tables/*_arrows.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/<module>/grid_tables/*_selected_grids.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/<module>/spot_tables/*_bdy_spot_arrow_assignments.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/<module>/plots/*BoundaryRelated*.pdf
# - output/<sample>/12_boundary_related_lsgi_arrows/strategy_summary.csv 汇总所有 module。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "R", "boundarygrad_core.R"))
source(file.path(script_dir, "R", "spatial_plot_core.R"))
load_required_packages(c("readr", "dplyr", "tibble", "ggplot2", "png", "grid", "jsonlite", "patchwork"))

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

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path, call. = FALSE)
  }
  path
}

coalesce_chr <- function(x, fallback = "") {
  ifelse(is.na(x), fallback, x)
}

cli_opts <- parse_cli_options(commandArgs(trailingOnly = TRUE))
plot_mode <- bg_validate_plot_mode(get_opt(cli_opts, "plot-mode", default = params$plot_mode %||% "full"))

lsgi_dir <- file.path(paths$output, "11_lsgi_gradient")
out_dir <- file.path(paths$output, "12_boundary_related_lsgi_arrows")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
legacy_root_output_dirs <- file.path(out_dir, c("arrow_tables", "grid_tables", "spot_tables", "plots", "plots_by_strategy"))
invisible(lapply(legacy_root_output_dirs, bg_assert_safe_child, root = paths$output,
                 label = "legacy output cleanup target"))
unlink(legacy_root_output_dirs, recursive = TRUE, force = TRUE)

component_methods <- parse_chr_vec(
  get_opt(cli_opts, "component-methods", default = params$lsgi_component_methods %||% "cell_component,nmf,marker_module,pathway,single_gene"),
  default = c("cell_component", "nmf", "marker_module", "pathway", "single_gene")
)
if (length(component_methods) == 1 && identical(tolower(component_methods), "all")) {
  method_summary_path <- file.path(lsgi_dir, "component_method_summary.csv")
  if (file.exists(method_summary_path)) {
    component_methods <- readr::read_csv(method_summary_path, show_col_types = FALSE)$method_id
  } else {
    method_dirs <- list.dirs(lsgi_dir, recursive = FALSE, full.names = FALSE)
    component_methods <- method_dirs[file.exists(file.path(lsgi_dir, method_dirs, "arrow_tables", "arrows_by_grid.csv"))]
  }
}
component_methods <- unique(component_methods)
if (length(component_methods) == 0) {
  stop("No LSGI component methods requested or discoverable.", call. = FALSE)
}

partition_path <- require_file(file.path(lsgi_dir, "grid_partition_membership.csv"), "grid partition membership")
local_path <- require_file(file.path(lsgi_dir, "grid_local_spot_membership.csv"), "grid local spot membership")
grid_info_path <- require_file(file.path(lsgi_dir, "grid_info.csv"), "grid info")
spot_matrix_candidates <- c(
  file.path(paths$output, "07", "spot_matrix_pre_lsgi.csv"),
  file.path(paths$output, "07_spatial_deconvolution", "spot_matrix_pre_lsgi.csv")
)
existing_spot_matrix <- spot_matrix_candidates[file.exists(spot_matrix_candidates)]
spot_matrix_path <- if (length(existing_spot_matrix) > 0) existing_spot_matrix[[1]] else spot_matrix_candidates[[1]]

partition_df <- readr::read_csv(partition_path, show_col_types = FALSE)
local_df <- readr::read_csv(local_path, show_col_types = FALSE)
grid_info_df <- readr::read_csv(grid_info_path, show_col_types = FALSE)

spot_matrix_07 <- NULL
spot_decon_07 <- NULL
if (file.exists(spot_matrix_path)) {
  spot_matrix_07 <- readr::read_csv(spot_matrix_path, show_col_types = FALSE)
  decon_cols <- setdiff(colnames(spot_matrix_07), c("cell_ID", "X", "Y", "Location"))
  spot_decon_07 <- spot_matrix_07[, c("cell_ID", decon_cols), drop = FALSE]
} else {
  warning(
    "07 spot matrix not found in output; spot-arrow assignment tables will not include deconvolution columns. Tried: ",
    paste(spot_matrix_candidates, collapse = ", "),
    call. = FALSE
  )
}

include_unfiltered_arrows <- as_bool(get_opt(cli_opts, "include-unfiltered-arrows", default = FALSE))
arrow_length_scale <- as.numeric(get_opt(cli_opts, "arrow-length-scale", default = params$lsgi_arrow_length_scale %||% 1.4))
arrow_linewidth <- as.numeric(get_opt(cli_opts, "arrow-linewidth", default = params$lsgi_arrow_linewidth %||% 1.0))
arrow_head_cm <- as.numeric(get_opt(cli_opts, "arrow-head-cm", default = params$lsgi_arrow_head_cm %||% 0.20))
arrow_head_angle <- as.numeric(get_opt(cli_opts, "arrow-head-angle", default = params$lsgi_arrow_head_angle %||% 30))
arrow_closed <- as_bool(get_opt(cli_opts, "arrow-closed", default = params$lsgi_arrow_closed %||% TRUE))
plot_nonselected_arrows <- as_bool(get_opt(cli_opts, "plot-nonselected-arrows", default = TRUE))
strategy_filter <- strsplit(as.character(get_opt(cli_opts, "strategy", default = params$boundary_lsgi_arrow_strategy %||% "local_broad")), ",", fixed = TRUE)[[1]]
strategy_filter <- trimws(strategy_filter[nzchar(strategy_filter)])
arrow_type <- if (arrow_closed) "closed" else "open"
if (!is.finite(arrow_head_angle) || arrow_head_angle <= 0 || arrow_head_angle >= 180) {
  stop("--arrow-head-angle must be a finite value between 0 and 180 degrees.", call. = FALSE)
}

required_arrow_cols <- c("grid", "component", "X", "Y", "vx.u", "vy.u", "included_in_filtered_output")
boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")

read_method_arrows <- function(method_id) {
  arrow_path <- require_file(
    file.path(lsgi_dir, method_id, "arrow_tables", "arrows_by_grid.csv"),
    paste0("LSGI arrow table for module ", method_id)
  )
  arrows <- readr::read_csv(arrow_path, show_col_types = FALSE)
  missing_arrow_cols <- setdiff(required_arrow_cols, colnames(arrows))
  if (length(missing_arrow_cols) > 0) {
    stop(
      "Missing columns in arrow table for module ", method_id, ": ",
      paste(missing_arrow_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"arrow_head_angle_mapped" %in% colnames(arrows)) {
    arrows$arrow_head_angle_mapped <- arrow_head_angle
  }
  arrows %>%
    dplyr::mutate(
      component_method = method_id,
      component_key = paste(method_id, component, sep = "::"),
      arrow_source_path = arrow_path
    )
}

make_grid_boundary_summary <- function(membership_df, grid_col, method_name) {
  membership_df %>%
    dplyr::rename(grid = dplyr::all_of(grid_col)) %>%
    dplyr::mutate(Location = factor(Location, levels = c("Mal", "Bdy", "nMal"))) %>%
    dplyr::group_by(grid) %>%
    dplyr::summarise(
      method = method_name,
      n_total = dplyr::n(),
      n_Mal = sum(Location == "Mal", na.rm = TRUE),
      n_Bdy = sum(Location == "Bdy", na.rm = TRUE),
      n_nMal = sum(Location == "nMal", na.rm = TRUE),
      frac_Bdy = n_Bdy / n_total,
      .groups = "drop"
    ) %>%
    dplyr::right_join(grid_info_df[, c("grid", "grid_index", "X", "Y"), drop = FALSE], by = "grid") %>%
    dplyr::mutate(
      method = coalesce_chr(method, method_name),
      n_total = ifelse(is.na(n_total), 0L, n_total),
      n_Mal = ifelse(is.na(n_Mal), 0L, n_Mal),
      n_Bdy = ifelse(is.na(n_Bdy), 0L, n_Bdy),
      n_nMal = ifelse(is.na(n_nMal), 0L, n_nMal),
      frac_Bdy = ifelse(is.na(frac_Bdy), 0, frac_Bdy)
    ) %>%
    dplyr::arrange(grid_index)
}

partition_summary <- make_grid_boundary_summary(partition_df, "nearest_grid", "partition")
local_summary <- make_grid_boundary_summary(local_df, "grid", "local")

grid_metrics <- partition_summary %>%
  dplyr::select(
    grid,
    grid_index,
    grid_center_X = X,
    grid_center_Y = Y,
    partition_n_total = n_total,
    partition_n_Mal = n_Mal,
    partition_n_Bdy = n_Bdy,
    partition_n_nMal = n_nMal,
    partition_frac_Bdy = frac_Bdy
  ) %>%
  dplyr::left_join(
    local_summary %>%
      dplyr::select(
        grid,
        local_n_total = n_total,
        local_n_Mal = n_Mal,
        local_n_Bdy = n_Bdy,
        local_n_nMal = n_nMal,
        local_frac_Bdy = frac_Bdy
      ),
    by = "grid"
  )
readr::write_csv(grid_metrics, file.path(out_dir, "grid_boundary_metrics.csv"))

# 内置试验策略：
# partition_*：grid 的 nearest partition 区域中 Bdy spot 数量/比例达标。
# local_*：grid 的 LSGI local regression spots 中 Bdy spot 数量/比例达标。
# consensus/union：同时或任一满足推荐的 partition/local 主阈值。
strategies <- tibble::tribble(
  ~strategy_id, ~strategy_type, ~min_partition_n_Bdy, ~min_partition_frac_Bdy, ~min_local_n_Bdy, ~min_local_frac_Bdy, ~description,
  "partition_any", "partition", 1, 0.00, NA_real_, NA_real_, "Any partition Bdy spot under the grid; intentionally broad.",
  "partition_relaxed", "partition", 2, 0.20, NA_real_, NA_real_, "Relaxed partition boundary grid.",
  "partition_primary", "partition", 3, 0.30, NA_real_, NA_real_, "Recommended CRC1 main definition: grid territory enriched for Bdy spots.",
  "partition_strict", "partition", 5, 0.40, NA_real_, NA_real_, "Strict partition boundary grid.",
  "local_broad", "local", NA_real_, NA_real_, 5, 0.10, "Broad local-regression definition.",
  "local_relaxed", "local", NA_real_, NA_real_, 8, 0.15, "Relaxed local-regression definition.",
  "local_primary", "local", NA_real_, NA_real_, 10, 0.20, "Recommended CRC1 local-regression definition.",
  "local_strict", "local", NA_real_, NA_real_, 12, 0.25, "Strict local-regression definition.",
  "primary_consensus", "intersection", 3, 0.30, 10, 0.20, "Grid passes both partition_primary and local_primary.",
  "primary_union", "union", 3, 0.30, 10, 0.20, "Grid passes either partition_primary or local_primary."
)

if (!identical(strategy_filter, "all")) {
  missing_strategies <- setdiff(strategy_filter, strategies$strategy_id)
  if (length(missing_strategies) > 0) {
    stop("Unknown strategy: ", paste(missing_strategies, collapse = ", "), call. = FALSE)
  }
  strategies <- strategies[strategies$strategy_id %in% strategy_filter, , drop = FALSE]
}

passes_partition <- bg_passes_partition
passes_local <- bg_passes_local
selected_grid_flags <- bg_selected_grid_flags

make_plot_arrow_df <- function(df, component_levels) {
  if (!"arrow_head_angle_mapped" %in% colnames(df)) {
    df$arrow_head_angle_mapped <- arrow_head_angle
  }
  out <- sp_prepare_arrows(
    df,
    arrow_length_scale = arrow_length_scale,
    flip_y = TRUE
  )
  out %>%
    dplyr::mutate(
      arrow_head_angle_mapped = ifelse(is.finite(arrow_head_angle_mapped), arrow_head_angle_mapped, arrow_head_angle),
      component = factor(component, levels = component_levels)
    )
}

add_selected_arrows <- function(p, selected_arrow_df, background_arrow_df = NULL, component_palette, method_label, color_drop = FALSE) {
  out <- p
  if (plot_nonselected_arrows && !is.null(background_arrow_df) && nrow(background_arrow_df) > 0) {
    background_angles <- sort(unique(background_arrow_df$arrow_head_angle_mapped[is.finite(background_arrow_df$arrow_head_angle_mapped)]))
    if (length(background_angles) == 0) {
      background_angles <- arrow_head_angle
    }
    for (angle_value in background_angles) {
      layer_data <- background_arrow_df[is.finite(background_arrow_df$arrow_head_angle_mapped) & background_arrow_df$arrow_head_angle_mapped == angle_value, , drop = FALSE]
      if (nrow(layer_data) == 0) next
      out <- out +
        ggplot2::geom_segment(
          data = layer_data,
          ggplot2::aes(x = X, y = Y, xend = X_end, yend = Y_end),
          inherit.aes = FALSE,
          color = "grey70",
          alpha = 0.25,
          linewidth = arrow_linewidth,
          lineend = "round",
          na.rm = TRUE,
          arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), angle = angle_value, type = arrow_type)
        )
    }
  }
  if (nrow(selected_arrow_df) == 0) {
    return(out + ggplot2::labs(subtitle = "No boundary-related LSGI arrows under this strategy"))
  }

  selected_angles <- sort(unique(selected_arrow_df$arrow_head_angle_mapped[is.finite(selected_arrow_df$arrow_head_angle_mapped)]))
  if (length(selected_angles) == 0) {
    selected_angles <- arrow_head_angle
  }
  for (angle_value in selected_angles) {
    layer_data <- selected_arrow_df[is.finite(selected_arrow_df$arrow_head_angle_mapped) & selected_arrow_df$arrow_head_angle_mapped == angle_value, , drop = FALSE]
    if (nrow(layer_data) == 0) next
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
  }
  out +
    ggplot2::scale_color_manual(values = component_palette, drop = color_drop) +
    ggplot2::labs(color = method_label)
}

point_df <- partition_df %>%
  dplyr::select(cell_ID, X, Y, Location) %>%
  dplyr::distinct(cell_ID, .keep_all = TRUE) %>%
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))

plot_context <- sp_read_context(paths$spaceranger, load_image = !identical(plot_mode, "none"))
has_image <- plot_context$has_image
point_spot_polygons <- if (!identical(plot_mode, "none")) {
  sp_make_spot_polygons(
    point_df,
    plot_context,
    segments = params$spatial_plot_spot_segments %||% 24L
  )
} else NULL

save_boundary_arrow_plot <- function(plot, filename) {
  sp_save_plot(
    plot,
    filename,
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    map_width = params$spatial_plot_map_width %||% 6.45,
    legend_width = params$spatial_plot_legend_width %||% 1.55
  )
}

process_method <- function(method_id) {
  message("Screening boundary-related LSGI arrows for module: ", method_id)
  module_out_dir <- file.path(out_dir, method_id)
  arrow_out_dir <- file.path(module_out_dir, "arrow_tables")
  grid_out_dir <- file.path(module_out_dir, "grid_tables")
  spot_out_dir <- file.path(module_out_dir, "spot_tables")
  plot_out_dir <- file.path(module_out_dir, "plots")
  output_dirs <- c(module_out_dir, arrow_out_dir, grid_out_dir, spot_out_dir)
  if (!identical(plot_mode, "none")) output_dirs <- c(output_dirs, plot_out_dir)
  invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  arrows_all <- read_method_arrows(method_id)
  candidate_arrows <- arrows_all
  if (!include_unfiltered_arrows) {
    candidate_arrows <- candidate_arrows[candidate_arrows$included_in_filtered_output, , drop = FALSE]
  }

  component_levels <- sort(unique(arrows_all$component))
  component_palette <- setNames(
    grDevices::hcl.colors(length(component_levels), palette = "Dark 3"),
    component_levels
  )
  background_arrows <- if (!identical(plot_mode, "none")) {
    make_plot_arrow_df(candidate_arrows, component_levels)
  } else NULL
  readr::write_csv(grid_metrics, file.path(module_out_dir, "grid_boundary_metrics.csv"))

  strategy_summary_rows <- vector("list", nrow(strategies))
  for (i in seq_len(nrow(strategies))) {
    strategy <- strategies[i, , drop = FALSE]
    strategy_id <- strategy$strategy_id[[1]]
    grid_flags <- selected_grid_flags(grid_metrics, strategy)
    selected_grids <- grid_metrics[grid_flags, , drop = FALSE]

    selected_arrows <- candidate_arrows %>%
      dplyr::filter(grid %in% selected_grids$grid) %>%
      dplyr::left_join(grid_metrics, by = "grid") %>%
      dplyr::mutate(
        strategy_id = strategy_id,
        strategy_type = strategy$strategy_type[[1]],
        strategy_description = strategy$description[[1]]
      ) %>%
      dplyr::select(strategy_id, strategy_type, strategy_description, component_method, component_key, dplyr::everything())

    selected_grids_out <- selected_grids %>%
      dplyr::mutate(
        component_method = method_id,
        strategy_id = strategy_id,
        strategy_type = strategy$strategy_type[[1]],
        strategy_description = strategy$description[[1]]
      ) %>%
      dplyr::select(component_method, strategy_id, strategy_type, strategy_description, dplyr::everything())

    readr::write_csv(selected_arrows, file.path(arrow_out_dir, paste0(strategy_id, "_arrows.csv")))
    readr::write_csv(selected_grids_out, file.path(grid_out_dir, paste0(strategy_id, "_selected_grids.csv")))

    partition_spot_arrow_base <- partition_df %>%
      dplyr::filter(Location == "Bdy", nearest_grid %in% selected_grids$grid) %>%
      dplyr::rename(source_grid = nearest_grid) %>%
      dplyr::mutate(assignment_method = "partition")
    local_spot_arrow_base <- local_df %>%
      dplyr::filter(Location == "Bdy", grid %in% selected_grids$grid) %>%
      dplyr::rename(source_grid = grid) %>%
      dplyr::mutate(assignment_method = "local")
    if (strategy$strategy_type[[1]] == "partition") {
      spot_arrow_base <- partition_spot_arrow_base
    } else if (strategy$strategy_type[[1]] %in% c("local", "intersection")) {
      spot_arrow_base <- local_spot_arrow_base
    } else if (strategy$strategy_type[[1]] == "union") {
      spot_arrow_base <- dplyr::bind_rows(partition_spot_arrow_base, local_spot_arrow_base) %>%
        dplyr::distinct(assignment_method, source_grid, cell_ID, .keep_all = TRUE)
    } else {
      stop("Unsupported strategy_type for spot assignment: ", strategy$strategy_type[[1]], call. = FALSE)
    }
    spot_arrow_join <- selected_arrows %>%
      dplyr::rename(source_grid = grid)

    spot_arrow_assignments <- spot_arrow_base %>%
      dplyr::inner_join(
        spot_arrow_join,
        by = "source_grid",
        suffix = c("_spot", "_arrow"),
        relationship = "many-to-many"
      ) %>%
      dplyr::mutate(strategy_id = strategy_id, component_method = method_id) %>%
      dplyr::select(component_method, strategy_id, dplyr::everything())
    if (!is.null(spot_decon_07)) {
      spot_arrow_assignments <- spot_arrow_assignments %>%
        dplyr::left_join(spot_decon_07, by = "cell_ID")
    }
    readr::write_csv(spot_arrow_assignments, file.path(spot_out_dir, paste0(strategy_id, "_bdy_spot_arrow_assignments.csv")))

    if (!identical(plot_mode, "none")) {
    selected_plot_arrows <- make_plot_arrow_df(selected_arrows, component_levels)
    background_plot_arrows <- background_arrows %>%
      dplyr::filter(!grid %in% selected_grids$grid)

    boundary_base <- sp_boundary_base(
      point_spot_polygons,
      plot_context,
      boundary_cols,
      title = paste0(sample_name, " ", method_id, " boundary-related LSGI arrows: ", strategy_id),
      show_image = FALSE
    )

    boundary_plot <- add_selected_arrows(
      boundary_base,
      selected_plot_arrows,
      background_plot_arrows,
      component_palette,
      method_id
    )
    save_boundary_arrow_plot(
      boundary_plot,
      file.path(plot_out_dir, paste0(sample_name, "_", method_id, "_", strategy_id, "_BoundaryRelated_LSGIGradient.pdf"))
    )

    if (has_image) {
      he_base <- sp_boundary_base(
        point_spot_polygons,
        plot_context,
        boundary_cols,
        title = paste0(sample_name, " HE ", method_id, " boundary-related LSGI arrows: ", strategy_id),
        show_image = TRUE
      )

      he_plot <- add_selected_arrows(
        he_base,
        selected_plot_arrows,
        background_plot_arrows,
        component_palette,
        method_id
      )
      save_boundary_arrow_plot(
        he_plot,
        file.path(plot_out_dir, paste0(sample_name, "_", method_id, "_", strategy_id, "_BoundaryRelated_HE_LSGIGradient.pdf"))
      )
    }

    component_files <- sanitize_filename(component_levels, prefix = "component")
    if (identical(plot_mode, "full")) for (component_idx in seq_along(component_levels)) {
      component_id <- component_levels[[component_idx]]
      component_file <- component_files[[component_idx]]
      selected_component_arrows <- selected_plot_arrows %>%
        dplyr::filter(as.character(component) == component_id)
      background_component_arrows <- background_plot_arrows %>%
        dplyr::filter(as.character(component) == component_id)

      component_boundary_plot <- add_selected_arrows(
        boundary_base + ggplot2::ggtitle(paste0(sample_name, " ", method_id, " boundary-related LSGI arrow: ", strategy_id, " / ", component_id)),
        selected_component_arrows,
        background_component_arrows,
        component_palette[component_id],
        method_id,
        color_drop = TRUE
      )
      save_boundary_arrow_plot(
        component_boundary_plot,
        file.path(plot_out_dir, paste0(sample_name, "_", method_id, "_", strategy_id, "_", component_file, "_BoundaryRelated_LSGIGradient.pdf"))
      )

      if (has_image) {
        component_he_plot <- add_selected_arrows(
          he_base + ggplot2::ggtitle(paste0(sample_name, " HE ", method_id, " boundary-related LSGI arrow: ", strategy_id, " / ", component_id)),
          selected_component_arrows,
          background_component_arrows,
          component_palette[component_id],
          method_id,
          color_drop = TRUE
        )
        save_boundary_arrow_plot(
          component_he_plot,
          file.path(plot_out_dir, paste0(sample_name, "_", method_id, "_", strategy_id, "_", component_file, "_BoundaryRelated_HE_LSGIGradient.pdf"))
        )
      }
    }
    }

    component_counts <- sort(table(selected_arrows$component), decreasing = TRUE)
    strategy_summary_rows[[i]] <- data.frame(
      component_method = method_id,
      strategy_id = strategy_id,
      strategy_type = strategy$strategy_type[[1]],
      min_partition_n_Bdy = strategy$min_partition_n_Bdy[[1]],
      min_partition_frac_Bdy = strategy$min_partition_frac_Bdy[[1]],
      min_local_n_Bdy = strategy$min_local_n_Bdy[[1]],
      min_local_frac_Bdy = strategy$min_local_frac_Bdy[[1]],
      n_selected_grids = nrow(selected_grids),
      n_candidate_arrows = nrow(candidate_arrows),
      n_selected_arrows = nrow(selected_arrows),
      n_selected_bdy_spot_arrow_rows = nrow(spot_arrow_assignments),
      components = paste(names(component_counts), as.integer(component_counts), sep = "=", collapse = ";"),
      description = strategy$description[[1]],
      stringsAsFactors = FALSE
    )
  }

  strategy_summary <- dplyr::bind_rows(strategy_summary_rows)
  readr::write_csv(strategy_summary, file.path(module_out_dir, "strategy_summary.csv"))

  sink(file.path(module_out_dir, "run_summary_12.txt"))
  cat("Boundary-related LSGI arrow strategy screen\n")
  cat("==========================================\n\n")
  cat("Sample:", sample_name, "\n")
  cat("Module/method:", method_id, "\n")
  cat("Total arrows:", nrow(arrows_all), "\n")
  cat("Candidate arrows:", nrow(candidate_arrows), "\n")
  cat("Include unfiltered arrows:", include_unfiltered_arrows, "\n")
  cat("Partition spots:", nrow(partition_df), "\n")
  cat("Local membership rows:", nrow(local_df), "\n")
  cat("07 spot matrix available:", !is.null(spot_matrix_07), "\n")
  cat("07 spot matrix path:", if (file.exists(spot_matrix_path)) spot_matrix_path else "not found", "\n")
  cat("H&E overlay available:", has_image, "\n")
  cat("Low-resolution spot diameter px:", plot_context$spot_diameter_lowres, "\n")
  cat("Arrow length scale:", arrow_length_scale, "\n")
  cat("Spatial canvas inches:", params$spatial_plot_width %||% 8, "x", params$spatial_plot_height %||% 7, "\n\n")
  cat("Strategies:\n")
  print(strategy_summary)
  cat("\nOutputs written to:\n")
  cat(module_out_dir, "\n")
  sink()

  strategy_summary
}

strategy_summary_all <- dplyr::bind_rows(lapply(component_methods, process_method))
readr::write_csv(strategy_summary_all, file.path(out_dir, "strategy_summary.csv"))

sink(file.path(out_dir, "run_summary_12.txt"))
cat("Boundary-related LSGI arrow strategy screen\n")
cat("==========================================\n\n")
cat("Sample:", sample_name, "\n")
cat("Modules/methods:", paste(component_methods, collapse = ", "), "\n")
cat("Include unfiltered arrows:", include_unfiltered_arrows, "\n")
cat("Partition spots:", nrow(partition_df), "\n")
cat("Local membership rows:", nrow(local_df), "\n")
cat("07 spot matrix available:", !is.null(spot_matrix_07), "\n")
cat("07 spot matrix path:", if (file.exists(spot_matrix_path)) spot_matrix_path else "not found", "\n")
cat("H&E overlay available:", has_image, "\n")
cat("Low-resolution spot diameter px:", plot_context$spot_diameter_lowres, "\n")
cat("Arrow length scale:", arrow_length_scale, "\n")
cat("Spatial canvas inches:", params$spatial_plot_width %||% 8, "x", params$spatial_plot_height %||% 7, "\n\n")
cat("Strategy summaries:\n")
print(strategy_summary_all)
cat("\nModule output directories:\n")
for (method_id in component_methods) {
  cat("- ", file.path(out_dir, method_id), "\n", sep = "")
}
sink()
