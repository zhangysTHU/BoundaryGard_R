# 12：试验性筛选边界相关 LSGI arrows，并按 11 的风格绘制边界/H&E 叠加图。
# 输入：
# - output/<sample>/11_lsgi_gradient/arrow_tables/cell_component_arrows_by_grid.csv
# - output/<sample>/11_lsgi_gradient/grid_partition_membership.csv
# - output/<sample>/11_lsgi_gradient/grid_local_spot_membership.csv
# - output/<sample>/07/spot_matrix_pre_lsgi.csv 或 output/<sample>/07_spatial_deconvolution/spot_matrix_pre_lsgi.csv：
#   07 现在写在 output 中；这里用于补充 spot 反卷积信息。
# 输出：
# - output/<sample>/12_boundary_related_lsgi_arrows/strategy_summary.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/arrow_tables/*_arrows.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/grid_tables/*_selected_grids.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/spot_tables/*_bdy_spot_arrow_assignments.csv
# - output/<sample>/12_boundary_related_lsgi_arrows/plots/*BoundaryRelated*.pdf
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("readr", "dplyr", "tibble", "ggplot2", "png", "grid"))

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

lsgi_dir <- file.path(paths$output, "11_lsgi_gradient")
out_dir <- file.path(paths$output, "12_boundary_related_lsgi_arrows")
arrow_out_dir <- file.path(out_dir, "arrow_tables")
grid_out_dir <- file.path(out_dir, "grid_tables")
spot_out_dir <- file.path(out_dir, "spot_tables")
plot_out_dir <- file.path(out_dir, "plots")
invisible(lapply(c(out_dir, arrow_out_dir, grid_out_dir, spot_out_dir, plot_out_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

arrow_path <- require_file(file.path(lsgi_dir, "arrow_tables", "cell_component_arrows_by_grid.csv"), "LSGI arrow table")
partition_path <- require_file(file.path(lsgi_dir, "grid_partition_membership.csv"), "grid partition membership")
local_path <- require_file(file.path(lsgi_dir, "grid_local_spot_membership.csv"), "grid local spot membership")
grid_info_path <- require_file(file.path(lsgi_dir, "grid_info.csv"), "grid info")
spot_matrix_candidates <- c(
  file.path(paths$output, "07", "spot_matrix_pre_lsgi.csv"),
  file.path(paths$output, "07_spatial_deconvolution", "spot_matrix_pre_lsgi.csv")
)
existing_spot_matrix <- spot_matrix_candidates[file.exists(spot_matrix_candidates)]
spot_matrix_path <- if (length(existing_spot_matrix) > 0) existing_spot_matrix[[1]] else spot_matrix_candidates[[1]]

arrows_all <- readr::read_csv(arrow_path, show_col_types = FALSE)
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

required_arrow_cols <- c("grid", "component", "X", "Y", "vx.u", "vy.u", "included_in_filtered_output")
missing_arrow_cols <- setdiff(required_arrow_cols, colnames(arrows_all))
if (length(missing_arrow_cols) > 0) {
  stop("Missing columns in arrow table: ", paste(missing_arrow_cols, collapse = ", "), call. = FALSE)
}

include_unfiltered_arrows <- as_bool(get_opt(cli_opts, "include-unfiltered-arrows", default = FALSE))
arrow_length_scale <- as.numeric(get_opt(cli_opts, "arrow-length-scale", default = params$lsgi_arrow_length_scale %||% 1.4))
arrow_linewidth <- as.numeric(get_opt(cli_opts, "arrow-linewidth", default = params$lsgi_arrow_linewidth %||% 1.0))
arrow_head_cm <- as.numeric(get_opt(cli_opts, "arrow-head-cm", default = params$lsgi_arrow_head_cm %||% 0.20))
arrow_closed <- as_bool(get_opt(cli_opts, "arrow-closed", default = params$lsgi_arrow_closed %||% TRUE))
plot_nonselected_arrows <- as_bool(get_opt(cli_opts, "plot-nonselected-arrows", default = TRUE))
strategy_filter <- strsplit(as.character(get_opt(cli_opts, "strategy", default = "all")), ",", fixed = TRUE)[[1]]
strategy_filter <- trimws(strategy_filter[nzchar(strategy_filter)])
arrow_type <- if (arrow_closed) "closed" else "open"

candidate_arrows <- arrows_all
if (!include_unfiltered_arrows) {
  candidate_arrows <- candidate_arrows[candidate_arrows$included_in_filtered_output, , drop = FALSE]
}

boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
component_palette <- setNames(
  grDevices::hcl.colors(length(unique(arrows_all$component)), palette = "Dark 3"),
  sort(unique(arrows_all$component))
)

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

passes_partition <- function(metrics, min_n, min_frac) {
  metrics$partition_n_Bdy >= min_n & metrics$partition_frac_Bdy >= min_frac
}

passes_local <- function(metrics, min_n, min_frac) {
  metrics$local_n_Bdy >= min_n & metrics$local_frac_Bdy >= min_frac
}

selected_grid_flags <- function(metrics, strategy) {
  type <- strategy$strategy_type[[1]]
  if (identical(type, "partition")) {
    return(passes_partition(metrics, strategy$min_partition_n_Bdy[[1]], strategy$min_partition_frac_Bdy[[1]]))
  }
  if (identical(type, "local")) {
    return(passes_local(metrics, strategy$min_local_n_Bdy[[1]], strategy$min_local_frac_Bdy[[1]]))
  }
  partition_ok <- passes_partition(metrics, strategy$min_partition_n_Bdy[[1]], strategy$min_partition_frac_Bdy[[1]])
  local_ok <- passes_local(metrics, strategy$min_local_n_Bdy[[1]], strategy$min_local_frac_Bdy[[1]])
  if (identical(type, "intersection")) return(partition_ok & local_ok)
  if (identical(type, "union")) return(partition_ok | local_ok)
  stop("Unsupported strategy_type: ", type, call. = FALSE)
}

make_plot_arrow_df <- function(df) {
  df %>%
    dplyr::mutate(
      X_end = X + vx.u * arrow_length_scale,
      Y_end = Y + vy.u * arrow_length_scale,
      component = factor(component, levels = names(component_palette))
    )
}

add_selected_arrows <- function(p, selected_arrow_df, background_arrow_df = NULL) {
  out <- p
  if (plot_nonselected_arrows && !is.null(background_arrow_df) && nrow(background_arrow_df) > 0) {
    out <- out +
      ggplot2::geom_segment(
        data = background_arrow_df,
        ggplot2::aes(x = X, y = Y, xend = X_end, yend = Y_end),
        inherit.aes = FALSE,
        color = "grey70",
        alpha = 0.25,
        linewidth = max(0.2, arrow_linewidth * 0.55),
        lineend = "round",
        arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm * 0.75, "cm"), type = arrow_type)
      )
  }
  if (nrow(selected_arrow_df) == 0) {
    return(out + ggplot2::labs(subtitle = "No boundary-related LSGI arrows under this strategy"))
  }
  out +
    ggplot2::geom_segment(
      data = selected_arrow_df,
      ggplot2::aes(x = X, y = Y, xend = X_end, yend = Y_end, color = component),
      inherit.aes = FALSE,
      linewidth = arrow_linewidth,
      lineend = "round",
      arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), type = arrow_type)
    ) +
    ggplot2::scale_color_manual(values = component_palette, drop = FALSE) +
    ggplot2::labs(color = "LSGI component")
}

point_df <- partition_df %>%
  dplyr::select(cell_ID, X, Y, Location) %>%
  dplyr::distinct(cell_ID, .keep_all = TRUE) %>%
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))

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
}

background_arrows <- make_plot_arrow_df(candidate_arrows)
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
    dplyr::select(strategy_id, strategy_type, strategy_description, dplyr::everything())

  selected_grids_out <- selected_grids %>%
    dplyr::mutate(
      strategy_id = strategy_id,
      strategy_type = strategy$strategy_type[[1]],
      strategy_description = strategy$description[[1]]
    ) %>%
    dplyr::select(strategy_id, strategy_type, strategy_description, dplyr::everything())

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
    dplyr::inner_join(spot_arrow_join, by = "source_grid", suffix = c("_spot", "_arrow")) %>%
    dplyr::mutate(strategy_id = strategy_id) %>%
    dplyr::select(strategy_id, dplyr::everything())
  if (!is.null(spot_decon_07)) {
    spot_arrow_assignments <- spot_arrow_assignments %>%
      dplyr::left_join(spot_decon_07, by = "cell_ID")
  }
  readr::write_csv(spot_arrow_assignments, file.path(spot_out_dir, paste0(strategy_id, "_bdy_spot_arrow_assignments.csv")))

  selected_plot_arrows <- make_plot_arrow_df(selected_arrows)
  background_plot_arrows <- background_arrows %>%
    dplyr::filter(!grid %in% selected_grids$grid)

  boundary_base <- ggplot2::ggplot(point_df, ggplot2::aes(x = X, y = Y, fill = Location)) +
    ggplot2::geom_point(shape = 21, size = 1.8, stroke = 0.1, color = "grey25", alpha = 0.9) +
    ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::ggtitle(paste0(sample_name, " boundary-related LSGI arrows: ", strategy_id))

  boundary_plot <- add_selected_arrows(boundary_base, selected_plot_arrows, background_plot_arrows)
  ggplot2::ggsave(
    file.path(plot_out_dir, paste0(sample_name, "_", strategy_id, "_BoundaryRelated_LSGIGradient.pdf")),
    boundary_plot,
    width = 8,
    height = 7
  )

  if (has_image) {
    point_df_he <- point_df %>%
      dplyr::mutate(Y = -Y)
    selected_plot_arrows_he <- selected_plot_arrows %>%
      dplyr::mutate(Y = -Y, Y_end = -Y_end)
    background_plot_arrows_he <- background_plot_arrows %>%
      dplyr::mutate(Y = -Y, Y_end = -Y_end)

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
      ggplot2::ggtitle(paste0(sample_name, " HE boundary-related LSGI arrows: ", strategy_id))

    he_plot <- add_selected_arrows(he_base, selected_plot_arrows_he, background_plot_arrows_he)
    ggplot2::ggsave(
      file.path(plot_out_dir, paste0(sample_name, "_", strategy_id, "_BoundaryRelated_HE_LSGIGradient.pdf")),
      he_plot,
      width = 8,
      height = 7
    )
  }

  component_counts <- sort(table(selected_arrows$component), decreasing = TRUE)
  strategy_summary_rows[[i]] <- data.frame(
    strategy_id = strategy_id,
    strategy_type = strategy$strategy_type[[1]],
    min_partition_n_Bdy = strategy$min_partition_n_Bdy[[1]],
    min_partition_frac_Bdy = strategy$min_partition_frac_Bdy[[1]],
    min_local_n_Bdy = strategy$min_local_n_Bdy[[1]],
    min_local_frac_Bdy = strategy$min_local_frac_Bdy[[1]],
    n_selected_grids = nrow(selected_grids),
    n_selected_arrows = nrow(selected_arrows),
    n_selected_bdy_spot_arrow_rows = nrow(spot_arrow_assignments),
    components = paste(names(component_counts), as.integer(component_counts), sep = "=", collapse = ";"),
    description = strategy$description[[1]]
  )
}

strategy_summary <- dplyr::bind_rows(strategy_summary_rows)
readr::write_csv(strategy_summary, file.path(out_dir, "strategy_summary.csv"))
readr::write_csv(grid_metrics, file.path(out_dir, "grid_boundary_metrics.csv"))

sink(file.path(out_dir, "run_summary_12.txt"))
cat("Boundary-related LSGI arrow strategy screen\n")
cat("==========================================\n\n")
cat("Sample:", sample_name, "\n")
cat("Candidate arrows:", nrow(candidate_arrows), "\n")
cat("Include unfiltered arrows:", include_unfiltered_arrows, "\n")
cat("Partition spots:", nrow(partition_df), "\n")
cat("Local membership rows:", nrow(local_df), "\n")
cat("07 spot matrix available:", !is.null(spot_matrix_07), "\n")
cat("07 spot matrix path:", if (file.exists(spot_matrix_path)) spot_matrix_path else "not found", "\n")
cat("H&E overlay available:", has_image, "\n\n")
cat("Strategies:\n")
print(strategy_summary)
cat("\nOutputs written to:\n")
cat(out_dir, "\n")
sink()
