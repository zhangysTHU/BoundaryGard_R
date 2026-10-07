# BoundaryGrad callable numerical core.
#
# This file has no package-loading, file-writing, plotting, or command-line side
# effects. Pipeline scripts 04--06 and the benchmark adapter source this exact
# implementation so equivalence checks exercise the evaluated code path.

boundarygrad_core_api_version <- "1.0.0"

bg_validate_plot_mode <- function(x) {
  if (is.null(x)) x <- "full"
  mode <- tolower(trimws(as.character(x)[[1]]))
  allowed <- c("none", "summary", "full")
  if (!mode %in% allowed) {
    stop("plot_mode must be one of none, summary, or full; got ", sQuote(mode),
         call. = FALSE)
  }
  mode
}

bg_assert_safe_child <- function(path, root, label = "path") {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(root, winslash = "/", mustWork = FALSE)
  prefix <- paste0(sub("/+$", "", root_norm), "/")
  if (identical(path_norm, root_norm) || !startsWith(path_norm, prefix)) {
    stop(label, " must be a child of the configured root: ", path_norm,
         " (root: ", root_norm, ")", call. = FALSE)
  }
  invisible(path_norm)
}

bg_fit_lsgi_local_linear <- function(grids, dist_to_grid, embeddings,
                                     spatial_coords, n_cells_per_meta,
                                     calc_local_linear, optimize_arrow) {
  embeddings <- as.matrix(embeddings)
  storage.mode(embeddings) <- "numeric"
  embeddings <- embeddings[, colSums(abs(embeddings), na.rm = TRUE) > 0,
                           drop = FALSE]
  if (ncol(embeddings) < 1L) {
    stop("No non-zero embedding columns available for LSGI.", call. = FALSE)
  }
  embeddings <- embeddings[rownames(spatial_coords), , drop = FALSE]
  fit <- calc_local_linear(
    grids = grids,
    dist.to.grid = dist_to_grid,
    latent.embeddings = embeddings,
    spatial_coords = spatial_coords,
    n.cells.per.meta = n_cells_per_meta
  )
  grid_info <- cbind(grids, fit$coeff, fit$rsquared.list, fit$assignment.list)
  colnames(grid_info) <- c("X", "Y", "vx", "vy", "R_squared", "Assignment")
  grid_info$Assignment <- factor(grid_info$Assignment, levels = colnames(embeddings))
  grid_info <- optimize_arrow(grid.info = grid_info)
  list(
    grid.info = grid_info,
    local.linear.info = fit,
    grids = grids,
    dist.to.grid = dist_to_grid,
    spatial_coords = spatial_coords,
    embeddings = embeddings
  )
}

bg_passes_boundary_threshold <- function(n_bdy, frac_bdy, min_n, min_frac) {
  n_bdy >= min_n & frac_bdy >= min_frac
}

bg_passes_partition <- function(metrics, min_n, min_frac) {
  bg_passes_boundary_threshold(
    metrics$partition_n_Bdy, metrics$partition_frac_Bdy, min_n, min_frac
  )
}

bg_passes_local <- function(metrics, min_n, min_frac) {
  bg_passes_boundary_threshold(
    metrics$local_n_Bdy, metrics$local_frac_Bdy, min_n, min_frac
  )
}

bg_selected_grid_flags <- function(metrics, strategy) {
  type <- strategy$strategy_type[[1]]
  if (identical(type, "partition")) {
    return(bg_passes_partition(
      metrics,
      strategy$min_partition_n_Bdy[[1]],
      strategy$min_partition_frac_Bdy[[1]]
    ))
  }
  if (identical(type, "local")) {
    return(bg_passes_local(
      metrics,
      strategy$min_local_n_Bdy[[1]],
      strategy$min_local_frac_Bdy[[1]]
    ))
  }
  partition_ok <- bg_passes_partition(
    metrics,
    strategy$min_partition_n_Bdy[[1]],
    strategy$min_partition_frac_Bdy[[1]]
  )
  local_ok <- bg_passes_local(
    metrics,
    strategy$min_local_n_Bdy[[1]],
    strategy$min_local_frac_Bdy[[1]]
  )
  if (identical(type, "intersection")) return(partition_ok & local_ok)
  if (identical(type, "union")) return(partition_ok | local_ok)
  stop("Unsupported strategy_type: ", type, call. = FALSE)
}

bg_boundary_support_score <- function(r_squared, local_frac_bdy) {
  pmax(r_squared, 0) * local_frac_bdy
}

bg_project_axis <- function(x, y, origin_x, origin_y, direction_x, direction_y) {
  direction_norm <- sqrt(direction_x^2 + direction_y^2)
  if (length(direction_norm) != 1L || !is.finite(direction_norm) ||
      direction_norm <= 0) {
    stop("Axis direction must be one finite non-zero 2D vector.", call. = FALSE)
  }
  unit_v <- c(direction_x / direction_norm, direction_y / direction_norm)
  unit_t <- c(-unit_v[[2]], unit_v[[1]])
  rel_x <- as.numeric(x) - origin_x
  rel_y <- as.numeric(y) - origin_y
  list(
    direction_norm = direction_norm,
    unit_v = unit_v,
    unit_t = unit_t,
    axis_s = rel_x * unit_v[[1]] + rel_y * unit_v[[2]],
    lateral_signed = rel_x * unit_t[[1]] + rel_y * unit_t[[2]]
  )
}

bg_select_zero_spot <- function(location, lateral_abs, axis_s, finite_extension,
                                local_rank, width_modes, width_values) {
  if (length(width_modes) != length(width_values)) {
    stop("width_modes and width_values must have equal length.", call. = FALSE)
  }
  zero_idx <- NA_integer_
  zero_mode <- NA_character_
  zero_width <- NA_real_
  zero_candidate_count <- 0L
  for (i in seq_along(width_modes)) {
    candidates <- which(
      location == "Bdy" & lateral_abs <= width_values[[i]] & finite_extension
    )
    zero_candidate_count <- zero_candidate_count + length(candidates)
    if (!length(candidates)) next
    ord <- order(
      lateral_abs[candidates],
      abs(axis_s[candidates]),
      as.numeric(local_rank[candidates])
    )
    zero_idx <- candidates[ord[[1]]]
    zero_mode <- as.character(width_modes[[i]])
    zero_width <- as.numeric(width_values[[i]])
    break
  }
  list(
    found = !is.na(zero_idx),
    index = zero_idx,
    mode = zero_mode,
    width = zero_width,
    candidate_count = zero_candidate_count
  )
}

bg_center_axis_on_zero <- function(axis_s, zero_index) {
  zero_s <- if (!is.na(zero_index)) axis_s[[zero_index]] else 0
  list(zero_s = zero_s, profile_x = axis_s - zero_s)
}
