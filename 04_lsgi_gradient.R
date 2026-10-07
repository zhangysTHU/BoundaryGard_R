# 04: LSGI gradient analysis for multiple embedding/component sources.
# Common grid and spot-membership outputs are written under output/<sample>/11_lsgi_gradient/.
# Method-specific arrows, plots, and summaries are written under method subdirectories, e.g.
# output/<sample>/11_lsgi_gradient/cell_component/ and output/<sample>/11_lsgi_gradient/nmf/.
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "R", "boundarygrad_core.R"))
source(file.path(script_dir, "R", "spatial_plot_core.R"))
load_required_packages(c("Seurat", "readr", "dplyr", "tibble", "ggplot2", "png", "grid", "jsonlite", "patchwork", "viridis", "ComplexHeatmap", "reshape2", "magrittr", "Matrix", "singlet", "msigdbr"))

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

parse_int_vec <- function(x, default) {
  values <- parse_chr_vec(x, as.character(default))
  parsed <- suppressWarnings(as.integer(values))
  if (any(is.na(parsed))) {
    stop("Expected a comma-separated integer vector; got: ", as.character(x), call. = FALSE)
  }
  parsed
}

sanitize_component_names <- function(x, prefix = "component") {
  out <- make.names(as.character(x), unique = TRUE)
  empty <- !nzchar(out) | is.na(out)
  out[empty] <- paste0(prefix, seq_len(sum(empty)))
  out
}

sanitize_filename <- function(x, prefix = "item") {
  out <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  empty <- !nzchar(out) | is.na(out)
  out[empty] <- paste0(prefix, seq_len(sum(empty)))
  make.unique(out, sep = "_")
}

write_matrix_csv <- function(mat, path, row_id = "id") {
  df <- data.frame(mat, check.names = FALSE)
  df <- tibble::rownames_to_column(df, row_id)
  utils::write.csv(df, path, row.names = FALSE)
}

make_cache_key <- function(...) {
  paste(unlist(list(...), use.names = TRUE), collapse = "|")
}

make_matrix_signature <- function(mat, digits = 6) {
  mat <- as.matrix(mat)
  paste(
    nrow(mat),
    ncol(mat),
    paste(rownames(mat) %||% character(), collapse = ";"),
    paste(colnames(mat) %||% character(), collapse = ";"),
    paste(format(round(as.numeric(mat), digits = digits), scientific = FALSE, trim = TRUE), collapse = ","),
    sep = "|"
  )
}

same_numeric_matrix <- function(x, y, tolerance = 1e-8) {
  if (is.null(x) || is.null(y)) {
    return(FALSE)
  }
  x <- as.matrix(x)
  y <- as.matrix(y)
  identical(dim(x), dim(y)) &&
    identical(rownames(x), rownames(y)) &&
    identical(colnames(x), colnames(y)) &&
    isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}

cache_is_compatible <- function(lsgi_res, embedding_result, grid_cache_key, allow_legacy = FALSE) {
  cached_key <- lsgi_res$component_method_cache_key
  cached_grid_key <- lsgi_res$grid_cache_key

  embedding_ok <- if (is.null(cached_key)) {
    isTRUE(allow_legacy) &&
      same_numeric_matrix(lsgi_res$embeddings, embedding_result$embeddings, tolerance = 1e-10)
  } else {
    identical(cached_key, embedding_result$cache_key)
  }
  grid_ok <- if (is.null(cached_grid_key)) {
    identical(nrow(lsgi_res$grid.info), length(grid_ids)) &&
      same_numeric_matrix(lsgi_res$grids, base_grids) &&
      same_numeric_matrix(lsgi_res$spatial_coords, spatial_coords)
  } else {
    identical(cached_grid_key, grid_cache_key)
  }

  embedding_ok && grid_ok
}

get_assay_layer <- function(object, assay, layer) {
  layer_result <- try(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    silent = TRUE
  )
  if (!inherits(layer_result, "try-error")) {
    return(layer_result)
  }
  slot_result <- try(
    Seurat::GetAssayData(object, assay = assay, slot = layer),
    silent = TRUE
  )
  if (!inherits(slot_result, "try-error")) {
    return(slot_result)
  }
  stop("Cannot read assay ", assay, " layer/slot ", layer, " from TumorST.", call. = FALSE)
}

cli_opts <- parse_cli_options(commandArgs(trailingOnly = TRUE))
plot_mode <- bg_validate_plot_mode(get_opt(cli_opts, "plot-mode", default = params$plot_mode %||% "full"))

out_dir <- file.path(paths$output, "11_lsgi_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lsgi_source_override <- Sys.getenv("LSGI_SOURCE_DIR", unset = "")
lsgi_root_candidates <- c(
  if (nzchar(lsgi_source_override)) lsgi_source_override else character(),
  file.path(script_dir, "..", "LSGI-master"),
  file.path(script_dir, "..", "Cottrazm-main", "LSGI-master"),
  file.path(script_dir, "..", "Cottrazm-main-archived", "LSGI-master"),
  file.path(script_dir, "..", "Cottrazm-main-archived", "LSGI_original_work")
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
plot_context <- sp_read_context(paths$spaceranger, load_image = !identical(plot_mode, "none"))

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

scale_factor <- plot_context$tissue_lowres_scalef
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

n_grids_scale <- as.numeric(get_opt(cli_opts, "n-grids-scale", default = params$lsgi_n_grids_scale %||% 10))
n_cells_per_meta <- as.integer(as.numeric(get_opt(cli_opts, "n-cells-per-meta", default = params$lsgi_n_cells_per_meta %||% min(50, nrow(spatial_coords)))))
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
write_legacy_cell_component_outputs <- as_bool(get_opt(cli_opts, "write-legacy-cell-component-outputs", default = TRUE))
plot_common_grid_plots <- identical(plot_mode, "full") &&
  as_bool(get_opt(cli_opts, "plot-common-grid-plots", default = TRUE))
component_methods <- parse_chr_vec(
  get_opt(cli_opts, "component-methods", default = params$lsgi_component_methods %||% "cell_component,nmf,marker_module,pathway,single_gene"),
  default = c("cell_component", "nmf", "marker_module", "pathway", "single_gene")
)

nmf_ranks <- parse_int_vec(get_opt(cli_opts, "nmf-ranks", default = params$lsgi_nmf_ranks %||% "6,7,8,9,10"), default = 6:10)
nmf_cv_replicates <- as.integer(as.numeric(get_opt(cli_opts, "nmf-cv-replicates", default = params$lsgi_nmf_cv_replicates %||% 3)))
nmf_cv_tol <- as.numeric(get_opt(cli_opts, "nmf-cv-tol", default = params$lsgi_nmf_cv_tol %||% 1e-3))
nmf_cv_maxit <- as.integer(as.numeric(get_opt(cli_opts, "nmf-cv-maxit", default = params$lsgi_nmf_cv_maxit %||% 100)))
nmf_final_tol <- as.numeric(get_opt(cli_opts, "nmf-final-tol", default = params$lsgi_nmf_final_tol %||% 1e-4))
nmf_final_maxit <- as.integer(as.numeric(get_opt(cli_opts, "nmf-final-maxit", default = params$lsgi_nmf_final_maxit %||% 300)))
nmf_test_density <- as.numeric(get_opt(cli_opts, "nmf-test-density", default = params$lsgi_nmf_test_density %||% 0.05))
nmf_l1 <- as.numeric(get_opt(cli_opts, "nmf-l1", default = params$lsgi_nmf_l1 %||% 0.01))
nmf_l2 <- as.numeric(get_opt(cli_opts, "nmf-l2", default = params$lsgi_nmf_l2 %||% 0))
nmf_threads <- as.integer(as.numeric(get_opt(cli_opts, "nmf-threads", default = params$lsgi_nmf_threads %||% 1)))
nmf_precision <- as.character(get_opt(cli_opts, "nmf-precision", default = params$lsgi_nmf_precision %||% "double"))
nmf_assay <- as.character(get_opt(cli_opts, "nmf-assay", default = params$lsgi_nmf_assay %||% "Spatial"))
nmf_layer <- as.character(get_opt(cli_opts, "nmf-layer", default = params$lsgi_nmf_layer %||% "counts"))
nmf_top_genes <- as.integer(as.numeric(get_opt(cli_opts, "nmf-top-genes", default = params$lsgi_nmf_top_genes %||% 2000)))
nmf_min_gene_spots <- as.integer(as.numeric(get_opt(cli_opts, "nmf-min-gene-spots", default = params$lsgi_nmf_min_gene_spots %||% 10)))
nmf_scale_factor <- as.numeric(get_opt(cli_opts, "nmf-scale-factor", default = params$lsgi_nmf_scale_factor %||% 10000))
nmf_rank_error_tolerance <- as.numeric(get_opt(cli_opts, "nmf-rank-error-tolerance", default = params$lsgi_nmf_rank_error_tolerance %||% 0.01))
nmf_seed <- as.integer(as.numeric(get_opt(cli_opts, "nmf-seed", default = params$lsgi_nmf_seed %||% 666)))

embedding_catalog_dir <- as.character(get_opt(
  cli_opts,
  "embedding-catalog-dir",
  default = params$lsgi_embedding_catalog_dir %||% file.path(paths$resources, "lsgi_embedding_catalog")
))
lsgi_expression_assay <- as.character(get_opt(cli_opts, "expression-assay", default = params$lsgi_expression_assay %||% "Spatial"))
lsgi_expression_layer <- as.character(get_opt(cli_opts, "expression-layer", default = params$lsgi_expression_layer %||% "counts"))
lsgi_expression_scale_factor <- as.numeric(get_opt(
  cli_opts,
  "expression-scale-factor",
  default = params$lsgi_expression_scale_factor %||% 10000
))
lsgi_min_feature_genes <- as.integer(as.numeric(get_opt(
  cli_opts,
  "min-feature-genes",
  default = params$lsgi_min_feature_genes %||% 2
)))
marker_module_ids <- parse_chr_vec(
  get_opt(cli_opts, "marker-module-ids", default = params$lsgi_marker_module_ids %||% character()),
  default = character()
)
pathway_ids <- parse_chr_vec(
  get_opt(cli_opts, "pathway-ids", default = params$lsgi_pathway_ids %||% character()),
  default = character()
)
single_gene_ids <- parse_chr_vec(
  get_opt(cli_opts, "single-gene-ids", default = params$lsgi_single_gene_ids %||% character()),
  default = character()
)

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
if (n_cells_per_meta < 3 || n_cells_per_meta > nrow(spatial_coords)) {
  stop("--n-cells-per-meta must be between 3 and the number of matched spots.", call. = FALSE)
}

allowed_methods <- c("cell_component", "nmf", "marker_module", "pathway", "single_gene")
unknown_methods <- setdiff(component_methods, allowed_methods)
if (length(unknown_methods) > 0) {
  stop("Unknown component method(s): ", paste(unknown_methods, collapse = ", "), call. = FALSE)
}
catalog_methods <- c("marker_module", "pathway", "single_gene")
if (any(component_methods %in% catalog_methods) && !dir.exists(embedding_catalog_dir)) {
  stop("Cannot find LSGI embedding catalog directory: ", embedding_catalog_dir, call. = FALSE)
}
if (!is.finite(lsgi_expression_scale_factor) || lsgi_expression_scale_factor <= 0) {
  stop("--expression-scale-factor must be a finite positive value.", call. = FALSE)
}
if (is.na(lsgi_min_feature_genes) || lsgi_min_feature_genes < 1) {
  stop("--min-feature-genes must be an integer >= 1.", call. = FALSE)
}

spot_base_df <- spatial_df %>%
  dplyr::left_join(DeconData, by = "cell_ID")

boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
point_df <- spot_base_df %>%
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))

has_image <- plot_context$has_image
point_spot_polygons <- if (!identical(plot_mode, "none")) {
  sp_make_spot_polygons(
    point_df,
    plot_context,
    segments = params$spatial_plot_spot_segments %||% 24L
  )
} else NULL

message("Building common LSGI grid from spatial coordinates.")
base_grids <- get.grid.coords(spatial_coords = spatial_coords, n.grids.scale = n_grids_scale)
colnames(base_grids) <- c("X", "Y")
grid_ids <- paste0("grid_", seq_len(nrow(base_grids)))
base_dist_to_grid <- spa_vectorized_pdist(A = as.matrix(spatial_coords), B = as.matrix(base_grids))
base_dist_to_grid <- as.matrix(base_dist_to_grid)
rownames(base_dist_to_grid) <- rownames(spatial_coords)
colnames(base_dist_to_grid) <- grid_ids
grid_cache_key <- make_cache_key(
  n_grids_scale = n_grids_scale,
  n_cells_per_meta = n_cells_per_meta,
  spatial_coords = make_matrix_signature(spatial_coords),
  grids = make_matrix_signature(base_grids)
)

grid_info_common <- base_grids %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    grid = grid_ids,
    grid_index = seq_len(dplyr::n())
  ) %>%
  dplyr::select(grid, grid_index, X, Y)
utils::write.csv(grid_info_common, file.path(out_dir, "grid_info.csv"), row.names = FALSE)

grid_palette <- setNames(
  grDevices::hcl.colors(length(grid_ids), palette = "Dynamic"),
  grid_ids
)

grid_centers_df <- grid_info_common %>%
  dplyr::transmute(
    grid,
    grid_center_X = X,
    grid_center_Y = Y
  )

local_membership_list <- vector("list", length(grid_ids))
for (i in seq_along(grid_ids)) {
  cells <- rownames(base_dist_to_grid)[order(base_dist_to_grid[, i])[seq_len(n_cells_per_meta)]]
  local_membership_list[[i]] <- tibble::tibble(
    grid = grid_ids[i],
    cell_ID = cells,
    local_rank = seq_along(cells),
    dist_to_grid = as.numeric(base_dist_to_grid[cells, i])
  )
}

local_membership_df <- dplyr::bind_rows(local_membership_list) %>%
  dplyr::left_join(spot_base_df[, c("cell_ID", "X", "Y", "Location"), drop = FALSE], by = "cell_ID") %>%
  dplyr::left_join(grid_centers_df, by = "grid") %>%
  dplyr::arrange(grid, local_rank)
utils::write.csv(local_membership_df, file.path(out_dir, "grid_local_spot_membership.csv"), row.names = FALSE)

nearest_idx <- max.col(-base_dist_to_grid, ties.method = "first")
partition_df <- tibble::tibble(
  cell_ID = rownames(base_dist_to_grid),
  nearest_grid = grid_ids[nearest_idx],
  nearest_grid_distance = base_dist_to_grid[cbind(seq_len(nrow(base_dist_to_grid)), nearest_idx)]
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

grid_table_df <- grid_info_common %>%
  dplyr::left_join(grid_local_summary, by = "grid") %>%
  dplyr::left_join(grid_partition_summary, by = "grid") %>%
  dplyr::mutate(
    n_local_spots = ifelse(is.na(n_local_spots), 0L, n_local_spots),
    local_spots = dplyr::coalesce(local_spots, ""),
    n_partition_spots = ifelse(is.na(n_partition_spots), 0L, n_partition_spots),
    partition_spots = dplyr::coalesce(partition_spots, "")
  )
utils::write.csv(grid_table_df, file.path(out_dir, "grid_spot_summary.csv"), row.names = FALSE)

make_lsgi_result <- function(embeddings) {
  bg_fit_lsgi_local_linear(
    grids = base_grids,
    dist_to_grid = base_dist_to_grid,
    embeddings = embeddings,
    spatial_coords = spatial_coords,
    n_cells_per_meta = n_cells_per_meta,
    calc_local_linear = calc.local.linear,
    optimize_arrow = optimize.arrow
  )
}

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

calc_gradient_distance_fast <- function(lin_res) {
  distance_input <- lin_res %>%
    dplyr::filter(rsquared > r_squared_thresh) %>%
    dplyr::group_by(fctr) %>%
    dplyr::filter(dplyr::n() >= minimum_fctr) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fctr = as.character(fctr))
  if (nrow(distance_input) < 2) {
    return(NULL)
  }
  factors <- unique(distance_input$fctr)
  rows <- vector("list", length(factors) * length(factors))
  row_i <- 0L
  for (source_factor in factors) {
    source_df <- distance_input[distance_input$fctr == source_factor, , drop = FALSE]
    source_coords <- as.matrix(source_df[, c("X", "Y"), drop = FALSE])
    for (target_factor in factors) {
      target_df <- distance_input[distance_input$fctr == target_factor, , drop = FALSE]
      target_coords <- as.matrix(target_df[, c("X", "Y"), drop = FALSE])
      pair_dist <- spa_vectorized_pdist(A = source_coords, B = target_coords)
      if (identical(source_factor, target_factor)) {
        diag(pair_dist) <- NA_real_
      }
      nearest <- apply(pair_dist, 1, function(x) {
        finite_x <- x[is.finite(x)]
        if (length(finite_x) == 0) NA_real_ else min(finite_x)
      })
      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        grid.1 = source_factor,
        grid.2 = target_factor,
        distance = mean(nearest, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  dist_df <- dplyr::bind_rows(rows[seq_len(row_i)])
  dist_df <- dist_df[order_by_numeric_suffix(dist_df$grid.1), , drop = FALSE]
  dist_df
}

write_distance_heatmap <- function(dist_mat, output_path, label) {
  if (is.null(dist_mat) || nrow(dist_mat) == 0 || !"distance" %in% colnames(dist_mat)) {
    return(FALSE)
  }
  heatmap_df <- dist_mat[is.finite(dist_mat$distance) & dist_mat$distance > 0, , drop = FALSE]
  if (nrow(heatmap_df) == 0) {
    warning("LSGI distance heatmap skipped for ", label, ": no finite positive distances.", call. = FALSE)
    return(FALSE)
  }
  ok <- tryCatch(
    {
      mat2plot <- reshape2::acast(heatmap_df, formula = grid.1 ~ grid.2, value.var = "distance", fun.aggregate = mean)
      mat2plot <- log2(mat2plot)
      mat2plot[!is.finite(mat2plot)] <- NA_real_
      if (!any(is.finite(mat2plot))) {
        stop("all log2 distances are non-finite")
      }
      row_idx <- order_by_numeric_suffix(rownames(mat2plot))
      col_idx <- order_by_numeric_suffix(colnames(mat2plot))
      mat2plot <- mat2plot[row_idx, col_idx, drop = FALSE]
      grDevices::pdf(output_path, width = 7, height = 6)
      ComplexHeatmap::draw(
        ComplexHeatmap::Heatmap(
          mat2plot,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = rev(viridis::viridis(30)),
          heatmap_legend_param = list(title = "log2(distance)")
        )
      )
      grDevices::dev.off()
      TRUE
    },
    error = function(e) {
      if (grDevices::dev.cur() > 1) {
        grDevices::dev.off()
      }
      warning("LSGI distance heatmap skipped for ", label, ": ", conditionMessage(e), call. = FALSE)
      FALSE
    }
  )
  isTRUE(ok)
}

plot_grid_memberships <- function() {
  partition_polygons <- sp_make_spot_polygons(
    partition_df,
    plot_context,
    segments = params$spatial_plot_spot_segments %||% 24L
  )
  grid_centers_plot <- grid_info_common %>% dplyr::mutate(Y = -Y)

  make_partition_plot <- function(show_image, title_text) {
    sp_spatial_canvas(plot_context, show_image = show_image, title = title_text) +
      sp_spot_layer(
        partition_polygons,
        fill_col = "nearest_grid",
        colour = NA,
        alpha = 0.9,
        linewidth = 0
      ) +
    ggplot2::geom_point(
      data = grid_centers_plot,
      ggplot2::aes(x = X, y = Y),
      inherit.aes = FALSE,
      shape = 4,
      size = 1.8,
      stroke = 0.8,
      color = "black"
    ) +
      ggplot2::scale_fill_manual(values = grid_palette, guide = "none")
  }

  partition_plot <- make_partition_plot(FALSE, paste0(sample_name, " grid partition (nearest-grid assignment)"))
  partition_overlay_plot <- make_partition_plot(has_image, paste0(sample_name, " grid partition overlay"))
  sp_save_plot(
    partition_plot,
    file.path(out_dir, "grid_partition.pdf"),
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    map_width = params$spatial_plot_map_width %||% 6.45,
    legend_width = params$spatial_plot_legend_width %||% 1.55
  )
  sp_save_plot(
    partition_overlay_plot,
    file.path(out_dir, "grid_partition_overlay.pdf"),
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    map_width = params$spatial_plot_map_width %||% 6.45,
    legend_width = params$spatial_plot_legend_width %||% 1.55
  )

  make_local_plot <- function(grid_id, overlay = FALSE) {
    local_df <- local_membership_df[local_membership_df$grid == grid_id, , drop = FALSE]
    center_df <- grid_info_common[grid_info_common$grid == grid_id, , drop = FALSE]
    center_df$Y <- -center_df$Y
    grid_color <- grid_palette[[grid_id]]
    title_text <- paste0(sample_name, " ", grid_id, " local spots (n=", nrow(local_df), ")")
    local_polygons <- sp_make_spot_polygons(
      local_df,
      plot_context,
      segments = params$spatial_plot_spot_segments %||% 24L
    )

    sp_spatial_canvas(
      plot_context,
      show_image = isTRUE(overlay) && has_image,
      title = if (isTRUE(overlay)) paste0(title_text, " overlay") else title_text
    ) +
      sp_spot_layer(
        point_spot_polygons,
        fill = "grey82",
        colour = NA,
        alpha = if (isTRUE(overlay)) 0.55 else 0.65,
        linewidth = 0
      ) +
      sp_spot_layer(
        local_polygons,
        fill = grid_color,
        colour = "grey20",
        alpha = 0.95,
        linewidth = 0.1
      ) +
      ggplot2::geom_point(
        data = center_df,
        ggplot2::aes(x = X, y = Y),
        inherit.aes = FALSE,
        shape = 4,
        size = 2.3,
        stroke = 0.9,
        color = "black"
      )
  }

  grDevices::cairo_pdf(
    file.path(out_dir, "grid_local_spots.pdf"),
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    onefile = TRUE
  )
  for (grid_id in grid_ids) {
    print(sp_fixed_layout(
      make_local_plot(grid_id, overlay = FALSE),
      map_width = params$spatial_plot_map_width %||% 6.45,
      legend_width = params$spatial_plot_legend_width %||% 1.55
    ))
  }
  grDevices::dev.off()

  grDevices::cairo_pdf(
    file.path(out_dir, "grid_local_spots_overlay.pdf"),
    width = params$spatial_plot_width %||% 8,
    height = params$spatial_plot_height %||% 7,
    onefile = TRUE
  )
  for (grid_id in grid_ids) {
    print(sp_fixed_layout(
      make_local_plot(grid_id, overlay = TRUE),
      map_width = params$spatial_plot_map_width %||% 6.45,
      legend_width = params$spatial_plot_legend_width %||% 1.55
    ))
  }
  grDevices::dev.off()
}

collapse_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "" else paste(unique(x), collapse = ";")
}

split_catalog_genes <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(character())
  }
  genes <- unlist(strsplit(as.character(x[[1]]), "[,;]", perl = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  unique(genes[nzchar(genes)])
}

read_lsgi_catalog <- function(file_name, required_cols) {
  path <- file.path(embedding_catalog_dir, file_name)
  if (!file.exists(path)) {
    stop("Cannot find LSGI embedding catalog file: ", path, call. = FALSE)
  }
  catalog <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  missing_cols <- setdiff(required_cols, colnames(catalog))
  if (length(missing_cols) > 0) {
    stop("Catalog ", path, " is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  attr(catalog, "catalog_path") <- normalizePath(path, winslash = "/", mustWork = FALSE)
  catalog
}

select_catalog_rows <- function(catalog, id_col, selected_ids, catalog_label) {
  catalog[[id_col]] <- trimws(as.character(catalog[[id_col]]))
  if (any(!nzchar(catalog[[id_col]]) | is.na(catalog[[id_col]]))) {
    stop(catalog_label, " catalog contains empty IDs in column ", id_col, ".", call. = FALSE)
  }
  duplicated_ids <- unique(catalog[[id_col]][duplicated(catalog[[id_col]])])
  if (length(duplicated_ids) > 0) {
    stop(catalog_label, " catalog contains duplicated ID(s): ", paste(duplicated_ids, collapse = ", "), call. = FALSE)
  }
  selected_ids <- unique(trimws(as.character(selected_ids)))
  selected_ids <- selected_ids[nzchar(selected_ids)]
  if (length(selected_ids) == 0) {
    return(catalog)
  }
  missing_ids <- setdiff(selected_ids, catalog[[id_col]])
  if (length(missing_ids) > 0) {
    stop(
      "Unknown ", catalog_label, " ID(s): ", paste(missing_ids, collapse = ", "),
      ". Edit 00_config.R or ", attr(catalog, "catalog_path"), ".",
      call. = FALSE
    )
  }
  catalog[match(selected_ids, catalog[[id_col]]), , drop = FALSE]
}

match_expression_genes <- function(genes, expression_genes) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  if (length(genes) == 0) {
    return(data.frame(
      requested_gene = character(),
      matched_gene = character(),
      gene_index = integer(),
      present = logical(),
      stringsAsFactors = FALSE
    ))
  }
  idx <- match(genes, expression_genes)
  missing_exact <- is.na(idx)
  if (any(missing_exact)) {
    expression_upper <- toupper(expression_genes)
    first_upper <- !duplicated(expression_upper)
    upper_to_idx <- stats::setNames(which(first_upper), expression_upper[first_upper])
    idx[missing_exact] <- unname(upper_to_idx[toupper(genes[missing_exact])])
  }
  present <- !is.na(idx)
  data.frame(
    requested_gene = genes,
    matched_gene = ifelse(present, expression_genes[idx], NA_character_),
    gene_index = as.integer(idx),
    present = present,
    stringsAsFactors = FALSE
  )
}

lsgi_expression_matrix_cache <- NULL
get_lsgi_expression_matrix <- function() {
  if (!is.null(lsgi_expression_matrix_cache)) {
    return(lsgi_expression_matrix_cache)
  }
  spot_names <- rownames(spatial_coords)
  counts <- get_assay_layer(TumorST, assay = lsgi_expression_assay, layer = lsgi_expression_layer)
  missing_spots <- setdiff(spot_names, colnames(counts))
  if (length(missing_spots) > 0) {
    stop(
      "Expression input assay is missing matched spots: ",
      paste(head(missing_spots, 10), collapse = ", "),
      call. = FALSE
    )
  }
  counts <- counts[, spot_names, drop = FALSE]
  counts <- as(as(as(counts, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  rownames(counts) <- rownames(get_assay_layer(TumorST, assay = lsgi_expression_assay, layer = lsgi_expression_layer))
  colnames(counts) <- spot_names
  counts@x[counts@x < 0] <- 0
  lib_size <- Matrix::colSums(counts)
  if (any(lib_size <= 0)) {
    stop("Expression input contains spots with zero library size.", call. = FALSE)
  }
  norm_mat <- counts %*% Matrix::Diagonal(x = lsgi_expression_scale_factor / lib_size)
  norm_mat <- as(norm_mat, "CsparseMatrix")
  norm_mat@x <- log1p(norm_mat@x)
  rownames(norm_mat) <- rownames(counts)
  colnames(norm_mat) <- spot_names
  lsgi_expression_matrix_cache <<- norm_mat
  norm_mat
}

zscore_selected_genes <- function(norm_mat, gene_indices) {
  gene_indices <- unique(as.integer(gene_indices[!is.na(gene_indices)]))
  selected <- as.matrix(norm_mat[gene_indices, , drop = FALSE])
  gene_means <- rowMeans(selected)
  gene_sds <- sqrt(pmax(rowMeans(selected^2) - gene_means^2, 0))
  gene_sds[!is.finite(gene_sds) | gene_sds <= 0] <- 1
  z <- sweep(selected, 1, gene_means, "-")
  z <- sweep(z, 1, gene_sds, "/")
  rownames(z) <- rownames(norm_mat)[gene_indices]
  colnames(z) <- colnames(norm_mat)
  z
}

sparse_rank_percentiles <- function(norm_mat, gene_indices) {
  norm_mat <- as(norm_mat, "CsparseMatrix")
  gene_indices <- unique(as.integer(gene_indices[!is.na(gene_indices)]))
  n_genes <- nrow(norm_mat)
  out <- matrix(
    NA_real_,
    nrow = length(gene_indices),
    ncol = ncol(norm_mat),
    dimnames = list(rownames(norm_mat)[gene_indices], colnames(norm_mat))
  )
  if (length(gene_indices) == 0) {
    return(out)
  }

  for (j in seq_len(ncol(norm_mat))) {
    start <- norm_mat@p[j] + 1L
    end <- norm_mat@p[j + 1L]
    nnz <- end - start + 1L
    if (end < start) {
      nnz <- 0L
    }
    zero_count <- n_genes - nnz
    if (zero_count > 0) {
      out[, j] <- ((1 + zero_count) / 2) / n_genes
    }
    if (nnz == 0L) {
      next
    }
    idx <- norm_mat@i[start:end] + 1L
    vals <- norm_mat@x[start:end]
    ord <- order(vals)
    sorted_vals <- vals[ord]
    sorted_idx <- idx[ord]
    value_runs <- rle(sorted_vals)
    run_ends <- cumsum(value_runs$lengths)
    run_starts <- run_ends - value_runs$lengths + 1L
    run_avg_ranks <- zero_count + (run_starts + run_ends) / 2
    sorted_ranks <- rep(run_avg_ranks, value_runs$lengths)
    selected_pos <- match(gene_indices, sorted_idx)
    present <- !is.na(selected_pos)
    out[present, j] <- sorted_ranks[selected_pos[present]] / n_genes
  }
  out
}

score_gene_set_embeddings <- function(feature_defs, method_dir, method_id, source_label, catalog_path = NULL) {
  required_cols <- c("feature_id", "label", "scorer", "genes_up", "genes_down")
  missing_cols <- setdiff(required_cols, colnames(feature_defs))
  if (length(missing_cols) > 0) {
    stop("Internal feature definition is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  norm_mat <- get_lsgi_expression_matrix()
  expression_genes <- rownames(norm_mat)
  supported_scorers <- c("zscore_mean", "zscore_contrast", "rank_mean", "rank_contrast")
  feature_defs$feature_id <- trimws(as.character(feature_defs$feature_id))
  feature_defs$label <- trimws(as.character(feature_defs$label))
  feature_defs$scorer <- tolower(trimws(as.character(feature_defs$scorer)))
  unsupported_scorers <- setdiff(unique(feature_defs$scorer), supported_scorers)
  if (length(unsupported_scorers) > 0) {
    stop(
      "Unsupported scorer(s) for ", method_id, ": ",
      paste(unsupported_scorers, collapse = ", "),
      ". Supported scorers: ", paste(supported_scorers, collapse = ", "),
      call. = FALSE
    )
  }

  up_matches <- vector("list", nrow(feature_defs))
  down_matches <- vector("list", nrow(feature_defs))
  up_indices <- vector("list", nrow(feature_defs))
  down_indices <- vector("list", nrow(feature_defs))
  missing_rows <- list()
  for (i in seq_len(nrow(feature_defs))) {
    up_matches[[i]] <- match_expression_genes(split_catalog_genes(feature_defs$genes_up[i]), expression_genes)
    down_matches[[i]] <- match_expression_genes(split_catalog_genes(feature_defs$genes_down[i]), expression_genes)
    up_indices[[i]] <- up_matches[[i]]$gene_index[up_matches[[i]]$present]
    down_indices[[i]] <- down_matches[[i]]$gene_index[down_matches[[i]]$present]
    if (any(!up_matches[[i]]$present)) {
      missing_rows[[length(missing_rows) + 1L]] <- data.frame(
        feature_id = feature_defs$feature_id[i],
        direction = "up",
        requested_gene = up_matches[[i]]$requested_gene[!up_matches[[i]]$present],
        stringsAsFactors = FALSE
      )
    }
    if (any(!down_matches[[i]]$present)) {
      missing_rows[[length(missing_rows) + 1L]] <- data.frame(
        feature_id = feature_defs$feature_id[i],
        direction = "down",
        requested_gene = down_matches[[i]]$requested_gene[!down_matches[[i]]$present],
        stringsAsFactors = FALSE
      )
    }
  }

  uses_contrast <- grepl("_contrast$", feature_defs$scorer)
  usable <- vapply(up_indices, length, integer(1)) >= lsgi_min_feature_genes
  usable <- usable & (!uses_contrast | vapply(down_indices, length, integer(1)) >= 1L)
  component_names <- sanitize_component_names(feature_defs$feature_id, prefix = paste0(method_id, "_"))
  skip_reason <- ifelse(
    usable,
    "",
    ifelse(
      vapply(up_indices, length, integer(1)) < lsgi_min_feature_genes,
      paste0("fewer_than_", lsgi_min_feature_genes, "_up_genes_present"),
      "no_down_genes_present_for_contrast"
    )
  )
  feature_summary <- data.frame(
    feature_id = feature_defs$feature_id,
    component = component_names,
    label = feature_defs$label,
    scorer = feature_defs$scorer,
    n_up_requested = vapply(up_matches, nrow, integer(1)),
    n_up_present = vapply(up_indices, length, integer(1)),
    n_down_requested = vapply(down_matches, nrow, integer(1)),
    n_down_present = vapply(down_indices, length, integer(1)),
    genes_up_present = vapply(up_matches, function(x) collapse_chr(x$matched_gene[x$present]), character(1)),
    genes_up_missing = vapply(up_matches, function(x) collapse_chr(x$requested_gene[!x$present]), character(1)),
    genes_down_present = vapply(down_matches, function(x) collapse_chr(x$matched_gene[x$present]), character(1)),
    genes_down_missing = vapply(down_matches, function(x) collapse_chr(x$requested_gene[!x$present]), character(1)),
    used_in_lsgi = usable,
    skip_reason = skip_reason,
    stringsAsFactors = FALSE
  )
  extra_cols <- setdiff(colnames(feature_defs), c("feature_id", "label", "scorer", "genes_up", "genes_down"))
  if (length(extra_cols) > 0) {
    feature_summary <- cbind(feature_summary, feature_defs[, extra_cols, drop = FALSE])
  }
  utils::write.csv(feature_summary, file.path(method_dir, "selected_features.csv"), row.names = FALSE)

  missing_df <- if (length(missing_rows) == 0) {
    data.frame(feature_id = character(), direction = character(), requested_gene = character(), stringsAsFactors = FALSE)
  } else {
    dplyr::bind_rows(missing_rows)
  }
  utils::write.csv(missing_df, file.path(method_dir, "missing_features.csv"), row.names = FALSE)

  if (!any(usable)) {
    stop("No usable ", method_id, " features after matching genes to expression matrix.", call. = FALSE)
  }

  score_mat <- matrix(
    NA_real_,
    nrow = ncol(norm_mat),
    ncol = sum(usable),
    dimnames = list(colnames(norm_mat), component_names[usable])
  )

  z_feature_idx <- which(usable & grepl("^zscore_", feature_defs$scorer))
  rank_feature_idx <- which(usable & grepl("^rank_", feature_defs$scorer))
  z_scores <- NULL
  rank_scores <- NULL
  if (length(z_feature_idx) > 0) {
    z_gene_indices <- unique(unlist(c(up_indices[z_feature_idx], down_indices[z_feature_idx]), use.names = FALSE))
    z_scores <- zscore_selected_genes(norm_mat, z_gene_indices)
  }
  if (length(rank_feature_idx) > 0) {
    rank_gene_indices <- unique(unlist(c(up_indices[rank_feature_idx], down_indices[rank_feature_idx]), use.names = FALSE))
    message("Computing sparse rank-percentile scores for ", length(rank_gene_indices), " genes across ", ncol(norm_mat), " spots.")
    rank_scores <- sparse_rank_percentiles(norm_mat, rank_gene_indices)
  }

  for (i in which(usable)) {
    source_scores <- if (grepl("^zscore_", feature_defs$scorer[i])) z_scores else rank_scores
    up_gene_names <- rownames(norm_mat)[up_indices[[i]]]
    score <- colMeans(source_scores[up_gene_names, , drop = FALSE])
    if (uses_contrast[i]) {
      down_gene_names <- rownames(norm_mat)[down_indices[[i]]]
      score <- score - colMeans(source_scores[down_gene_names, , drop = FALSE])
    }
    score_mat[, component_names[i]] <- as.numeric(score)
  }
  if (any(!is.finite(score_mat))) {
    stop(method_id, " scorer produced non-finite embedding values.", call. = FALSE)
  }
  nonzero_score <- colSums(abs(score_mat), na.rm = TRUE) > 0
  if (any(!nonzero_score)) {
    zero_components <- colnames(score_mat)[!nonzero_score]
    zero_idx <- match(zero_components, feature_summary$component)
    feature_summary$used_in_lsgi[zero_idx] <- FALSE
    feature_summary$skip_reason[zero_idx] <- "zero_variance_or_zero_score"
    score_mat <- score_mat[, nonzero_score, drop = FALSE]
  }
  if (ncol(score_mat) == 0) {
    stop("No non-zero ", method_id, " embeddings available for LSGI.", call. = FALSE)
  }
  utils::write.csv(feature_summary, file.path(method_dir, "selected_features.csv"), row.names = FALSE)
  write_matrix_csv(score_mat, file.path(method_dir, "embedding_matrix.csv"), row_id = "cell_ID")
  utils::write.csv(
    data.frame(
      method_id = method_id,
      assay = lsgi_expression_assay,
      layer = lsgi_expression_layer,
      scale_factor = lsgi_expression_scale_factor,
      min_feature_genes = lsgi_min_feature_genes,
      catalog_path = catalog_path %||% "",
      n_requested_features = nrow(feature_defs),
      n_used_features = ncol(score_mat),
      scorers = paste(unique(feature_defs$scorer[usable]), collapse = ","),
      stringsAsFactors = FALSE
    ),
    file.path(method_dir, "score_parameters.csv"),
    row.names = FALSE
  )

  list(
    embeddings = score_mat,
    metadata = list(
      source = source_label,
      selected_ids = feature_summary$feature_id[feature_summary$used_in_lsgi],
      selected_features = ncol(score_mat),
      scorer = paste(unique(feature_summary$scorer[feature_summary$used_in_lsgi]), collapse = ", "),
      catalog_path = catalog_path
    ),
    cache_key = make_cache_key(
      method = method_id,
      source = source_label,
      assay = lsgi_expression_assay,
      layer = lsgi_expression_layer,
      scale_factor = lsgi_expression_scale_factor,
      min_feature_genes = lsgi_min_feature_genes,
      features = paste(
        feature_summary$feature_id,
        feature_summary$scorer,
        feature_summary$genes_up_present,
        feature_summary$genes_down_present,
        feature_summary$used_in_lsgi,
        sep = ":",
        collapse = "|"
      )
    )
  )
}

make_cell_component_embeddings <- function() {
  embeddings <- as.matrix(DeconData[, decon_cols, drop = FALSE])
  storage.mode(embeddings) <- "numeric"
  rownames(embeddings) <- DeconData$cell_ID
  embeddings <- embeddings[, colSums(abs(embeddings), na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(embeddings) < 1) {
    stop("No non-zero cell-component columns found in 07_DeconData.rds.gz.", call. = FALSE)
  }
  colnames(embeddings) <- sanitize_component_names(colnames(embeddings), prefix = "cell_component_")
  list(
    embeddings = embeddings,
    metadata = list(source = "07_DeconData.rds.gz"),
    cache_key = make_cache_key(
      method = "cell_component",
      n_spots = nrow(embeddings),
      n_components = ncol(embeddings),
      components = paste(colnames(embeddings), collapse = ";")
    )
  )
}

make_nmf_embeddings <- function(method_dir) {
  spot_names <- rownames(spatial_coords)
  counts <- get_assay_layer(TumorST, assay = nmf_assay, layer = nmf_layer)
  missing_spots <- setdiff(spot_names, colnames(counts))
  if (length(missing_spots) > 0) {
    stop("NMF input assay is missing matched spots: ", paste(head(missing_spots, 10), collapse = ", "), call. = FALSE)
  }
  counts <- counts[, spot_names, drop = FALSE]
  counts <- as(as(as(counts, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  colnames(counts) <- spot_names
  counts@x[counts@x < 0] <- 0

  detected_spots <- Matrix::rowSums(counts > 0)
  total_counts <- Matrix::rowSums(counts)
  keep_genes <- detected_spots >= nmf_min_gene_spots & total_counts > 0
  counts <- counts[keep_genes, , drop = FALSE]
  colnames(counts) <- spot_names
  if (nrow(counts) < max(nmf_ranks)) {
    stop("Too few genes remain for NMF after filtering.", call. = FALSE)
  }

  lib_size <- Matrix::colSums(counts)
  if (any(lib_size <= 0)) {
    stop("NMF input contains spots with zero library size after gene filtering.", call. = FALSE)
  }
  norm_mat <- counts %*% Matrix::Diagonal(x = nmf_scale_factor / lib_size)
  colnames(norm_mat) <- spot_names
  rownames(norm_mat) <- rownames(counts)
  norm_mat@x <- log1p(norm_mat@x)

  gene_means <- Matrix::rowMeans(norm_mat)
  gene_vars <- Matrix::rowMeans(norm_mat^2) - gene_means^2
  gene_vars[!is.finite(gene_vars)] <- NA_real_
  ordered_genes <- names(sort(gene_vars, decreasing = TRUE, na.last = NA))
  selected_genes <- head(ordered_genes, min(nmf_top_genes, length(ordered_genes)))
  if (length(selected_genes) < max(nmf_ranks)) {
    stop("Too few variable genes available for NMF.", call. = FALSE)
  }
  A <- norm_mat[selected_genes, , drop = FALSE]
  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  rownames(A) <- selected_genes
  colnames(A) <- spot_names

  nmf_input_summary <- data.frame(
    assay = nmf_assay,
    layer = nmf_layer,
    input_genes_before_filter = nrow(get_assay_layer(TumorST, assay = nmf_assay, layer = nmf_layer)),
    input_spots = ncol(A),
    genes_after_min_spot_filter = nrow(counts),
    selected_variable_genes = nrow(A),
    min_gene_spots = nmf_min_gene_spots,
    scale_factor = nmf_scale_factor
  )
  utils::write.csv(nmf_input_summary, file.path(method_dir, "nmf_input_summary.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(
      gene = selected_genes,
      variance = gene_vars[selected_genes],
      detected_spots = detected_spots[selected_genes],
      total_counts = total_counts[selected_genes],
      row.names = NULL
    ),
    file.path(method_dir, "nmf_selected_genes.csv"),
    row.names = FALSE
  )

  set.seed(nmf_seed)
  message("Running singlet NMF cross-validation for ranks: ", paste(nmf_ranks, collapse = ", "))
  cv_data <- singlet::cross_validate_nmf(
    A,
    ranks = nmf_ranks,
    n_replicates = nmf_cv_replicates,
    tol = nmf_cv_tol,
    maxit = nmf_cv_maxit,
    verbose = 1,
    L1 = nmf_l1,
    L2 = nmf_l2,
    threads = nmf_threads,
    test_density = nmf_test_density,
    precision = nmf_precision
  )
  cv_df <- as.data.frame(cv_data)
  utils::write.csv(cv_df, file.path(method_dir, "nmf_cv_rank_errors.csv"), row.names = FALSE)

  cv_summary <- cv_df %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(
      mean_test_error = mean(test_error, na.rm = TRUE),
      sd_test_error = stats::sd(test_error, na.rm = TRUE),
      n_replicates = sum(is.finite(test_error)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(k)
  finite_cv_summary <- cv_summary[is.finite(cv_summary$mean_test_error), , drop = FALSE]
  if (nrow(finite_cv_summary) == 0) {
    stop("NMF cross-validation did not return any finite test reconstruction errors.", call. = FALSE)
  }
  min_error <- min(finite_cv_summary$mean_test_error)
  eligible_ranks <- finite_cv_summary$k[finite_cv_summary$mean_test_error <= min_error * (1 + nmf_rank_error_tolerance)]
  best_k <- min(eligible_ranks)
  cv_summary <- cv_summary %>%
    dplyr::mutate(
      selected_rank = k == best_k,
      rank_error_tolerance = nmf_rank_error_tolerance
    )
  utils::write.csv(cv_summary, file.path(method_dir, "nmf_cv_rank_summary.csv"), row.names = FALSE)

  cv_plot <- ggplot2::ggplot(cv_df, ggplot2::aes(x = k, y = test_error, color = rep, group = rep)) +
    ggplot2::geom_point(size = 1.6, alpha = 0.85) +
    ggplot2::geom_line(alpha = 0.55) +
    ggplot2::geom_line(
      data = cv_summary,
      ggplot2::aes(x = k, y = mean_test_error, group = 1),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = cv_summary,
      ggplot2::aes(x = k, y = mean_test_error),
      inherit.aes = FALSE,
      color = "black",
      size = 2
    ) +
    ggplot2::geom_vline(xintercept = best_k, linetype = "dashed", color = "red") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "NMF rank",
      y = "test reconstruction error",
      color = "CV rep",
      caption = paste0("selected k = ", best_k)
    )
  if (!identical(plot_mode, "none")) {
    ggplot2::ggsave(file.path(method_dir, "nmf_cv_rank_plot.pdf"), cv_plot, width = 5.5, height = 4.5)
  }

  message("Running final singlet NMF with k = ", best_k)
  set.seed(nmf_seed)
  nmf_model <- singlet::run_nmf(
    A,
    rank = best_k,
    tol = nmf_final_tol,
    maxit = nmf_final_maxit,
    verbose = FALSE,
    L1 = nmf_l1,
    L2 = nmf_l2,
    threads = nmf_threads,
    precision = nmf_precision
  )
  factor_names <- paste0("nmf_", seq_len(best_k))
  rownames(nmf_model$w) <- rownames(A)
  colnames(nmf_model$w) <- factor_names
  rownames(nmf_model$h) <- factor_names
  colnames(nmf_model$h) <- colnames(A)

  embeddings <- t(nmf_model$h)
  rownames(embeddings) <- colnames(A)
  colnames(embeddings) <- factor_names

  write_matrix_csv(embeddings, file.path(method_dir, "nmf_embeddings.csv"), row_id = "cell_ID")
  write_matrix_csv(nmf_model$w, file.path(method_dir, "nmf_gene_loadings.csv"), row_id = "gene")
  readr::write_rds(nmf_model, file.path(method_dir, "nmf_model.rds.gz"), compress = "gz")
  writeLines(as.character(best_k), con = file.path(method_dir, "nmf_selected_k.txt"))

  top_gene_rows <- lapply(factor_names, function(factor_name) {
    values <- nmf_model$w[, factor_name]
    top_idx <- head(order(values, decreasing = TRUE), min(50, length(values)))
    data.frame(
      component = factor_name,
      rank = seq_along(top_idx),
      gene = rownames(nmf_model$w)[top_idx],
      loading = values[top_idx],
      row.names = NULL
    )
  })
  utils::write.csv(dplyr::bind_rows(top_gene_rows), file.path(method_dir, "nmf_top_genes_by_factor.csv"), row.names = FALSE)

  list(
    embeddings = embeddings,
    metadata = list(
      source = paste0(nmf_assay, ":", nmf_layer),
      selected_k = best_k,
      selected_genes = length(selected_genes),
      cv_summary = cv_summary
    ),
    cache_key = make_cache_key(
      method = "nmf",
      source = paste0(nmf_assay, ":", nmf_layer),
      selected_k = best_k,
      selected_genes = paste(selected_genes, collapse = ";"),
      ranks = paste(nmf_ranks, collapse = ","),
      cv_replicates = nmf_cv_replicates,
      cv_tol = nmf_cv_tol,
      cv_maxit = nmf_cv_maxit,
      final_tol = nmf_final_tol,
      final_maxit = nmf_final_maxit,
      test_density = nmf_test_density,
      l1 = nmf_l1,
      l2 = nmf_l2,
      precision = nmf_precision,
      seed = nmf_seed
    )
  )
}

make_marker_module_embeddings <- function(method_dir) {
  catalog <- read_lsgi_catalog(
    "marker_modules.tsv",
    required_cols = c("module_id", "label", "category", "scorer", "genes_up", "genes_down", "description")
  )
  selected <- select_catalog_rows(catalog, "module_id", marker_module_ids, "marker module")
  feature_defs <- data.frame(
    feature_id = selected$module_id,
    label = selected$label,
    category = selected$category,
    scorer = selected$scorer,
    genes_up = selected$genes_up,
    genes_down = selected$genes_down,
    description = selected$description,
    stringsAsFactors = FALSE
  )
  score_gene_set_embeddings(
    feature_defs = feature_defs,
    method_dir = method_dir,
    method_id = "marker_module",
    source_label = "resources/lsgi_embedding_catalog/marker_modules.tsv",
    catalog_path = attr(catalog, "catalog_path")
  )
}

msigdbr_table_cache <- NULL
get_msigdbr_table <- function() {
  if (!is.null(msigdbr_table_cache)) {
    return(msigdbr_table_cache)
  }
  msig <- msigdbr::msigdbr(species = "Homo sapiens")
  required_cols <- c("gs_name", "gene_symbol")
  missing_cols <- setdiff(required_cols, colnames(msig))
  if (length(missing_cols) > 0) {
    stop("msigdbr result is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  msigdbr_table_cache <<- msig
  msig
}

make_pathway_embeddings <- function(method_dir) {
  catalog <- read_lsgi_catalog(
    "pathways.tsv",
    required_cols = c("pathway_id", "label", "source", "collection", "name", "scorer", "description")
  )
  selected <- select_catalog_rows(catalog, "pathway_id", pathway_ids, "pathway")
  unsupported_sources <- setdiff(tolower(selected$source), "msigdbr")
  if (length(unsupported_sources) > 0) {
    stop("Unsupported pathway source(s): ", paste(unique(unsupported_sources), collapse = ", "), call. = FALSE)
  }

  msig <- get_msigdbr_table()
  msig_names <- unique(msig$gs_name)
  missing_gene_sets <- setdiff(selected$name, msig_names)
  if (length(missing_gene_sets) > 0) {
    stop("Pathway gene set(s) not found in msigdbr: ", paste(missing_gene_sets, collapse = ", "), call. = FALSE)
  }
  gene_sets <- split(msig$gene_symbol, msig$gs_name)
  pathway_genes <- lapply(selected$name, function(gs_name) unique(gene_sets[[gs_name]]))
  names(pathway_genes) <- selected$pathway_id

  pathway_gene_rows <- dplyr::bind_rows(lapply(seq_len(nrow(selected)), function(i) {
    data.frame(
      pathway_id = selected$pathway_id[i],
      gs_name = selected$name[i],
      gene_symbol = pathway_genes[[selected$pathway_id[i]]],
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(pathway_gene_rows, file.path(method_dir, "pathway_gene_sets.csv"), row.names = FALSE)

  feature_defs <- data.frame(
    feature_id = selected$pathway_id,
    label = selected$label,
    source = selected$source,
    collection = selected$collection,
    gs_name = selected$name,
    scorer = selected$scorer,
    genes_up = vapply(pathway_genes, function(x) paste(x, collapse = ","), character(1)),
    genes_down = "",
    description = selected$description,
    stringsAsFactors = FALSE
  )
  result <- score_gene_set_embeddings(
    feature_defs = feature_defs,
    method_dir = method_dir,
    method_id = "pathway",
    source_label = "msigdbr Homo sapiens",
    catalog_path = attr(catalog, "catalog_path")
  )
  result$metadata$msigdbr_db_version <- if ("db_version" %in% colnames(msig)) {
    paste(unique(msig$db_version), collapse = ", ")
  } else {
    "unknown"
  }
  result
}

make_single_gene_embeddings <- function(method_dir) {
  catalog <- read_lsgi_catalog(
    "single_genes.tsv",
    required_cols = c("gene_id", "gene_symbol", "panel", "scorer", "label", "description")
  )
  selected <- select_catalog_rows(catalog, "gene_id", single_gene_ids, "single gene")
  selected$scorer <- tolower(trimws(as.character(selected$scorer)))
  unsupported_scorers <- setdiff(unique(selected$scorer), "lognorm_zscore")
  if (length(unsupported_scorers) > 0) {
    stop("Unsupported single-gene scorer(s): ", paste(unsupported_scorers, collapse = ", "), call. = FALSE)
  }

  norm_mat <- get_lsgi_expression_matrix()
  matches <- match_expression_genes(selected$gene_symbol, rownames(norm_mat))
  match_lookup <- matches[match(selected$gene_symbol, matches$requested_gene), , drop = FALSE]
  component_names <- sanitize_component_names(selected$gene_id, prefix = "single_gene_")
  usable <- match_lookup$present
  selected_summary <- data.frame(
    gene_id = selected$gene_id,
    component = component_names,
    gene_symbol = selected$gene_symbol,
    matched_gene = match_lookup$matched_gene,
    panel = selected$panel,
    scorer = selected$scorer,
    label = selected$label,
    description = selected$description,
    used_in_lsgi = usable,
    skip_reason = ifelse(usable, "", "gene_not_found_in_expression_matrix"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(selected_summary, file.path(method_dir, "selected_features.csv"), row.names = FALSE)
  missing_df <- selected_summary[!selected_summary$used_in_lsgi, c("gene_id", "gene_symbol"), drop = FALSE]
  if (nrow(missing_df) == 0) {
    missing_df <- data.frame(gene_id = character(), gene_symbol = character(), stringsAsFactors = FALSE)
  }
  utils::write.csv(missing_df, file.path(method_dir, "missing_features.csv"), row.names = FALSE)
  if (!any(usable)) {
    stop("No usable single genes after matching genes to expression matrix.", call. = FALSE)
  }

  z_scores <- zscore_selected_genes(norm_mat, match_lookup$gene_index[usable])
  embeddings <- matrix(
    NA_real_,
    nrow = ncol(norm_mat),
    ncol = sum(usable),
    dimnames = list(colnames(norm_mat), component_names[usable])
  )
  for (i in which(usable)) {
    embeddings[, component_names[i]] <- as.numeric(z_scores[match_lookup$matched_gene[i], ])
  }
  if (any(!is.finite(embeddings))) {
    stop("single_gene scorer produced non-finite embedding values.", call. = FALSE)
  }
  nonzero_score <- colSums(abs(embeddings), na.rm = TRUE) > 0
  if (any(!nonzero_score)) {
    zero_components <- colnames(embeddings)[!nonzero_score]
    zero_idx <- match(zero_components, selected_summary$component)
    selected_summary$used_in_lsgi[zero_idx] <- FALSE
    selected_summary$skip_reason[zero_idx] <- "zero_variance_or_zero_score"
    embeddings <- embeddings[, nonzero_score, drop = FALSE]
  }
  if (ncol(embeddings) == 0) {
    stop("No non-zero single-gene embeddings available for LSGI.", call. = FALSE)
  }
  utils::write.csv(selected_summary, file.path(method_dir, "selected_features.csv"), row.names = FALSE)
  write_matrix_csv(embeddings, file.path(method_dir, "embedding_matrix.csv"), row_id = "cell_ID")
  utils::write.csv(
    data.frame(
      method_id = "single_gene",
      assay = lsgi_expression_assay,
      layer = lsgi_expression_layer,
      scale_factor = lsgi_expression_scale_factor,
      catalog_path = attr(catalog, "catalog_path"),
      n_requested_features = nrow(selected),
      n_used_features = ncol(embeddings),
      scorers = paste(unique(selected$scorer[usable]), collapse = ","),
      stringsAsFactors = FALSE
    ),
    file.path(method_dir, "score_parameters.csv"),
    row.names = FALSE
  )

  list(
    embeddings = embeddings,
    metadata = list(
      source = paste0(lsgi_expression_assay, ":", lsgi_expression_layer),
      selected_ids = selected_summary$gene_id[selected_summary$used_in_lsgi],
      selected_features = ncol(embeddings),
      scorer = paste(unique(selected_summary$scorer[selected_summary$used_in_lsgi]), collapse = ", "),
      catalog_path = attr(catalog, "catalog_path")
    ),
    cache_key = make_cache_key(
      method = "single_gene",
      source = paste0(lsgi_expression_assay, ":", lsgi_expression_layer),
      scale_factor = lsgi_expression_scale_factor,
      features = paste(
        selected_summary$gene_id,
        selected_summary$gene_symbol,
        selected_summary$matched_gene,
        selected_summary$used_in_lsgi,
        sep = ":",
        collapse = "|"
      )
    )
  )
}

write_component_outputs <- function(method_id, method_label, embedding_result) {
  method_dir <- file.path(out_dir, method_id)
  arrow_dir <- file.path(method_dir, "arrow_tables")
  plot_dir <- file.path(method_dir, "plots")
  dir.create(arrow_dir, recursive = TRUE, showWarnings = FALSE)
  if (!identical(plot_mode, "none")) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  lsgi_result_path <- file.path(paths$intermediate, paste0("11_lsgi_", method_id, "_result.rds.gz"))
  legacy_lsgi_result_path <- file.path(paths$intermediate, "11_lsgi_cell_component_result.rds.gz")
  did_reuse_lsgi <- FALSE
  lsgi_res <- NULL
  if (reuse_lsgi && file.exists(lsgi_result_path)) {
    cached_lsgi_res <- readr::read_rds(lsgi_result_path)
    if (cache_is_compatible(cached_lsgi_res, embedding_result, grid_cache_key, allow_legacy = identical(method_id, "cell_component"))) {
      lsgi_res <- cached_lsgi_res
      did_reuse_lsgi <- TRUE
    } else {
      message("Ignoring stale LSGI cache for method ", method_id, " because its embedding/grid/spot signature does not match.")
    }
  }
  if (is.null(lsgi_res) && reuse_lsgi && identical(method_id, "cell_component") && file.exists(legacy_lsgi_result_path)) {
    cached_lsgi_res <- readr::read_rds(legacy_lsgi_result_path)
    if (cache_is_compatible(cached_lsgi_res, embedding_result, grid_cache_key, allow_legacy = TRUE)) {
      lsgi_res <- cached_lsgi_res
      did_reuse_lsgi <- TRUE
    } else {
      message("Ignoring stale legacy LSGI cache for cell_component because its embedding/grid/spot signature does not match.")
    }
  }
  if (is.null(lsgi_res)) {
    lsgi_res <- make_lsgi_result(embedding_result$embeddings)
    lsgi_res$component_method <- method_id
    lsgi_res$component_method_cache_key <- embedding_result$cache_key
    lsgi_res$grid_cache_key <- grid_cache_key
    readr::write_rds(lsgi_res, lsgi_result_path, compress = "gz")
    if (identical(method_id, "cell_component")) {
      readr::write_rds(lsgi_res, legacy_lsgi_result_path, compress = "gz")
    }
  }

  method_grid_info_df <- lsgi_res$grid.info %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      grid = grid_ids,
      grid_index = seq_len(dplyr::n())
    ) %>%
    dplyr::select(grid, grid_index, X, Y, vx, vy, R_squared, Assignment, qsum, sf, vx.u, vy.u)
  utils::write.csv(method_grid_info_df, file.path(method_dir, "grid_info.csv"), row.names = FALSE)

  lin_res <- get.ind.rsqrs(lsgi_res)
  lin_res <- stats::na.omit(lin_res)
  lin_res <- rescale_arrow_length(lin_res)
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
  method_grid_info_df <- method_grid_info_df %>%
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
  utils::write.csv(method_grid_info_df, file.path(method_dir, "grid_info.csv"), row.names = FALSE)

  arrows_by_grid <- all_arrow_df[, c(
    "grid", "component", "X", "Y", "vx", "vy", "rsquared", "raw_length",
    "gradient_strength", "gradient_strength_norm", "arrow_length_multiplier",
    "arrow_head_angle_mapped",
    "sf", "original_vx.u", "original_vy.u", "original_scaled_length",
    "vx.u", "vy.u", "scaled_length", "plotted_length",
    "passes_r_squared", "passes_minimum_fctr", "included_in_filtered_output"
  )]
  utils::write.csv(arrows_by_grid, file.path(arrow_dir, "arrows_by_grid.csv"), row.names = FALSE)

  arrow_df <- all_arrow_df[all_arrow_df$included_in_filtered_output, c(
    "vx", "vy", "rsquared", "component", "grid", "X", "Y", "qsum", "sf",
    "gradient_strength", "gradient_strength_norm", "arrow_length_multiplier",
    "arrow_head_angle_mapped",
    "original_vx.u", "original_vy.u", "original_scaled_length",
    "vx.u", "vy.u", "scaled_length", "plotted_length"
  ), drop = FALSE]
  colnames(arrow_df)[colnames(arrow_df) == "component"] <- "fctr"
  utils::write.csv(arrow_df, file.path(method_dir, "gradient_arrows.csv"), row.names = FALSE)

  dist_mat <- tryCatch(
    calc_gradient_distance_fast(lin_res),
    error = function(e) {
      warning("LSGI distance calculation skipped for ", method_id, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (!is.null(dist_mat) && nrow(dist_mat) > 0) {
    utils::write.csv(dist_mat, file.path(method_dir, "gradient_distance.csv"), row.names = FALSE)
    if (!identical(plot_mode, "none")) {
      write_distance_heatmap(dist_mat, file.path(method_dir, "gradient_distance_heatmap.pdf"), method_id)
    }
  }

  if (!identical(plot_mode, "none")) {
  arrow_plot_df <- sp_prepare_arrows(
    arrow_df,
    arrow_length_scale = arrow_length_scale,
    flip_y = TRUE
  )

  add_gradient_arrows <- function(p, arrow_data = arrow_plot_df) {
    if (nrow(arrow_data) == 0) {
      return(p + ggplot2::labs(subtitle = paste0("No ", method_label, " gradients passed R2 > ", r_squared_thresh)))
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
    out + ggplot2::labs(color = method_label)
  }

  save_spatial_method_plot <- function(plot, filename) {
    sp_save_plot(
      plot,
      filename,
      width = params$spatial_plot_width %||% 8,
      height = params$spatial_plot_height %||% 7,
      map_width = params$spatial_plot_map_width %||% 6.45,
      legend_width = params$spatial_plot_legend_width %||% 1.55
    )
  }

  boundary_base <- sp_boundary_base(
    point_spot_polygons,
    plot_context,
    boundary_cols,
    title = paste0(sample_name, " boundary with LSGI ", method_label, " gradients"),
    show_image = FALSE
  )

  boundary_gradient <- add_gradient_arrows(boundary_base)
  boundary_pdf <- file.path(plot_dir, paste0(sample_name, "_BoundaryDefine_LSGIGradient.pdf"))
  save_spatial_method_plot(boundary_gradient, boundary_pdf)

  he_pdf <- NULL
  if (has_image) {
    he_base <- sp_boundary_base(
      point_spot_polygons,
      plot_context,
      boundary_cols,
      title = paste0(sample_name, " HE-boundary with LSGI ", method_label, " gradients"),
      show_image = TRUE
    )

    he_gradient <- add_gradient_arrows(he_base, arrow_plot_df)
    he_pdf <- file.path(plot_dir, paste0(sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf"))
    save_spatial_method_plot(he_gradient, he_pdf)
  }

  plain_gradient_base <- sp_spatial_canvas(
    plot_context,
    show_image = FALSE,
    title = paste0("LSGI ", method_label, " gradients")
  ) +
    sp_spot_layer(
      point_spot_polygons,
      fill = "lightgrey",
      colour = NA,
      alpha = 1,
      linewidth = 0
    )
  plain_pdf <- file.path(plot_dir, "gradients_plain_lsgi.pdf")
  save_spatial_method_plot(add_gradient_arrows(plain_gradient_base), plain_pdf)

  legacy_boundary_pdf <- file.path(method_dir, paste0(sample_name, "_BoundaryDefine_LSGIGradient.pdf"))
  legacy_he_pdf <- file.path(method_dir, paste0(sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf"))
  legacy_plain_pdf <- file.path(method_dir, "gradients_plain_lsgi.pdf")
  file.copy(boundary_pdf, legacy_boundary_pdf, overwrite = TRUE)
  if (!is.null(he_pdf)) {
    file.copy(he_pdf, legacy_he_pdf, overwrite = TRUE)
  }
  file.copy(plain_pdf, legacy_plain_pdf, overwrite = TRUE)

  component_levels <- colnames(lsgi_res$embeddings)
  component_files <- sanitize_filename(component_levels, prefix = "component")
  if (identical(plot_mode, "full")) for (component_idx in seq_along(component_levels)) {
    component_id <- component_levels[[component_idx]]
    component_file <- component_files[[component_idx]]
    component_arrow_df <- arrow_plot_df[as.character(arrow_plot_df$fctr) == component_id, , drop = FALSE]

    component_boundary <- add_gradient_arrows(
      boundary_base + ggplot2::ggtitle(paste0(sample_name, " boundary with LSGI ", method_label, " gradient: ", component_id)),
      component_arrow_df
    )
    save_spatial_method_plot(
      component_boundary,
      file.path(plot_dir, paste0(sample_name, "_", method_id, "_", component_file, "_BoundaryDefine_LSGIGradient.pdf"))
    )

    if (has_image) {
      component_he <- add_gradient_arrows(
        he_base + ggplot2::ggtitle(paste0(sample_name, " HE-boundary with LSGI ", method_label, " gradient: ", component_id)),
        component_arrow_df
      )
      save_spatial_method_plot(
        component_he,
        file.path(plot_dir, paste0(sample_name, "_", method_id, "_", component_file, "_BoundaryDefine_HE_LSGIGradient.pdf"))
      )
    }
  }
  }

  if (identical(method_id, "cell_component") && write_legacy_cell_component_outputs) {
    legacy_arrow_dir <- file.path(out_dir, "arrow_tables")
    dir.create(legacy_arrow_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(arrows_by_grid, file.path(legacy_arrow_dir, "cell_component_arrows_by_grid.csv"), row.names = FALSE)
    utils::write.csv(arrow_df, file.path(out_dir, "cell_component_gradient_arrows.csv"), row.names = FALSE)
    if (!is.null(dist_mat) && nrow(dist_mat) > 0) {
      utils::write.csv(dist_mat, file.path(out_dir, "cell_component_gradient_distance.csv"), row.names = FALSE)
      if (!identical(plot_mode, "none")) {
        write_distance_heatmap(dist_mat, file.path(out_dir, "cell_component_gradient_distance_heatmap.pdf"), "cell_component legacy")
      }
    }
    if (!identical(plot_mode, "none")) {
      file.copy(boundary_pdf, file.path(out_dir, paste0(sample_name, "_BoundaryDefine_LSGIGradient.pdf")), overwrite = TRUE)
      if (!is.null(he_pdf)) {
        file.copy(he_pdf, file.path(out_dir, paste0(sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf")), overwrite = TRUE)
      }
      file.copy(plain_pdf, file.path(out_dir, "cell_component_gradients_plain_lsgi.pdf"), overwrite = TRUE)
    }
  }

  sink(file.path(method_dir, "method_summary.txt"))
  cat("LSGI component method summary\n")
  cat("=============================\n\n")
  cat("Sample:", sample_name, "\n")
  cat("Method ID:", method_id, "\n")
  cat("Method label:", method_label, "\n")
  cat("Embedding source:", embedding_result$metadata$source %||% "unspecified", "\n")
  cat("Matched spots:", nrow(spatial_coords), "\n")
  cat("Components:", paste(colnames(lsgi_res$embeddings), collapse = ", "), "\n")
  cat("Number of components:", ncol(lsgi_res$embeddings), "\n")
  cat("Selected gradient arrows:", nrow(arrow_df), "\n")
  cat("Generated grids:", nrow(method_grid_info_df), "\n")
  cat("Reused LSGI result:", did_reuse_lsgi, "\n\n")
  if (!is.null(embedding_result$metadata$catalog_path)) {
    cat("Catalog:", embedding_result$metadata$catalog_path, "\n")
  }
  if (!is.null(embedding_result$metadata$scorer)) {
    cat("Scorer(s):", embedding_result$metadata$scorer, "\n")
  }
  if (!is.null(embedding_result$metadata$selected_features)) {
    cat("Selected features:", embedding_result$metadata$selected_features, "\n")
  }
  if (!is.null(embedding_result$metadata$selected_ids)) {
    cat("Selected IDs:", paste(embedding_result$metadata$selected_ids, collapse = ", "), "\n")
  }
  if (!is.null(embedding_result$metadata$msigdbr_db_version)) {
    cat("MSigDB version:", embedding_result$metadata$msigdbr_db_version, "\n")
  }
  if (
    !is.null(embedding_result$metadata$catalog_path) ||
      !is.null(embedding_result$metadata$scorer) ||
      !is.null(embedding_result$metadata$selected_features) ||
      !is.null(embedding_result$metadata$selected_ids) ||
      !is.null(embedding_result$metadata$msigdbr_db_version)
  ) {
    cat("\n")
  }
  if (identical(method_id, "nmf")) {
    cat("NMF selected k:", embedding_result$metadata$selected_k, "\n")
    cat("NMF selected genes:", embedding_result$metadata$selected_genes, "\n")
    cat("NMF ranks scanned:", paste(nmf_ranks, collapse = ", "), "\n")
    cat("NMF CV replicates:", nmf_cv_replicates, "\n")
    cat("NMF CV tol:", nmf_cv_tol, "\n")
    cat("NMF final tol:", nmf_final_tol, "\n")
    cat("NMF rank error tolerance:", nmf_rank_error_tolerance, "\n\n")
  }
  cat("Outputs:\n")
  cat("- grid_info.csv\n")
  cat("- arrow_tables/arrows_by_grid.csv\n")
  cat("- gradient_arrows.csv\n")
  cat("- gradient_distance.csv / gradient_distance_heatmap.pdf when calculable\n")
  cat("- gradients_plain_lsgi.pdf\n")
  cat("- *_BoundaryDefine_LSGIGradient.pdf\n")
  cat("- *_BoundaryDefine_HE_LSGIGradient.pdf when H&E image is available\n")
  sink()

  list(
    method_id = method_id,
    method_label = method_label,
    n_components = ncol(lsgi_res$embeddings),
    n_arrows_all = nrow(arrows_by_grid),
    n_arrows_filtered = nrow(arrow_df),
    lsgi_result_path = lsgi_result_path,
    method_dir = method_dir,
    did_reuse_lsgi = did_reuse_lsgi
  )
}

if (plot_common_grid_plots) {
  plot_grid_memberships()
} else {
  message("Skipping common grid membership PDF plots (plot_mode=", plot_mode, ").")
}

method_summaries <- list()
for (method_id in component_methods) {
  method_dir <- file.path(out_dir, method_id)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  if (identical(method_id, "cell_component")) {
    message("Preparing cell-component embeddings.")
    embedding_result <- make_cell_component_embeddings()
    method_label <- "cell component"
  } else if (identical(method_id, "nmf")) {
    message("Preparing NMF embeddings.")
    embedding_result <- make_nmf_embeddings(method_dir)
    method_label <- "NMF"
  } else if (identical(method_id, "marker_module")) {
    message("Preparing marker-module embeddings.")
    embedding_result <- make_marker_module_embeddings(method_dir)
    method_label <- "marker module"
  } else if (identical(method_id, "pathway")) {
    message("Preparing pathway embeddings.")
    embedding_result <- make_pathway_embeddings(method_dir)
    method_label <- "pathway"
  } else if (identical(method_id, "single_gene")) {
    message("Preparing single-gene embeddings.")
    embedding_result <- make_single_gene_embeddings(method_dir)
    method_label <- "single gene"
  } else {
    stop("Unsupported component method: ", method_id, call. = FALSE)
  }

  message("Running LSGI arrow analysis for method: ", method_id)
  method_summaries[[method_id]] <- write_component_outputs(method_id, method_label, embedding_result)
}

method_summary_df <- dplyr::bind_rows(lapply(method_summaries, as.data.frame))
utils::write.csv(method_summary_df, file.path(out_dir, "component_method_summary.csv"), row.names = FALSE)

run_summary_path <- file.path(out_dir, "run_summary_11.txt")
sink(run_summary_path)
cat("Cottrazm + LSGI multi-embedding gradient analysis\n")
cat("=================================================\n\n")
cat("Sample:", sample_name, "\n")
cat("Matched spots:", nrow(spatial_coords), "\n")
cat("Component methods:", paste(component_methods, collapse = ", "), "\n")
cat("Embedding catalog dir:", normalizePath(embedding_catalog_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Expression assay/layer:", lsgi_expression_assay, "/", lsgi_expression_layer, "\n", sep = "")
cat("Expression scale factor:", lsgi_expression_scale_factor, "\n")
cat("Minimum feature genes:", lsgi_min_feature_genes, "\n")
cat("Grid clustering backend:", grid_clustering_backend, "\n")
cat("n.grids.scale:", n_grids_scale, "\n")
cat("n.cells.per.meta:", n_cells_per_meta, "\n")
cat("Generated grids:", nrow(grid_info_common), "\n")
cat("H&E overlay available:", has_image, "\n")
cat("Low-resolution spot diameter px:", plot_context$spot_diameter_lowres, "\n")
cat("Spatial canvas inches:", params$spatial_plot_width %||% 8, "x", params$spatial_plot_height %||% 7, "\n\n")
cat("Common grid PDF plots:", plot_common_grid_plots, "\n\n")
cat("Arrow filters and styling:\n")
cat("R-squared threshold:", r_squared_thresh, "\n")
cat("Minimum arrows per component:", minimum_fctr, "\n")
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
cat("Arrow closed:", arrow_closed, "\n\n")
cat("Method summaries:\n")
print(method_summary_df)
cat("\nCommon outputs:\n")
cat("- output/11_lsgi_gradient/grid_info.csv\n")
cat("- output/11_lsgi_gradient/spot_grid_membership.csv\n")
cat("- output/11_lsgi_gradient/grid_spot_summary.csv\n")
cat("- output/11_lsgi_gradient/grid_local_spot_membership.csv\n")
cat("- output/11_lsgi_gradient/grid_partition_membership.csv\n")
cat("- output/11_lsgi_gradient/grid_partition.pdf\n")
cat("- output/11_lsgi_gradient/grid_partition_overlay.pdf\n")
cat("- output/11_lsgi_gradient/grid_local_spots.pdf\n")
cat("- output/11_lsgi_gradient/grid_local_spots_overlay.pdf\n")
cat("- output/11_lsgi_gradient/component_method_summary.csv\n\n")
cat("Method output directories:\n")
for (method_id in component_methods) {
  cat("- output/11_lsgi_gradient/", method_id, "/\n", sep = "")
}
if (write_legacy_cell_component_outputs && "cell_component" %in% component_methods) {
  cat("\nLegacy cell-component compatibility outputs were also refreshed for current step 12 defaults.\n")
}
sink()
file.copy(run_summary_path, file.path(out_dir, "run_summary_11b.txt"), overwrite = TRUE)

message("Done. LSGI outputs written to ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), ".")
