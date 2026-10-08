# Rotation-invariant spatial adjacency for standard hexagonal Visium arrays.
# Space Ranger array coordinates encode the six immediate neighbors exactly:
#   (row, col +/- 2) and (row +/- 1, col +/- 1).

visium_hex_neighbors_api_version <- "1.0.0"

vh_build_neighbors <- function(position) {
  required <- c("row", "col")
  missing_columns <- setdiff(required, colnames(position))
  if (length(missing_columns)) {
    stop("Position table is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (is.null(rownames(position)) || any(!nzchar(rownames(position))) || anyDuplicated(rownames(position))) {
    stop("Position table must have unique non-empty spot IDs as row names.", call. = FALSE)
  }

  row_value <- as.numeric(position$row)
  col_value <- as.numeric(position$col)
  if (any(!is.finite(row_value)) || any(!is.finite(col_value)) ||
      any(abs(row_value - round(row_value)) > 1e-8) ||
      any(abs(col_value - round(col_value)) > 1e-8)) {
    stop("Visium array row/column coordinates must be finite integers.", call. = FALSE)
  }
  row_value <- as.integer(round(row_value))
  col_value <- as.integer(round(col_value))
  coordinate_key <- paste(row_value, col_value, sep = ":")
  if (anyDuplicated(coordinate_key)) {
    stop("Duplicated Visium array row/column coordinate.", call. = FALSE)
  }

  spot_ids <- rownames(position)
  spot_by_coordinate <- stats::setNames(spot_ids, coordinate_key)
  offsets <- rbind(
    c(0L, -2L), c(0L, 2L),
    c(-1L, -1L), c(-1L, 1L),
    c(1L, -1L), c(1L, 1L)
  )
  neighbors <- lapply(seq_along(spot_ids), function(i) {
    candidate_keys <- paste(
      row_value[[i]] + offsets[, 1L],
      col_value[[i]] + offsets[, 2L],
      sep = ":"
    )
    candidates <- spot_by_coordinate[candidate_keys]
    unname(candidates[!is.na(candidates)])
  })
  names(neighbors) <- spot_ids
  neighbors
}

vh_neighbor_qc <- function(neighbors) {
  if (!is.list(neighbors) || is.null(names(neighbors))) {
    stop("neighbors must be a named list.", call. = FALSE)
  }
  degree <- lengths(neighbors)
  symmetric <- all(vapply(names(neighbors), function(id) {
    all(vapply(neighbors[[id]], function(other) {
      other %in% names(neighbors) && id %in% neighbors[[other]]
    }, logical(1)))
  }, logical(1)))
  data.frame(
    metric = c(
      "method", "spot_count", "undirected_edge_count", "isolated_spot_count",
      "isolated_spot_fraction", "minimum_degree", "median_degree", "mean_degree",
      "maximum_degree", "symmetric"
    ),
    value = c(
      "visium_array_hex_6_neighbor", length(degree), sum(degree) / 2,
      sum(degree == 0L), mean(degree == 0L), min(degree), stats::median(degree),
      mean(degree), max(degree), symmetric
    ),
    stringsAsFactors = FALSE
  )
}

vh_assert_neighbor_qc <- function(qc, maximum_isolated_fraction = 0.5) {
  values <- stats::setNames(as.character(qc$value), qc$metric)
  edge_count <- suppressWarnings(as.numeric(values[["undirected_edge_count"]]))
  isolated_fraction <- suppressWarnings(as.numeric(values[["isolated_spot_fraction"]]))
  maximum_degree <- suppressWarnings(as.numeric(values[["maximum_degree"]]))
  if (!identical(values[["symmetric"]], "TRUE")) {
    stop("Visium hex neighbor graph is not symmetric.", call. = FALSE)
  }
  if (!is.finite(edge_count) || edge_count <= 0 ||
      !is.finite(isolated_fraction) || isolated_fraction >= maximum_isolated_fraction) {
    stop("Visium hex neighbor graph is empty or has too many isolated spots.", call. = FALSE)
  }
  if (!is.finite(maximum_degree) || maximum_degree > 6) {
    stop("Visium hex neighbor graph has a spot with more than six neighbors.", call. = FALSE)
  }
  invisible(qc)
}
