# Pure Cottrazm CNV-to-boundary core shared by the production pipeline and
# benchmark adapters. This file deliberately has no file-system, plotting,
# Seurat, or evaluator-truth dependency.

cottrazm_boundary_core_api_version <- "1.0.0"

ct_latest_file <- function(directory, pattern) {
  candidates <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (!length(candidates)) return(NA_character_)
  candidates[[order(file.info(candidates)$mtime, basename(candidates),
                     decreasing = TRUE, method = "radix")[[1L]]]]
}

ct_read_infercnv_outputs <- function(out_dir, spot_ids, reference_ids,
                                     cnv_k = 4L, neutral_state = 3,
                                     seed = 666L) {
  spot_ids <- as.character(spot_ids)
  reference_ids <- intersect(as.character(reference_ids), spot_ids)
  group_file <- ct_latest_file(out_dir, "cell_groupings$")
  gene_file <- ct_latest_file(out_dir, "pred_cnv_genes\\.dat$")
  region_file <- ct_latest_file(out_dir, "pred_cnv_regions\\.dat$")
  if (is.na(group_file) || is.na(gene_file)) {
    stop("inferCNV HMM grouping/gene artifacts are missing", call. = FALSE)
  }
  groups <- utils::read.delim(group_file, check.names = FALSE,
                              stringsAsFactors = FALSE)
  genes <- utils::read.delim(gene_file, check.names = FALSE,
                             stringsAsFactors = FALSE)
  required_groups <- c("cell_group_name", "cell")
  required_genes <- c("cell_group_name", "state")
  if (!all(required_groups %in% names(groups)) ||
      !all(required_genes %in% names(genes))) {
    stop("inferCNV HMM artifact schema mismatch", call. = FALSE)
  }
  groups$cell <- as.character(groups$cell)
  groups$cell_group_name <- as.character(groups$cell_group_name)
  genes$cell_group_name <- as.character(genes$cell_group_name)
  genes$state <- as.numeric(genes$state)
  group_scores <- tapply(abs(genes$state - as.numeric(neutral_state)),
                         genes$cell_group_name, sum, na.rm = TRUE)
  groups$cnv_score <- as.numeric(group_scores[groups$cell_group_name])
  groups$cnv_score[!is.finite(groups$cnv_score)] <- 0
  observation_groups <- unique(groups$cell_group_name[!groups$cell %in% reference_ids])
  observation_scores <- tapply(
    groups$cnv_score[groups$cell_group_name %in% observation_groups],
    groups$cell_group_name[groups$cell_group_name %in% observation_groups],
    mean
  )
  observation_scores <- observation_scores[is.finite(observation_scores)]
  if (!length(observation_scores)) {
    stop("inferCNV produced no observation group score", call. = FALSE)
  }
  centers <- min(as.integer(cnv_k), length(unique(observation_scores)))
  if (centers <= 1L) {
    group_labels <- setNames(rep("1", length(observation_scores)),
                             names(observation_scores))
  } else {
    set.seed(as.integer(seed))
    fit <- stats::kmeans(as.matrix(observation_scores), centers = centers,
                         iter.max = 100L, nstart = 1L, algorithm = "Lloyd")
    group_labels <- setNames(as.character(fit$cluster), names(observation_scores))
  }
  groups$CNVLabel <- as.character(group_labels[groups$cell_group_name])
  groups$CNVLabel[groups$cell %in% reference_ids] <- "Normal"
  groups$CNVLabel[is.na(groups$CNVLabel)] <- "Normal"
  calls <- data.frame(
    cell_ID = spot_ids,
    CNVLabel = "Normal",
    cnv_score = 0,
    infercnv_role = ifelse(spot_ids %in% reference_ids, "Reference", "Observation"),
    stringsAsFactors = FALSE
  )
  matched <- match(calls$cell_ID, groups$cell)
  calls$CNVLabel[!is.na(matched)] <- groups$CNVLabel[matched[!is.na(matched)]]
  calls$cnv_score[!is.na(matched)] <- groups$cnv_score[matched[!is.na(matched)]]
  calls$CNVLabel[calls$infercnv_role == "Reference"] <- "Normal"
  calls$cnv_score[calls$infercnv_role == "Reference"] <- 0

  values <- data.frame(
    spot_id = character(), segment_id = character(), chromosome = character(),
    start = numeric(), end = numeric(), inferred_CNV_value = numeric(),
    inferred_cluster = character(), stringsAsFactors = FALSE
  )
  if (!is.na(region_file)) {
    regions <- utils::read.delim(region_file, check.names = FALSE,
                                 stringsAsFactors = FALSE)
    required_regions <- c("cell_group_name", "cnv_name", "state", "chr", "start", "end")
    if (all(required_regions %in% names(regions)) && nrow(regions)) {
      split_cells <- split(groups$cell, groups$cell_group_name)
      pieces <- lapply(seq_len(nrow(regions)), function(i) {
        cells <- split_cells[[as.character(regions$cell_group_name[[i]])]]
        cells <- intersect(as.character(cells), spot_ids)
        if (!length(cells)) return(NULL)
        data.frame(
          spot_id = cells,
          segment_id = as.character(regions$cnv_name[[i]]),
          chromosome = as.character(regions$chr[[i]]),
          start = as.numeric(regions$start[[i]]),
          end = as.numeric(regions$end[[i]]),
          inferred_CNV_value = as.numeric(regions$state[[i]]),
          inferred_cluster = as.character(group_labels[
            as.character(regions$cell_group_name[[i]])
          ]),
          stringsAsFactors = FALSE
        )
      })
      pieces <- pieces[!vapply(pieces, is.null, logical(1))]
      if (length(pieces)) values <- do.call(rbind, pieces)
    }
  }
  list(
    calls = calls,
    inferred_values = values,
    artifacts = data.frame(
      artifact = c("cell_groupings", "pred_cnv_genes", "pred_cnv_regions"),
      path = c(group_file, gene_file, region_file),
      stringsAsFactors = FALSE
    )
  )
}

ct_require_named_neighbors <- function(neighbors, spot_ids) {
  if (!is.list(neighbors) || is.null(names(neighbors)) ||
      !setequal(names(neighbors), spot_ids)) {
    stop("neighbors must be a named list with exactly one entry per spot", call. = FALSE)
  }
  unknown <- setdiff(unique(as.character(unlist(neighbors, use.names = FALSE))), spot_ids)
  if (length(unknown)) {
    stop("neighbors contain unknown spot IDs: ", paste(head(unknown, 3L), collapse = ","),
         call. = FALSE)
  }
  invisible(TRUE)
}

ct_select_malignant_labels <- function(cnv_label, cnv_score, infercnv_role = NULL,
                                       malignant_label_n = 2L,
                                       forbidden_labels = c("Normal", "Filtered")) {
  cnv_label <- as.character(cnv_label)
  cnv_score <- as.numeric(cnv_score)
  malignant_label_n <- as.integer(malignant_label_n)
  if (length(cnv_label) != length(cnv_score) || malignant_label_n < 1L) {
    stop("invalid CNV label/score inputs", call. = FALSE)
  }
  role_ok <- if (is.null(infercnv_role)) {
    rep(TRUE, length(cnv_label))
  } else {
    as.character(infercnv_role) == "Observation"
  }
  valid <- is.finite(cnv_score) & !is.na(cnv_label) & nzchar(cnv_label) &
    !cnv_label %in% forbidden_labels & role_ok
  label_scores <- vapply(
    split(cnv_score[valid], cnv_label[valid]),
    stats::median,
    FUN.VALUE = numeric(1),
    na.rm = TRUE
  )
  label_scores <- label_scores[is.finite(label_scores)]
  if (length(label_scores) < malignant_label_n) {
    stop("Fewer than ", malignant_label_n,
         " valid Observation CNV labels for malignant seeding", call. = FALSE)
  }
  ordered <- order(-label_scores, names(label_scores), method = "radix")
  names(label_scores)[ordered][seq_len(malignant_label_n)]
}

ct_frontier_update <- function(frontier_ids, neighbors, embedding, normal_ids,
                               boundary_ids, malignant_ids, round_index,
                               expand_malignant_radius = 0.8) {
  pieces <- lapply(frontier_ids, function(id) {
    adjacent <- neighbors[[id]]
    current <- embedding[id, , drop = FALSE]
    adjacent_embedding <- embedding[adjacent, , drop = FALSE]
    malignant_adjacent <- adjacent[adjacent %in% malignant_ids]
    malignant_center <- colMeans(rbind(
      adjacent_embedding[malignant_adjacent, , drop = FALSE], current
    ))
    malignant_radius <- vapply(c(malignant_adjacent, id), function(candidate) {
      sqrt(sum((embedding[candidate, ] - malignant_center)^2))
    }, numeric(1))
    boundary_adjacent <- adjacent[adjacent %in% boundary_ids]
    candidates <- adjacent[
      !adjacent %in% malignant_ids & !adjacent %in% boundary_ids &
        !adjacent %in% normal_ids
    ]
    if (!length(candidates)) {
      return(data.frame(spot_id = character(), p1 = numeric(), p2 = numeric(),
                        proposed = character(), stringsAsFactors = FALSE))
    }
    p1 <- vapply(candidates, function(candidate) {
      sqrt(sum((embedding[candidate, ] - malignant_center)^2))
    }, numeric(1))
    if (length(boundary_adjacent) <= 1L) {
      p2 <- p1
      proposed <- ifelse(
        p1 <= expand_malignant_radius * max(malignant_radius), "Mal", "Bdy"
      )
    } else {
      boundary_center <- colMeans(adjacent_embedding[boundary_adjacent, , drop = FALSE])
      boundary_radius <- vapply(boundary_adjacent, function(candidate) {
        sqrt(sum((embedding[candidate, ] - boundary_center)^2))
      }, numeric(1))
      p2 <- vapply(candidates, function(candidate) {
        sqrt(sum((embedding[candidate, ] - boundary_center)^2))
      }, numeric(1))
      proposed <- ifelse(
        p1 <= expand_malignant_radius * max(malignant_radius) &
          p2 > max(boundary_radius),
        "Mal", "Bdy"
      )
    }
    data.frame(spot_id = candidates, p1 = p1, p2 = p2, proposed = proposed,
               stringsAsFactors = FALSE)
  })
  proposals <- do.call(rbind, pieces)
  if (is.null(proposals) || !nrow(proposals)) {
    return(data.frame(spot_id = character(), location = character(),
                      stringsAsFactors = FALSE))
  }
  proposed_by_spot <- split(proposals$proposed, proposals$spot_id)
  data.frame(
    spot_id = names(proposed_by_spot),
    location = vapply(proposed_by_spot, function(value) {
      if ("Mal" %in% value) paste0("Mal", round_index) else "Bdy"
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

ct_define_boundary <- function(spot_ids, cnv_label, cnv_score, infercnv_role,
                               cluster, normal_score, embedding, neighbors,
                               malignant_labels = NULL, malignant_label_n = 2L,
                               malignant_cluster_fraction = 0.3,
                               umap_malignant_ratio = 0.5,
                               expand_malignant_radius = 1,
                               maximum_rounds = 6L) {
  spot_ids <- as.character(spot_ids)
  n <- length(spot_ids)
  if (!n || anyDuplicated(spot_ids)) stop("spot_ids must be unique and non-empty", call. = FALSE)
  vectors <- list(cnv_label, cnv_score, infercnv_role, cluster, normal_score)
  if (any(vapply(vectors, length, integer(1)) != n)) {
    stop("all spot-level inputs must have the same length", call. = FALSE)
  }
  if (!is.matrix(embedding) && !is.data.frame(embedding)) {
    stop("embedding must be a numeric matrix or data frame", call. = FALSE)
  }
  embedding <- as.matrix(embedding)
  storage.mode(embedding) <- "double"
  if (nrow(embedding) != n || ncol(embedding) < 2L || any(!is.finite(embedding))) {
    stop("embedding must contain finite rows for every spot", call. = FALSE)
  }
  rownames(embedding) <- spot_ids
  ct_require_named_neighbors(neighbors, spot_ids)
  if (!is.finite(malignant_cluster_fraction) || malignant_cluster_fraction <= 0 ||
      malignant_cluster_fraction > 1 || !is.finite(umap_malignant_ratio) ||
      umap_malignant_ratio <= 0 || !is.finite(expand_malignant_radius) ||
      expand_malignant_radius <= 0 || as.integer(maximum_rounds) < 1L) {
    stop("invalid Cottrazm boundary parameter", call. = FALSE)
  }

  cnv_label <- as.character(cnv_label)
  cnv_score <- as.numeric(cnv_score)
  infercnv_role <- as.character(infercnv_role)
  cluster <- factor(as.character(cluster))
  normal_score <- as.numeric(normal_score)
  if (is.null(malignant_labels)) {
    malignant_labels <- ct_select_malignant_labels(
      cnv_label, cnv_score, infercnv_role, malignant_label_n
    )
  } else {
    malignant_labels <- setdiff(as.character(malignant_labels), c("Normal", "Filtered"))
    if (!length(malignant_labels)) {
      stop("malignant_labels contains only forbidden labels", call. = FALSE)
    }
  }

  cluster_levels <- levels(cluster)
  cluster_mean_normal <- vapply(cluster_levels, function(level) {
    mean(normal_score[cluster == level])
  }, numeric(1))
  if (!any(is.finite(cluster_mean_normal))) {
    stop("no finite cluster normal score", call. = FALSE)
  }
  normal_cluster <- cluster_levels[order(-cluster_mean_normal, cluster_levels,
                                          method = "radix")][[1L]]
  normal_ids <- spot_ids[cluster == normal_cluster]
  malignant_seed_ids <- spot_ids[cnv_label %in% malignant_labels]

  selected_clusters <- cluster_levels[vapply(cluster_levels, function(level) {
    members <- cluster == level
    sum(cnv_label[members] %in% malignant_labels) >
      sum(members) * malignant_cluster_fraction
  }, logical(1))]
  if (!length(selected_clusters)) {
    stop("no method-side cluster passed malignant CNV fraction", call. = FALSE)
  }

  malignant_centers <- lapply(selected_clusters, function(level) {
    ids <- intersect(malignant_seed_ids, spot_ids[cluster == level])
    if (!length(ids)) return(rep(NA_real_, ncol(embedding)))
    colMeans(embedding[ids, , drop = FALSE])
  })
  names(malignant_centers) <- selected_clusters
  normal_center <- colMeans(embedding[normal_ids, , drop = FALSE])

  filtered_malignant <- unlist(lapply(selected_clusters, function(level) {
    ids <- intersect(malignant_seed_ids, spot_ids[cluster == level])
    center <- malignant_centers[[level]]
    ids[vapply(ids, function(id) {
      tumor_distance <- sqrt(sum((embedding[id, ] - center)^2))
      normal_distance <- sqrt(sum((embedding[id, ] - normal_center)^2))
      tumor_distance < umap_malignant_ratio * normal_distance
    }, logical(1))]
  }), use.names = FALSE)
  filtered_malignant <- unique(as.character(filtered_malignant))
  if (!length(filtered_malignant)) {
    stop("no malignant CNV spot passed the embedding ratio", call. = FALSE)
  }

  isolated_malignant <- filtered_malignant[vapply(filtered_malignant, function(id) {
    !any(neighbors[[id]] %in% filtered_malignant)
  }, logical(1))]
  if (!length(isolated_malignant)) {
    # The legacy code used degree zero. Dense, correctly recovered simulated
    # lesions can have no such spot, so extend the same topology rule to the
    # minimum within-malignant degree without consulting geometry truth.
    malignant_degree <- vapply(filtered_malignant, function(id) {
      sum(neighbors[[id]] %in% filtered_malignant)
    }, integer(1))
    isolated_malignant <- filtered_malignant[
      malignant_degree == min(malignant_degree)
    ]
  }
  initial_candidates <- unique(unlist(lapply(isolated_malignant, function(id) {
    neighbors[[id]][!neighbors[[id]] %in% c(filtered_malignant, normal_ids)]
  }), use.names = FALSE))
  initial_rows <- lapply(isolated_malignant, function(id) {
    level <- as.character(cluster[match(id, spot_ids)])
    candidates <- neighbors[[id]][neighbors[[id]] %in% initial_candidates]
    if (!length(candidates)) return(NULL)
    center <- malignant_centers[[level]]
    tumor_distance <- vapply(candidates, function(candidate) {
      sqrt(sum((embedding[candidate, ] - center)^2))
    }, numeric(1))
    normal_distance <- vapply(candidates, function(candidate) {
      sqrt(sum((embedding[candidate, ] - normal_center)^2))
    }, numeric(1))
    data.frame(
      spot_id = candidates,
      proposed = ifelse(tumor_distance < umap_malignant_ratio * normal_distance,
                        "Mal", "Bdy"),
      stringsAsFactors = FALSE
    )
  })
  initial_rows <- do.call(rbind, initial_rows)
  if (is.null(initial_rows) || !nrow(initial_rows)) {
    stop("initial Cottrazm frontier has no unassigned neighbor", call. = FALSE)
  }
  initial_by_spot <- split(initial_rows$proposed, initial_rows$spot_id)
  initial <- data.frame(
    spot_id = names(initial_by_spot),
    location = vapply(initial_by_spot, function(value) {
      if ("Mal" %in% value) "Mal" else "Bdy"
    }, character(1)),
    stringsAsFactors = FALSE
  )
  malignant_ids <- unique(c(filtered_malignant,
                            initial$spot_id[initial$location == "Mal"]))
  boundary_ids <- unique(initial$spot_id[initial$location == "Bdy"])
  frontier <- malignant_ids
  round_states <- list(initial = initial)

  for (round_index in seq_len(as.integer(maximum_rounds))) {
    if (length(frontier) < 3L) break
    unseen <- unique(unlist(lapply(frontier, function(id) {
      neighbors[[id]][!neighbors[[id]] %in% c(malignant_ids, normal_ids, boundary_ids)]
    }), use.names = FALSE))
    if (length(unseen) < 3L) break
    update <- ct_frontier_update(
      frontier, neighbors, embedding, normal_ids, boundary_ids, malignant_ids,
      round_index, expand_malignant_radius
    )
    if (!nrow(update)) break
    round_states[[paste0("round_", round_index)]] <- update
    new_malignant <- update$spot_id[startsWith(update$location, "Mal")]
    new_boundary <- update$spot_id[update$location == "Bdy"]
    malignant_ids <- unique(c(malignant_ids, new_malignant))
    boundary_ids <- unique(c(boundary_ids, new_boundary))
    frontier <- unique(new_malignant)
  }

  normal_boundary <- unique(unlist(lapply(malignant_ids, function(id) {
    neighbors[[id]][!neighbors[[id]] %in% c(boundary_ids, malignant_ids)]
  }), use.names = FALSE))
  boundary_ids <- unique(c(boundary_ids, normal_boundary))
  location <- rep("nMal", n)
  names(location) <- spot_ids
  location[boundary_ids] <- "Bdy"
  location[malignant_ids] <- "Mal"

  list(
    location = unname(location[spot_ids]),
    selected_malignant_labels = malignant_labels,
    normal_cluster = normal_cluster,
    malignant_ids = malignant_ids,
    boundary_ids = boundary_ids,
    normal_ids = normal_ids,
    round_states = round_states
  )
}
