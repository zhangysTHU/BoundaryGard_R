# decon_helpers.R：07_spatial_deconvolution.R 的 marker enrichment 和 DWLS 反卷积 helper。
# 主要数据方向：filter_expr / filter_log_expr 是 gene x spot，filter_sig 是 gene x celltype。

get_enrich_matrix <- function(filter_sig, clustermarkers_list) {
  # 构造 gene x celltype 的 0/1 marker 矩阵：1 表示该基因是该细胞类型 marker。
  enrich_matrix <- matrix(0, nrow = nrow(filter_sig), ncol = ncol(filter_sig))
  rownames(enrich_matrix) <- rownames(filter_sig)
  colnames(enrich_matrix) <- colnames(filter_sig)
  for (cluster in colnames(enrich_matrix)) {
    feature <- intersect(clustermarkers_list[[cluster]], rownames(enrich_matrix))
    enrich_matrix[feature, cluster] <- 1
  }
  enrich_matrix
}

enrich_analysis <- function(filter_log_expr, enrich_matrix) {
  # 类似 z-score 的 marker enrichment：
  # 比较某细胞类型 marker 的平均 fold-change 是否高于该 spot 的背景基因变化。
  mean_gene_expr <- log2(rowMeans(2^filter_log_expr - 1, dims = 1) + 1)
  geneFold <- filter_log_expr - mean_gene_expr
  cellColMean <- apply(geneFold, 2, mean)
  cellColSd <- apply(geneFold, 2, stats::sd)
  enrichment <- matrix(data = NA, nrow = ncol(enrich_matrix), ncol = length(cellColMean))
  for (i in seq_len(ncol(enrich_matrix))) {
    signames <- rownames(enrich_matrix)[which(enrich_matrix[, i] == 1)]
    sigColMean <- apply(geneFold[signames, , drop = FALSE], 2, mean)
    m <- length(signames)
    enrichment[i, ] <- (sigColMean - cellColMean) * sqrt(m) / cellColSd
  }
  rownames(enrichment) <- colnames(enrich_matrix)
  colnames(enrichment) <- names(cellColMean)
  enrichment
}

solve_OLS_internal <- function(S, B) {
  # 非负约束 OLS 初始解：Signature 矩阵 S 乘比例向量，拟合 spot 表达 B。
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(ncol(S)))
  bzero <- rep(0, ncol(S))
  sc <- norm(D, "2")
  pd_D_mat <- Matrix::nearPD(D / sc)
  solution <- quadprog::solve.QP(as.matrix(pd_D_mat$mat), d / sc, A, bzero)$solution
  names(solution) <- colnames(S)
  solution
}

find_dampening_constant <- function(S, B, goldStandard) {
  # DWLS 中选择 dampening 常数：通过重复抽样评估解的稳定性，取方差最小的 j。
  solutionsSd <- NULL
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsScaledMinusInf <- wsScaled
  if (is.infinite(max(wsScaled))) {
    wsScaledMinusInf <- wsScaled[-which(is.infinite(wsScaled))]
  }
  for (j in seq_len(ceiling(log2(max(wsScaledMinusInf))))) {
    multiplier <- 2^(j - 1)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    solutions <- NULL
    for (i in 1:100) {
      set.seed(i)
      subset <- sample(length(ws), size = length(ws) * 0.5)
      fit <- stats::lm(B[subset] ~ -1 + S[subset, ], weights = wsDampened[subset])
      sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
      solutions <- cbind(solutions, sol)
    }
    solutionsSd <- cbind(solutionsSd, apply(solutions, 1, stats::sd))
  }
  which.min(colMeans(solutionsSd^2))
}

solve_dampened_WLSj <- function(S, B, goldStandard, j) {
  # 给高权重残差加上 dampening 上限，避免少数基因主导加权最小二乘。
  multiplier <- 2^(j - 1)
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsDampened <- wsScaled
  wsDampened[which(wsScaled > multiplier)] <- multiplier
  W <- diag(wsDampened)
  D <- t(S) %*% W %*% S
  d <- t(S) %*% W %*% B
  A <- cbind(diag(ncol(S)))
  bzero <- rep(0, ncol(S))
  sc <- norm(D, "2")
  pd_D_mat <- Matrix::nearPD(D / sc)
  solution <- quadprog::solve.QP(as.matrix(pd_D_mat$mat), d / sc, A, bzero)$solution
  names(solution) <- colnames(S)
  solution
}

optimize_solveDampenedWLS <- function(S, B, constant_J) {
  # 迭代 DWLS：反复求 dampened WLS，并和上一轮解做平滑，直到收敛或达到 1000 次。
  solution <- solve_OLS_internal(S, B)
  iterations <- 0
  change <- 1
  while (change > .01 && iterations < 1000) {
    newsolution <- solve_dampened_WLSj(S, B, solution, constant_J)
    solutionAverage <- rowMeans(cbind(newsolution, matrix(solution, nrow = length(solution), ncol = 4)))
    change <- norm(Matrix::as.matrix(solutionAverage - solution))
    solution <- solutionAverage
    iterations <- iterations + 1
  }
  solution / sum(solution)
}

optimize_deconvolute_dwls <- function(exp, Signature) {
  # 对一个 topic 内的多个 spot 执行 DWLS，返回 celltype x spot 的比例矩阵。
  Genes <- intersect(rownames(Signature), rownames(as.matrix(exp)))
  S <- Matrix::as.matrix(Signature[Genes, , drop = FALSE])
  subBulk <- as.matrix(exp[Genes, , drop = FALSE])
  allCounts_DWLS <- NULL
  all_exp <- rowMeans(as.matrix(exp))
  solution_all_exp <- solve_OLS_internal(S, all_exp[Genes])
  constant_J <- find_dampening_constant(S, all_exp[Genes], solution_all_exp)
  for (j in seq_len(ncol(subBulk))) {
    B <- subBulk[, j]
    solDWLS <- optimize_solveDampenedWLS(S, B, constant_J)
    allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
  }
  colnames(allCounts_DWLS) <- colnames(exp)
  allCounts_DWLS
}

spot_proportion_initial <- function(enrich_matrix, enrich_result, filter_expr, filter_sig,
                                    clustermarkers_list, meta_data,
                                    malignant_cluster, tissue_cluster, stromal_cluster) {
  # 第一阶段反卷积：
  # 先按 topic 的 marker score/enrichment 选择候选细胞类型，再对该 topic 内 spot 做 DWLS。
  topic <- meta_data$Decon_topics
  topic_sort <- sort(unique(topic))
  ct_sig0 <- lapply(names(clustermarkers_list), function(cluster) {
    cutoff_sig <- quantile(meta_data[, cluster], 0.75, na.rm = TRUE)
    test <- vapply(topic_sort, function(topic_i) median(meta_data[meta_data$Decon_topics == topic_i, cluster], na.rm = TRUE), numeric(1))
    as.character(topic_sort[which(test - cutoff_sig > 0)])
  })
  names(ct_sig0) <- names(clustermarkers_list)

  cutoff_enrich <- vapply(rownames(enrich_result), function(cluster) {
    cluster_enrich <- enrich_result[cluster, ]
    quantile(cluster_enrich[cluster_enrich > 0], 0.75, na.rm = TRUE)
  }, numeric(1))

  dwls_results <- matrix(0, nrow = ncol(enrich_matrix), ncol = ncol(filter_expr))
  rownames(dwls_results) <- colnames(enrich_matrix)
  colnames(dwls_results) <- colnames(filter_expr)

  for (i in seq_along(topic_sort)) {
    ct_sig <- names(clustermarkers_list)[vapply(names(clustermarkers_list), function(cluster) topic_sort[i] %in% ct_sig0[[cluster]], logical(1))]
    cluster_i_enrich <- as.matrix(enrich_result[, which(topic == topic_sort[i]), drop = FALSE])
    row_i_max <- Rfast::rowMaxs(cluster_i_enrich, value = TRUE)
    names(row_i_max) <- rownames(cluster_i_enrich)
    ct_sort <- sort(row_i_max[which(row_i_max - cutoff_enrich > 0)], decreasing = TRUE)[2]
    ct_enrich <- names(row_i_max[which(row_i_max - cutoff_enrich > 0)])[which(row_i_max[which(row_i_max - cutoff_enrich > 0)] >= ct_sort)]

    if (malignant_cluster %in% ct_sig && tissue_cluster %in% ct_sig) {
      ct_sig <- ct_sig[ct_sig != malignant_cluster]
    }
    ct <- if (length(ct_sig) == 0) ct_enrich else unique(c(ct_sig, ct_enrich))
    if (length(ct) == 0) ct <- names(sort(row_i_max - cutoff_enrich, decreasing = TRUE))[1:2]
    if (malignant_cluster %in% ct || tissue_cluster %in% ct) ct <- unique(c(ct, names(sort(row_i_max - cutoff_enrich, decreasing = TRUE)[1:3])))
    if (strsplit(topic_sort[i], "_")[[1]][1] == "Mal") ct <- unique(c(ct, malignant_cluster))
    if (median(meta_data[meta_data$Decon_topics == topic_sort[i], ]$nCount_Spatial, na.rm = TRUE) < 5000) ct <- unique(c(ct, stromal_cluster))
    ct <- as.character(stats::na.omit(ct))
    ct <- intersect(ct, colnames(enrich_matrix))
    if (length(ct) == 0) next
    message(topic_sort[i], ": ", paste(ct, collapse = ", "))

    ct_gene <- unique(unlist(lapply(ct, function(ct_i) rownames(enrich_matrix)[which(enrich_matrix[, ct_i] == 1)])))
    uniq_ct_gene <- intersect(unique(ct_gene), rownames(filter_expr))
    select_enrichig_exp <- as.matrix(filter_sig[uniq_ct_gene, ct, drop = FALSE])
    cluster_i_cell <- which(topic == topic_sort[i])
    cluster_cell_exp <- as.matrix(filter_expr[uniq_ct_gene, cluster_i_cell, drop = FALSE])
    cluster_i_dwls <- optimize_deconvolute_dwls(cluster_cell_exp, select_enrichig_exp)
    dwls_results[ct, cluster_i_cell] <- cluster_i_dwls
  }
  dwls_results[dwls_results < 0] <- 0
  dwls_results
}

spot_deconvolution <- function(expr, meta_data, ct_exp, enrich_matrix, binary_matrix) {
  # 第二阶段反卷积：
  # 用第一阶段比例 >= 0.01 的细胞类型作为每个 spot 的候选集合，再逐 spot 精修比例。
  topic <- meta_data$Decon_topics
  topic_sort <- sort(unique(topic))
  dwls_results <- matrix(0, nrow = ncol(ct_exp), ncol = ncol(expr))
  rownames(dwls_results) <- colnames(ct_exp)
  colnames(dwls_results) <- colnames(expr)
  for (topic_i in topic_sort) {
    cluster_i_matrix <- as.matrix(binary_matrix[, which(topic == topic_i), drop = FALSE])
    row_i_max <- Rfast::rowMaxs(cluster_i_matrix, value = TRUE)
    ct_i <- rownames(cluster_i_matrix)[which(row_i_max == 1)]
    ct_i <- intersect(ct_i, colnames(ct_exp))
    if (length(ct_i) == 0) next
    if (length(ct_i) == 1) {
      dwls_results[ct_i[1], which(topic == topic_i)] <- 1
    } else {
      ct_gene <- unique(unlist(lapply(ct_i, function(ct) rownames(enrich_matrix)[which(enrich_matrix[, ct] == 1)])))
      uniq_ct_gene <- intersect(rownames(expr), unique(ct_gene))
      select_enrichig_exp <- ct_exp[uniq_ct_gene, ct_i, drop = FALSE]
      cluster_i_cell <- which(topic == topic_i)
      cluster_cell_exp <- as.matrix(expr[uniq_ct_gene, cluster_i_cell, drop = FALSE])
      colnames(cluster_cell_exp) <- colnames(expr)[cluster_i_cell]
      all_exp <- rowMeans(cluster_cell_exp)
      solution_all_exp <- solve_OLS_internal(select_enrichig_exp, all_exp)
      constant_J <- find_dampening_constant(select_enrichig_exp, all_exp, solution_all_exp)
      for (k in seq_len(ncol(cluster_cell_exp))) {
        B <- Matrix::as.matrix(cluster_cell_exp[, k])
        ct_enrichpot_k <- rownames(cluster_i_matrix)[which(cluster_i_matrix[, k] == 1)]
        if (length(ct_enrichpot_k) == 1) {
          dwls_results[ct_enrichpot_k[1], colnames(cluster_cell_exp)[k]] <- 1
        } else {
          ct_k_gene <- unique(unlist(lapply(ct_enrichpot_k, function(ct) rownames(enrich_matrix)[which(enrich_matrix[, ct] == 1)])))
          uniq_ct_k_gene <- intersect(rownames(ct_exp), unique(ct_k_gene))
          S_k <- Matrix::as.matrix(ct_exp[uniq_ct_k_gene, ct_enrichpot_k, drop = FALSE])
          solDWLS <- optimize_solveDampenedWLS(S_k, B[uniq_ct_k_gene, ], constant_J)
          dwls_results[names(solDWLS), colnames(cluster_cell_exp)[k]] <- solDWLS
        }
      }
    }
  }
  dwls_results[dwls_results < 0] <- 0
  dwls_results
}
