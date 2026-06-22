# recon_helpers.R：08_spatial_reconstruction.R 的表达重构 helper。
# 目标是把一个 spot 的 gene count 按反卷积比例和 signature 权重拆成多个 celltype_spot 列。

get_feature_weight <- function(filter_sig) {
  # 对每个基因，把不同 celltype 的 signature 表达归一化成权重。
  W_st <- do.call(rbind, lapply(rownames(filter_sig), function(feature) {
    a <- filter_sig[feature, ]
    a / sum(a)
  }))
  rownames(W_st) <- rownames(filter_sig)
  W_st
}

get_obs <- function(DeconData, obs_ID) {
  # 为每个目标 spot 取出比例 > 0 的细胞类型及其比例。
  tibble::tibble(cell_ID = DeconData[DeconData$cell_ID %in% obs_ID, ]$cell_ID) |>
    dplyr::mutate(Decon = purrr::map(.x = cell_ID, .f = function(.x) {
      sub_Decon <- DeconData[DeconData$cell_ID == .x, 2:ncol(DeconData), drop = FALSE]
      sub_DeconF <- as.data.frame(sub_Decon[, colnames(sub_Decon)[which(sub_Decon > 0)], drop = FALSE])
      if (ncol(sub_DeconF) == 1) colnames(sub_DeconF) <- colnames(sub_Decon)[which(sub_Decon > 0)]
      sub_DeconF
    }))
}

get_recon_mtx <- function(TumorST, sig_exp, clustermarkers_list, DeconData, Location) {
  # 主重构函数：
  # 1) 选 Location 中的 spot；
  # 2) 选 marker/signature/空间表达交集基因；
  # 3) 对每个 spot 的每个基因按 signature 权重和反卷积比例分配 count。
  start_time <- Sys.time()
  SubID <- rownames(TumorST@meta.data[TumorST@meta.data$Location %in% Location, ])
  expr_values <- as.matrix(Seurat::GetAssayData(TumorST, assay = "Spatial", layer = "counts"))
  sub_expr <- expr_values[, SubID, drop = FALSE]
  ClusterMarkers <- unlist(clustermarkers_list)[grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", unlist(clustermarkers_list), invert = TRUE)]
  intersect_genes <- intersect(rownames(expr_values), unique(ClusterMarkers))
  filter_sig <- sig_exp[intersect(rownames(sig_exp), intersect_genes), , drop = FALSE]
  sub_expr_filter <- sub_expr[intersect(rownames(sig_exp), intersect_genes), , drop = FALSE]
  W_st <- get_feature_weight(filter_sig = filter_sig)
  obs <- get_obs(DeconData = DeconData, obs_ID = SubID)

  mtx <- c()
  pb <- txtProgressBar(style = 3)
  for (i in seq_along(obs$cell_ID)) {
    spot <- obs$cell_ID[i]
    decon <- as.data.frame(t(unlist(obs[obs$cell_ID == spot, ]$Decon)))
    sub_matrix <- matrix(0, nrow = nrow(sub_expr_filter), ncol = ncol(decon))
    colnames(sub_matrix) <- paste(colnames(decon), spot, sep = "_")
    rownames(sub_matrix) <- rownames(sub_expr_filter)
    for (feature in rownames(sub_matrix)) {
      # deno 是该基因在当前 spot 的 signature 加权期望；
      # deno 为 0 时无法按 signature 分配，退回复制该 spot 的原始 count。
      sub_wt <- W_st[feature, colnames(decon)]
      deno <- sum(sub_wt * decon)
      feature_modi <- if (deno == 0) rep(sub_expr_filter[feature, spot], ncol(sub_matrix)) else sub_expr_filter[feature, spot] / deno * sub_wt
      sub_matrix[feature, ] <- feature_modi
    }
    mtx <- cbind(mtx, sub_matrix)
    setTxtProgressBar(pb, i / length(obs$cell_ID))
  }
  close(pb)
  message(sprintf("Time to get reconstructed matrix was %s hours", Sys.time() - start_time))
  mtx
}
