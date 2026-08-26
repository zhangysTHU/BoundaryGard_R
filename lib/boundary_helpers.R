# boundary_helpers.R：05_define_boundary.R 的空间邻接和边界扩展 helper。
# 输入通常是 Space Ranger 坐标、UMAP 坐标和已标注的 Mal/Bdy/Normal spot 集合。

nbrs <- function(df_j, MalCellIDAdd, CellIDRaw) {
  # 对一组恶性边缘 spot，找出尚未出现在 CellIDRaw 中的空间邻居，作为下一轮候选。
  nbrs_of_Mal <- lapply(MalCellIDAdd, function(id) {
    nbs <- df_j[[id]]
    nbs[!nbs %in% CellIDRaw]
  })
  names(nbrs_of_Mal) <- MalCellIDAdd
  nbrs_of_Mal
}

load_spaceranger_positions <- function(spaceranger_dir, cells) {
  # 兼容 Space Ranger 新旧坐标文件：
  # tissue_positions.csv 有 header，tissue_positions_list.csv 通常无 header。
  # 统一输出 row/col/imagerow/imagecol，并按 Seurat 对象中的 cells 排序。
  spatial_dir <- file.path(spaceranger_dir, "spatial")
  pos_file <- file.path(spatial_dir, "tissue_positions.csv")
  if (!file.exists(pos_file)) {
    pos_file <- file.path(spatial_dir, "tissue_positions_list.csv")
  }
  if (!file.exists(pos_file)) {
    stop("Cannot find tissue_positions.csv or tissue_positions_list.csv in ", spatial_dir, call. = FALSE)
  }

  first_line <- readLines(pos_file, n = 1, warn = FALSE)
  has_header <- grepl("barcode|array_row|pxl_row", first_line, ignore.case = TRUE)
  pos_raw <- utils::read.csv(pos_file, header = has_header, stringsAsFactors = FALSE)
  if (!has_header) {
    colnames(pos_raw) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  }

  required <- c("barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  missing_cols <- setdiff(required, colnames(pos_raw))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", pos_file, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  pos <- data.frame(
    row = pos_raw$array_row,
    col = pos_raw$array_col,
    imagerow = pos_raw$pxl_row_in_fullres,
    imagecol = pos_raw$pxl_col_in_fullres,
    row.names = pos_raw$barcode
  )
  pos <- pos[cells, , drop = FALSE]
  if (any(!stats::complete.cases(pos))) {
    missing_cells <- cells[!stats::complete.cases(pos)]
    stop("Missing spaceranger coordinates for ", length(missing_cells), " cells; first missing: ", missing_cells[[1]], call. = FALSE)
  }
  pos
}

compute_interspot_distances <- function(position, scale.factor = 1.05) {
  # 用阵列 row/col 与图像像素 imagerow/imagecol 的线性关系估计相邻 spot 半径。
  cols <- c("row", "col", "imagerow", "imagecol")
  assertthat::assert_that(all(cols %in% colnames(position)))
  list(
    xdist = coef(lm(position$imagecol ~ position$col))[2],
    ydist = coef(lm(position$imagerow ~ position$row))[2]
  ) |>
    within(radius <- (abs(xdist) + abs(ydist)) * scale.factor)
}

find_neighbors <- function(position, radius, method = c("manhattan", "euclidean")) {
  # 根据像素坐标计算所有 spot 的邻居表；返回 named list，名字是 spot barcode。
  method <- match.arg(method)
  pdist <- as.matrix(stats::dist(as.matrix(position[, c("imagecol", "imagerow")]), method = method))
  neighbors <- (pdist <= radius & pdist > 0)
  df_j2 <- sapply(seq_len(nrow(position)), function(x) as.vector(names(neighbors[x, ])[which(neighbors[x, ])]))
  names(df_j2) <- rownames(position)
  df_j2
}

ClusterUpdate <- function(MalCellIDN, position, df_j, UMAPembeddings, NormalCellID, BdyCellID, MalCellID, x,
                          expand_mal_radius = 0.8) {
  # 第 x 轮边界扩展：
  # 对上一轮新 Mal spot 的未标注邻居，比较其到 Mal 中心和 Bdy 中心的 UMAP 距离，
  # 决定标为 Malx 还是 Bdy。
  P <- lapply(MalCellIDN, function(i) {
    ncellID <- df_j[[i]]
    cpos <- UMAPembeddings[i, ]
    npos <- UMAPembeddings[ncellID, , drop = FALSE]
    nMalID <- ncellID[ncellID %in% MalCellID]
    CiMal <- data.frame(t(apply(rbind(npos[nMalID, , drop = FALSE], cpos), 2, mean)))
    rMal <- lapply(c(nMalID, i), function(id) sqrt(sum((rbind(npos, cpos)[id, ] - CiMal)^2))) |> unlist()
    nBdyID <- ncellID[ncellID %in% BdyCellID]
    nbrsID <- ncellID[!ncellID %in% MalCellID & !ncellID %in% BdyCellID & !ncellID %in% NormalCellID]

    if (length(nbrsID) == 0) {
      return(data.frame(cellID = character(), p1 = numeric(), p2 = numeric(), cluster = character()))
    }

    if (length(nBdyID) <= 1) {
      p <- lapply(nbrsID, function(id) sqrt(sum((npos[id, ] - CiMal)^2)))
      d <- data.frame(cellID = nbrsID, p1 = unlist(p), p2 = unlist(p))
      d$cluster <- ifelse(d$p1 <= expand_mal_radius * max(rMal), "Mal", "Bdy")
    } else {
      CiBdy <- data.frame(t(apply(npos[nBdyID, , drop = FALSE], 2, mean)))
      rBdy <- lapply(nBdyID, function(id) sqrt(sum((npos[id, ] - CiBdy)^2))) |> unlist()
      p1 <- lapply(nbrsID, function(id) sqrt(sum((npos[id, ] - CiMal)^2)))
      p2 <- lapply(nbrsID, function(id) sqrt(sum((npos[id, ] - CiBdy)^2)))
      d <- data.frame(cellID = nbrsID, p1 = unlist(p1), p2 = unlist(p2))
      d$cluster <- ifelse(d$p1 <= expand_mal_radius * max(rMal) & d$p2 > max(rBdy), "Mal", "Bdy")
    }
    d
  })

  PDF <- do.call(rbind, lapply(P, data.frame))
  if (is.null(PDF) || nrow(PDF) == 0) {
    return(data.frame(CellID = character(), Location = character()))
  }
  rownames(PDF) <- NULL
  Cluster <- as.data.frame.array(table(PDF$cellID, PDF$cluster))
  if ("Mal" %in% colnames(Cluster)) {
    Cluster$Location <- ifelse(Cluster$Mal > 0, paste("Mal", x, sep = ""), "Bdy")
  } else {
    Cluster$Location <- rep("Bdy", nrow(Cluster))
  }
  data.frame(CellID = rownames(Cluster), Location = as.character(Cluster$Location))
}
