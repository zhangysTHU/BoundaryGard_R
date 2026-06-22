# 08：根据 07 的细胞比例和 06 的 signature，重构指定区域的“伪单细胞/亚型表达矩阵”。
# 输入：
# - intermediate/05_TumorST_boundary_defined.rds.gz：提供 Spatial counts 和 Location。
# - intermediate/06_sig_exp.rds.gz、06_clustermarkers_list.rds.gz：提供 signature 和重构基因集合。
# - intermediate/07_DeconData.rds.gz：提供每个 spot 的细胞类型比例。
# 输出：
# - intermediate/08_reconstructed_matrix.rds.gz：gene x reconstructed_cell 矩阵，列名形如 celltype_spotbarcode。
# - intermediate/08_TumorST_reconstructed.rds.gz：Seurat 对象，metadata 含 Subtypes、orig.ident、Location。
# 下游：当前 09/10 不依赖 08；08 主要供边界区重构表达单独分析。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "tibble", "dplyr", "purrr", "readr"))
source(file.path(paths$lib, "recon_helpers.R"))

out_dir <- file.path(paths$output, "08_spatial_reconstruction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TumorST <- readr::read_rds(file.path(paths$intermediate, "05_TumorST_boundary_defined.rds.gz"))
sig_exp <- readr::read_rds(file.path(paths$intermediate, "06_sig_exp.rds.gz"))
clustermarkers_list <- readr::read_rds(file.path(paths$intermediate, "06_clustermarkers_list.rds.gz"))
DeconData <- readr::read_rds(file.path(paths$intermediate, "07_DeconData.rds.gz"))

# smoke-test 或特殊数据中可能没有 Bdy；此时写空结果并退出，避免下游误以为重构成功。
selected_spots <- rownames(TumorST@meta.data)[TumorST@meta.data$Location %in% params$recon_locations]
if (length(selected_spots) == 0) {
  warning(
    "No spots found for params$recon_locations=",
    paste(params$recon_locations, collapse = ", "),
    "; writing empty reconstruction outputs."
  )
  readr::write_rds(matrix(numeric(0), nrow = 0, ncol = 0), file.path(paths$intermediate, "08_reconstructed_matrix.rds.gz"), compress = "gz")
  readr::write_rds(NULL, file.path(paths$intermediate, "08_TumorST_reconstructed.rds.gz"), compress = "gz")
  quit(save = "no", status = 0)
}

# 核心重构：把每个目标 spot 的 gene count 按 signature 权重和反卷积比例拆到不同细胞类型列。
mtx <- get_recon_mtx(
  TumorST = TumorST,
  sig_exp = sig_exp,
  clustermarkers_list = clustermarkers_list,
  DeconData = DeconData,
  Location = params$recon_locations
)

TumorSTRecon <- Seurat::CreateSeuratObject(counts = as.matrix(mtx), project = "TumorSTRecon", assay = "RNA")
split_ids <- as.data.frame(do.call(rbind, strsplit(rownames(TumorSTRecon@meta.data), split = "_")))
TumorSTRecon@meta.data$Subtypes <- gsub("\\.", "_", split_ids$V1)
TumorSTRecon@meta.data$orig.ident <- split_ids$V2
TumorSTRecon@meta.data$Location <- TumorST@meta.data$Location[match(TumorSTRecon@meta.data$orig.ident, rownames(TumorST@meta.data))]

readr::write_rds(mtx, file.path(paths$intermediate, "08_reconstructed_matrix.rds.gz"), compress = "gz")
readr::write_rds(TumorSTRecon, file.path(paths$intermediate, "08_TumorST_reconstructed.rds.gz"), compress = "gz")
