# 09：按空间区域 Mal / Bdy / nMal 做差异表达和富集。
# 输入：intermediate/05_TumorST_boundary_defined.rds.gz，使用 params$diff_assay 的 normalized data。
# 输出：
# - intermediate/09_DiffGenes.rds.gz：命名 list，每个 Location 一个差异表，供 10 火山图使用。
# - intermediate/09_Enrichment.rds.gz：GO/KEGG 富集结果。
# - output/09_diff_enrichment/DiffGenes_<Location>.xlsx：每个基因的 Diff、pvalue、Symbol、FDR。
# - output/09_diff_enrichment/DiffGenesSig.rds.gz：上调且显著的差异基因子集。
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "00_config.R"))
load_required_packages(c("Seurat", "magrittr", "tibble", "dplyr", "clusterProfiler", "org.Hs.eg.db", "readr", "openxlsx"))

out_dir <- file.path(paths$output, "09_diff_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TumorST <- readr::read_rds(file.path(paths$intermediate, "05_TumorST_boundary_defined.rds.gz"))
assay <- params$diff_assay
if (assay == "Spatial") {
  TumorST <- Seurat::NormalizeData(TumorST, assay = assay)
}

# t.test 可能因为某组为空或方差问题失败；失败时返回 NA，后面会过滤。
my_t_test <- function(...) {
  obj <- try(t.test(...), silent = TRUE)
  if (inherits(obj, "try-error")) NA else obj$p.value
}

TumorSTm <- Seurat::GetAssayData(TumorST, assay = assay, layer = "data")
TumorSTm <- data.frame(TumorSTm) |> tibble::rownames_to_column(var = "gene")
colnames(TumorSTm) <- c("gene", as.character(TumorST$Location))
MajorTypes <- data.frame(
  MajorTypes = data.frame(do.call(rbind, strsplit(colnames(TumorSTm[, 2:ncol(TumorSTm)]), "\\.")))$X1,
  SampleID = colnames(TumorSTm[2:ncol(TumorSTm)])
)

# 对每个 Location，做该 Location vs 其它所有 spot 的逐基因差异。
DiffGenes <- lapply(split(MajorTypes[, "SampleID"], MajorTypes$MajorTypes), function(ID) {
  OtherSampleID <- MajorTypes$SampleID[!(MajorTypes$SampleID) %in% ID]
  DiffGenes <- apply(TumorSTm[, colnames(TumorSTm) %in% MajorTypes$SampleID], 1, function(x) {
    CanExp <- as.numeric(x[ID])
    CanExp <- CanExp[!is.na(CanExp)]
    OtherExp <- as.numeric(x[OtherSampleID])
    OtherExp <- OtherExp[!is.na(OtherExp)]
    c(Diff = mean(CanExp) - mean(OtherExp), pvalue = my_t_test(CanExp, OtherExp))
  }) |> t() |> data.frame()
  DiffGenes$Symbol <- TumorSTm$gene
  DiffGenes <- DiffGenes[!is.na(DiffGenes$pvalue), ]
  DiffGenes$FDR <- p.adjust(DiffGenes$pvalue, method = "fdr")
  DiffGenes
})

# 富集只用上调且 FDR 达标的基因，并过滤免疫球蛋白/核糖体/线粒体等常见干扰基因。
DiffGenesSig <- lapply(DiffGenes, function(sub) {
  sub <- sub[sub$Diff >= params$diff_logfc_cutoff & sub$FDR <= params$diff_fdr_cutoff, ]
  sub <- sub[sub$Symbol %in% grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", sub$Symbol, invert = TRUE, value = TRUE), ]
  sub[!is.na(sub$Diff), ]
})

locations <- intersect(c("Bdy", "Mal", "nMal"), names(DiffGenesSig))
run_enrich_go <- function(x) {
  if (length(x) == 0) return(NULL)
  clusterProfiler::enrichGO(gene = x, keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.2)
}
run_enrich_kegg <- function(x) {
  if (length(x) == 0) return(NULL)
  geneL <- clusterProfiler::bitr(x, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Hs.eg.db")
  if (nrow(geneL) == 0) return(NULL)
  obj <- try(
    clusterProfiler::enrichKEGG(
      gene = geneL$ENTREZID,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      minGSSize = 3,
      pvalueCutoff = 0.2,
      qvalueCutoff = 0.2
    ),
    silent = TRUE
  )
  if (inherits(obj, "try-error")) {
    warning("KEGG enrichment skipped because KEGG REST was unavailable.", call. = FALSE)
    return(NULL)
  }
  obj
}
Enrichment <- tibble::tibble(Location = locations) |>
  dplyr::mutate(LocationDiffFeatures = purrr::map(Location, function(x) DiffGenesSig[[x]]$Symbol)) |>
  dplyr::mutate(GO = purrr::map(LocationDiffFeatures, run_enrich_go)) |>
  dplyr::mutate(KEGG = purrr::map(LocationDiffFeatures, run_enrich_kegg))

readr::write_rds(DiffGenes, file.path(paths$intermediate, "09_DiffGenes.rds.gz"), compress = "gz")
readr::write_rds(Enrichment, file.path(paths$intermediate, "09_Enrichment.rds.gz"), compress = "gz")
for (nm in names(DiffGenes)) {
  openxlsx::write.xlsx(DiffGenes[[nm]], file.path(out_dir, paste0("DiffGenes_", nm, ".xlsx")), overwrite = TRUE)
}
readr::write_rds(DiffGenesSig, file.path(out_dir, "DiffGenesSig.rds.gz"), compress = "gz")
