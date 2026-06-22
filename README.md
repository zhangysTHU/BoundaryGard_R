# BoundaryGrad R workflow

本仓库是 Cottrazm/BoundaryGrad 的 R 顺序脚本工作流。默认样本名为 `CRC1`，输入、输出和中间文件均按样本名分目录保存。R 版只有一条 inferCNV 路线，因此不使用 03a/03b 后缀。

## 目录放置规则

默认目录如下：

```text
BoundaryGard_R/
  input/
    CRC1/
      spaceranger_outs/
      single_cell/
  intermediate/
    CRC1/
  output/
    CRC1/
```

`input/CRC1/` 是样本输入目录。`intermediate/CRC1/` 是脚本之间传递的 RDS、TSV 等中间文件目录。`output/CRC1/` 是图、表、inferCNV 输出和人工检查结果目录。

仓库完整保留 `input/`，并上传 `output/` 中小于 GitHub 100 MB 单文件限制的结果。大型 inferCNV 对象、SME 矩阵、`intermediate/` 与本地 renv 缓存不纳入 Git；详见 `UPLOAD_EXCLUSIONS.md`。

## 空间转录组输入

Space Ranger 结果必须放在：

```text
input/CRC1/spaceranger_outs/
```

该目录至少需要包含以下内容之一：

```text
filtered_feature_bc_matrix.h5
```

或矩阵目录：

```text
filtered_feature_bc_matrix/
  matrix.mtx.gz
  features.tsv.gz
  barcodes.tsv.gz
```

同时必须包含空间图像和坐标目录：

```text
spatial/
  tissue_positions_list.csv
  scalefactors_json.json
  tissue_lowres_image.png
  tissue_hires_image.png
```

这些文件会被 `01_preprocess_st.R` 读取并转换为 Seurat 对象。样本名前缀来自 `sample_name`，默认是 `CRC1`。

## 单细胞参考输入

单细胞参考文件必须放在：

```text
input/CRC1/single_cell/
```

推荐直接提供整理好的 Cottrazm/vignette 风格参考文件：

```text
sig_exp.rds.gz
clustermarkers_list.rds.gz
```

`sig_exp.rds.gz` 要求：

- R 对象应为 matrix 或 data.frame。
- 行名是基因名。
- 列名是细胞类型名。
- 数值是该细胞类型的平均表达或 signature 表达。

`clustermarkers_list.rds.gz` 要求：

- R 对象应为命名 list。
- list 的名字是细胞类型名。
- 每个元素是该细胞类型的 marker gene 字符向量。

如果没有预先整理好的 signature，也可以提供：

```text
single_cell_seurat.rds
clustermarkers_list.rds.gz
define_types.txt
```

其中 `single_cell_seurat.rds` 应为 Seurat 对象；`define_types.txt` 用于指定或映射细胞类型。`06_prepare_single_cell_reference.R` 会据此生成：

```text
intermediate/CRC1/06_sig_exp.rds.gz
intermediate/CRC1/06_clustermarkers_list.rds.gz
```

## 运行方式

先在 R 中执行 `renv::restore()` 恢复依赖。PowerShell 脚本默认调用 PATH 中的 `Rscript`，也可通过 `-Rscript` 指定完整路径。

完整运行默认样本 `CRC1`：

```powershell
.\run_all.ps1
```

指定其他样本名：

```powershell
.\run_all.ps1 -SampleName CRC2
```

此时输入应放在：

```text
input/CRC2/
```

输出会写入：

```text
output/CRC2/
intermediate/CRC2/
```

也可以分步骤运行，例如只跑 1 到 5 步：

```powershell
.\run_all.ps1 -FromStep 1 -ToStep 5
```

## 关键输出

R 版 inferCNV 会生成统一 CNV 结果表：

```text
intermediate/<样本名>/03_cnv_calls.tsv
```

该文件包含：

```text
cell_ID    CNVLabel    cnv_score
```

04 之后的脚本主要依赖 `intermediate/<样本名>/` 中的中间对象和表格。

## 常见检查点

- 输入目录名必须和运行样本名一致，例如 `CRC1` 对应 `input/CRC1/`。
- 新样本请新建 `input/<新样本名>/`，不要覆盖已有样本。
- `output/<样本名>/` 和 `intermediate/<样本名>/` 可以删除后重跑。
- `input/<样本名>/` 是原始输入目录，不应被流程脚本写入或清空。
