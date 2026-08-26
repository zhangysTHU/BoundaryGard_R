# Cottrazm R 顺序脚本版说明

本目录是 Cottrazm 的 R 顺序脚本工作流。默认样本名为 `CRC1`，输入、输出和中间文件均按样本名分目录保存。R 版只有一条 inferCNV 路线，因此不使用 03a/03b 后缀。

## 流水线功能与生物学应用

这套脚本面向 10x Visium / Space Ranger 空间转录组样本，核心目标是把肿瘤组织切片上的 spot 按恶性核心、肿瘤-非肿瘤交界区和非恶性区域分层，并进一步解释这些区域中的细胞组成、差异表达、通路活性和空间梯度方向。其他项目调用该流水线时，可以把它理解为一条从原始 Visium 输出到“肿瘤边界生态位解析”的标准流程，适合回答肿瘤侵袭边缘、免疫排斥/浸润、基质屏障、血管/淋巴重塑、缺氧与 EMT 场、肿瘤-正常组织过渡状态等问题。

完整步骤如下：

- `01_preprocess_st.R` 读取 Space Ranger 表达矩阵、spot 坐标和 H&E 图像，构建 Seurat 空间对象并输出基础 QC 指标。
- `02_morphology_adjusted_cluster.R` 使用 stLearn/SME 结合表达和 H&E 形态信息生成 `Morph` assay，再做 PCA、UMAP 和聚类；同时用免疫/B 细胞 marker 计算 `NormalScore`，供 inferCNV 选择参考区域。
- `03_run_infercnv.R` 在空间 spot 上运行 inferCNV，基于正常参考 spot 推断 CNV 状态，并生成统一的 `03_cnv_calls.tsv`。
- `04_score_cnv.R` 将 CNV label 和 `cnv_score` 合并回 Seurat 对象，输出 CNV 的空间图、UMAP 图和分布图。
- `05_define_boundary.R` 是边界识别的关键步骤。脚本默认从 CNV score 中位数最高的有效 observation label 自动选取恶性 CNV 种子，也可以通过 `COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS` 手动指定。随后结合三类信息划分空间区域：CNV 恶性信号定义肿瘤种子；形态校正 UMAP 中恶性中心与正常中心的距离决定邻近 spot 更像恶性还是边界；真实空间邻接关系限制边界只能沿组织相邻 spot 逐圈扩展。最终在 `TumorST@meta.data$Location` 中写入 `Mal`、`Bdy`、`nMal` 三类标签，其中 `Mal` 表示恶性核心/扩展层，`Bdy` 表示肿瘤边缘及贴邻的正常侧过渡区，`nMal` 表示其余非恶性区域。该结果是 07、09、10、11、12 的共同分组基础。
- `06_prepare_single_cell_reference.R` 准备空间反卷积所需的单细胞 signature matrix 和 cell-type marker list。推荐直接提供项目适配后的 `sig_exp.rds.gz` 与 `clustermarkers_list.rds.gz`。
- `07_spatial_deconvolution.R` 使用单细胞参考对每个 spot 做细胞类型比例估计，输出 `07_DeconData.rds.gz` 和 `output/<样本名>/07/spot_matrix_pre_lsgi.csv`。这些结果可用于比较 `Mal/Bdy/nMal` 的细胞组成，也为 LSGI 梯度分析提供 cell-component embedding。
- `08_spatial_reconstruction.R` 按反卷积比例和 signature 将指定区域的 spot 表达拆解为伪单细胞/亚型表达矩阵。默认重构 `Bdy`，适合对边界区特定细胞类型做后续差异、通讯或轨迹分析。
- `09_diff_and_enrichment.R` 按 `Location` 做 `Mal`、`Bdy`、`nMal` 的区域差异表达与 GO/KEGG 富集，常用于识别边界特异上调基因、免疫反应、ECM、EMT、缺氧、细胞周期等生物过程。
- `10_plot_results.R` 汇总绘制区域细胞组成柱状图、H&E 上的 spot 反卷积 pie 图和区域差异基因火山图。
- `11_lsgi_gradient.R` 是空间梯度解析步骤。脚本先在组织空间上构建统一 grid，并为每个 grid 记录两类 spot 归属：nearest-grid partition 和局部回归使用的 local spots。随后对多种 embedding/component 来源运行 LSGI 局部线性梯度分析：`cell_component` 使用 07 的反卷积细胞类型比例；`nmf` 从空间表达中自动选择 NMF rank 并提取表达程序；`marker_module`、`pathway`、`single_gene` 从 `resources/lsgi_embedding_catalog/` 中读取预定义 marker 模块、MSigDB 通路或单基因。输出的 arrow 表记录每个 grid 上梯度所属 component、方向向量、R2、箭头长度和过滤状态，图形会把梯度箭头叠加到 `Mal/Bdy/nMal` 边界图或 H&E 图上。生物学上，这一步用于判断某类细胞比例、表达程序或通路活性在空间中从哪里增强、朝哪个方向变化，例如 CAF/ECM 屏障是否沿肿瘤边缘排列，IFN-γ/CXCL9/CXCL10 信号是否从免疫富集区指向肿瘤，缺氧/EMT/NMF 肿瘤程序是否由肿瘤核心向边界扩散。
- `12_boundary_related_lsgi_arrows.R` 将 11 的全局 LSGI arrows 与 05 定义的 `Bdy` spot 对齐，筛选“边界相关”的梯度箭头。脚本内置多种策略：按 nearest partition 中 Bdy spot 数量/比例筛选的 `partition_*`，按 local-regression spot 中 Bdy spot 数量/比例筛选的 `local_*`，以及推荐阈值的 `primary_consensus` 和 `primary_union`。每个 component method 都会输出选中 grid、选中 arrows、Bdy spot 到 arrow 的对应关系和可视化 PDF。该步骤适合把 11 中大量空间梯度收敛到肿瘤边界问题本身，例如筛选真正落在边界生态位的 CAF、巨噬细胞、T 细胞、ECM、抗原呈递、缺氧、血管生成或侵袭相关梯度。

跨项目复用时，建议至少完成三类适配。第一，检查 `00_config.R` 中样本名、聚类分辨率、inferCNV reference 参数、恶性 CNV label 策略和边界扩展参数；如果 CNV 自动选出的恶性标签与组织形态不符，应人工复核 04 的 CNV 图后设置 `COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS`。第二，确保单细胞参考的细胞类型命名与 `decon_malignant_cluster`、`decon_tissue_cluster`、`decon_stromal_cluster` 及 marker list 一致，否则 07 的反卷积和 11 的 cell-component 梯度会缺少关键成分。第三，根据项目假设修改 `resources/lsgi_embedding_catalog/marker_modules.tsv`、`pathways.tsv`、`single_genes.tsv` 或 `00_config.R` 中的 `lsgi_*_ids`，把 11/12 聚焦到当前癌种和生物学问题相关的程序、通路和基因。

调用下游结果时，最常用的入口文件是：

```text
intermediate/<样本名>/05_TumorST_boundary_defined.rds.gz
intermediate/<样本名>/07_DeconData.rds.gz
intermediate/<样本名>/09_DiffGenes.rds.gz
output/<样本名>/11_lsgi_gradient/component_method_summary.csv
output/<样本名>/11_lsgi_gradient/<component_method>/arrow_tables/arrows_by_grid.csv
output/<样本名>/12_boundary_related_lsgi_arrows/strategy_summary.csv
output/<样本名>/12_boundary_related_lsgi_arrows/<component_method>/arrow_tables/<strategy_id>_arrows.csv
output/<样本名>/12_boundary_related_lsgi_arrows/<component_method>/spot_tables/<strategy_id>_bdy_spot_arrow_assignments.csv
```

其中 `05_TumorST_boundary_defined.rds.gz` 是区域标签的权威对象；`07_DeconData.rds.gz` 是 spot 级细胞组成矩阵；`09_DiffGenes.rds.gz` 用于区域差异基因解释；`11_lsgi_gradient` 的 arrow 表用于全组织空间梯度分析；`12_boundary_related_lsgi_arrows` 的结果则用于只关注边界相关梯度。实际汇报时建议优先联合查看 `05_boundary` 的边界图、`10_plots` 的组成/差异图、`11_lsgi_gradient` 的全局梯度图和 `12_boundary_related_lsgi_arrows` 的边界筛选图，再决定哪些 arrows 或通路进入生物学验证。

## 目录放置规则

默认目录如下：

```text
scripts_format_R/
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

Linux 上处理样本时，默认方式是通过 `run_all_entry/runall_<sample>.sh` 以 `nohup` 后台提交，不再默认前台直跑。这样即使关闭终端，任务也会继续执行，日志统一写到 `run_all_entry/nohup_logs/`。

当前推荐入口：

```bash
bash run_all_entry/runall_CRC1.sh --keep-intermediate
bash run_all_entry/runall_BRCA1.sh --keep-intermediate
bash run_all_entry/runall_BRCA2.sh --keep-intermediate
```

这些 wrapper 默认就是后台模式；如果确实需要前台调试，再显式加 `--foreground`。

后台任务提交后，优先查看：

```text
run_all_entry/nohup_logs/<sample>_*.nohup.log
output/<样本名>/run_logs/
```

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

Linux 下如果不使用 wrapper，而是直接调用通用入口，则也建议显式放到后台：

```bash
nohup bash --noprofile --norc run_all_entry/run_all.sh --sample-name BRCA2 --keep-intermediate > run_all_entry/nohup_logs/BRCA2_manual.nohup.log 2>&1 < /dev/null &
```

## Linux 环境说明

当前项目默认使用统一的 Miniforge 环境 `BoundaryGrad`：

```text
/lulabdata3/huangkeyun/zhangys/tools/miniforge3/envs/BoundaryGrad
```

对应环境文件在项目根目录：

```text
../cottrazm_env_linux.yml
```

默认情况下不再自动启用项目内的 Windows `renv` 包库；若确实需要旧的 `renv` 路线，可显式设置：

```bash
export COTTRAZM_USE_RENV=1
```

正常在 Linux 上建议直接使用：

```bash
unset PYTHONHOME PYTHONPATH
export LANG=C
export LC_ALL=C
source /lulabdata3/huangkeyun/zhangys/tools/miniforge3/etc/profile.d/conda.sh
conda activate BoundaryGrad
```

如果后续从这个环境里直接启动 Jupyter，建议保留上面的 `unset PYTHONHOME PYTHONPATH` 与 `LANG/LC_ALL` 设置，避免宿主 shell 的 Python 或 locale 配置干扰 `ipykernel` / `IRkernel`。

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
