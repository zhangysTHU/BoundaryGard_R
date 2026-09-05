# Cottrazm R 整合模块版说明

本目录是 Cottrazm 的 R 顺序脚本工作流。默认样本名为 `CRC1`，输入、输出和中间文件均按样本名分目录保存。

当前默认流程已经从原来的 01-12 顺序脚本整理为 6 个主脚本。整合采用低影响方式：脚本数量减少，但关键中间文件名和下游输出目录保持旧路径，便于复用已有结果和下游脚本。

## 流水线功能与生物学应用

这套 R 工作流面向 10x Visium / Space Ranger 空间转录组切片，目标不是仅进行 spot 聚类或区域差异分析，而是把肿瘤组织解析为可计算的边界生态位。流程以原始空间表达矩阵、H&E 图像、Space Ranger 坐标和单细胞参考为输入，依次建立形态校正的空间对象、CNV 支持的肿瘤边界标签、spot 级细胞组成矩阵、多来源 LSGI 空间梯度方向场，以及以肿瘤上皮梯度为局部坐标轴的跨边界 profile。最终输出的核心对象是每个 spot 的 `Location = Mal / Bdy / nMal`、每个 grid 的 component-specific gradient arrow，以及每个边界 patch 上功能特征沿肿瘤方向轴的分布曲线。

该设计将 Cottrazm 的肿瘤边界定义思想与 LSGI 的局部空间梯度推断框架连接起来，并在二者之间加入 grid-aware boundary filtering 和 tumor-arrow-guided profile 两层下游分析。因而，使用者可以从三个层次解释肿瘤边界：第一，哪些 spot 构成恶性核心、边界和非恶性背景；第二，哪些细胞组分、转录程序或通路在边界附近形成方向性梯度；第三，这些功能特征在肿瘤上皮梯度指向的局部轴线上是边界富集、向肿瘤侧升高，还是向非恶性侧扩散。该框架适合研究肿瘤侵袭边缘、免疫浸润与免疫排斥、CAF/ECM 屏障、缺氧和 EMT 过渡、血管/淋巴重塑、脂质和 redox 保护生态位，以及治疗前后或不同疗效分组中的边界重塑。

### 0. 流程输入、配置和基础空间对象

`00_config.R` 不产生生物学结果，但固定了流程的可重复性和跨模块契约。它集中定义样本名、`input/intermediate/output/resources` 路径、inferCNV 参数、边界扩展阈值、反卷积细胞类型名称、LSGI embedding 来源、箭头可视化编码，以及第 6 模块使用的 tumor-arrow profile 参数。所有下游模块均从这里读取统一配置，因此同一套脚本可以通过环境变量或命令行参数切换样本、边界策略、feature panel 和 profile QC 标准。

`01_spatial_preprocess_cluster.R` 将 Space Ranger 表达矩阵、spot 坐标和 H&E 图像读入 Seurat，生成 `Spatial` assay 与 image slot，并输出空间 QC 图和表。随后脚本通过 stLearn/SME 生成融合局部组织形态信息的 `Morph` assay，在该形态校正表达空间中执行标准 Seurat 降维、邻居图、UMAP 和聚类。`NormalScore` 由免疫/B 细胞 marker 的平均表达近似估计，用于在下一步选择 CNV reference cluster。该模块的意义在于，后续边界定义不完全依赖原始表达聚类，而是同时继承了 H&E 形态和转录相似性的信息。

### 1. `02_boundary_definition.R`：CNV、形态表征和空间邻接共同约束的边界定义

`02_boundary_definition.R` 是整条流程的生物学坐标系构建模块。该模块的核心思想来源于 Cottrazm 对肿瘤恶性区和癌旁边界的定义策略，并在当前 R 版中被整合为 inferCNV 调用、CNV label 统一后处理、CNV score 写回、迭代边界扩展和空间可视化的一体化实现。其输入是 `01` 产生的形态校正 Seurat 对象，输出是带有 `CNVLabel`、`cnv_score` 和 `Location` 的 `05_TumorST_boundary_defined.rds.gz`。

算法首先以 `NormalScore` 最高的 Seurat cluster 作为 inferCNV reference group，对 `Spatial` assay 的原始计数矩阵运行 inferCNV，并从 HMM observation 表或 cell grouping/gene state 文件中重建统一的 `03_cnv_calls.tsv`。每个 spot 的 `cnv_score` 被定义为 HMM 状态相对于中性状态的偏离总量；若 dendrogram 可用，则根据树切分得到 `CNVLabel`，否则在 observation-level CNV score 上使用 k-means 生成 CNV label。为避免把 reference 或低质量过滤标签误当作恶性区域，脚本默认从有效 observation 标签中选择 `cnv_score` 中位数最高的两个标签作为 malignant CNV seeds，也允许用户通过 `boundary_malignant_cnv_labels` 手动指定经过病理或 marker 复核的标签。

在边界空间扩展阶段，脚本不把高 CNV spot 直接全部定义为肿瘤核心，而是进一步引入形态校正 UMAP 和 Visium 空间邻接关系。对于高 CNV 候选 spot，算法在每个恶性富集的 Seurat cluster 内计算恶性中心，同时计算正常参考中心；只有在 UMAP 空间中更接近恶性中心、并满足 `rt < boundary_umap_mal_ratio * rn` 的 spot 才作为初始 `Mal` seed。随后，脚本利用 Space Ranger 阵列坐标和像素坐标估计相邻 spot 半径，构建每个 spot 的空间邻居表，并从上一轮新增的 `Mal` spot 出发逐圈检查未标注邻居。每个候选邻居会根据其到局部恶性中心和边界中心的 UMAP 距离被更新为新的 `Mal` 扩展层或 `Bdy`，该过程最多迭代 `boundary_max_rounds` 轮。

最终标签被折叠为三类：`Mal` 表示 CNV 和形态表征共同支持的恶性核心及其扩展层；`Bdy` 表示与恶性区域空间邻接、但在 UMAP 表征中处于恶性与非恶性状态之间的边界 spot，并额外包含紧邻 `Mal` 的正常侧邻居；`nMal` 表示其余非恶性背景区域。因此，`Bdy` 不是简单的几何外圈，也不是单一 marker 阈值结果，而是由拷贝数异常、表达/形态相似性和空间连续性共同定义的过渡生态位。其生物学意义在于把后续分析锚定到真正发生肿瘤-正常组织接触、免疫进入或排斥、基质屏障形成和局部侵袭的空间区域。对应可视化包括每轮扩展图、最终 `BoundaryDefine.pdf` 以及 H&E 背景上的边界叠加图；结构化输出包括 `params.tsv`、`location_counts.tsv`、`04_TumorST_cnv_scored.rds.gz`、`05_TumorST_boundary_subset.rds.gz` 和 `05_TumorST_boundary_defined.rds.gz`。

### 2. 细胞组成与多来源 LSGI 梯度场

`03_spatial_deconvolution.R` 使用单细胞 signature matrix 和 marker list 对空间 spot 进行两阶段 DWLS 反卷积。第一阶段先按 `Location + seurat_clusters` 定义 topic，在 topic 层面结合 marker enrichment 和初始 DWLS 选择候选细胞类型；第二阶段再用初始比例二值化后的候选集合逐 spot 精修细胞组成比例。该模块输出 `07_DeconData.rds.gz` 和 `spot_matrix_pre_lsgi.csv`，为后续 `cell_component` 梯度提供 spot x cell type embedding。

`04_lsgi_gradient.R` 基于 LSGI 文献中的局部空间梯度推断思想，在组织切片上建立公共 meta-grid，并对每个 grid 取最近的 `n.cells.per.meta` 个 spot。对任意 `spot x component` embedding，LSGI 在每个局部邻域内拟合线性模型 `component ~ X + Y`，以回归系数 `vx` 和 `vy` 表示该 component 局部增加最快的方向，以 `sqrt(vx^2 + vy^2)` 表示局部梯度强度，以 `R_squared` 表示局部线性方向估计的可信度。当前 R 版将该算法扩展为多 embedding 来源框架：`cell_component` 来自反卷积比例，`nmf` 来自空间表达矩阵的 NMF 因子，`marker_module` 来自预设 marker module score，`pathway` 来自 MSigDB pathway rank score，`single_gene` 来自标准化单基因表达。所有方法共享同一套 grid、local spot membership 和 nearest-grid partition，使不同生物学层面的空间梯度可以在同一坐标系统中比较。

当前实现还对 LSGI 原始箭头可视化进行了二次编码。原始 LSGI 会把箭头缩放到近似等长，当前脚本保留原始方向向量，同时根据 `gradient_strength` 将箭头长度映射到相对旧默认长度的 0.5-2.0 倍；根据局部回归 `R_squared` 将箭头头部张角映射到 22.5-45 度。由此，箭头方向表示 feature 增加方向，箭头长度表示相对梯度强度，箭头头部张角表示拟合可信度，颜色表示 component 类型。输出包括 `grid_info.csv`、`grid_local_spot_membership.csv`、`grid_partition_membership.csv`、`arrow_tables/arrows_by_grid.csv`、`gradient_distance.csv`、`gradient_distance_heatmap.pdf`、不含背景的 `gradients_plain_lsgi.pdf`、Mal/Bdy/nMal 背景图和 H&E 叠加图。

### 3. `05_boundary_related_lsgi_arrows.R`：从全局方向场到边界相关梯度箭头

`05_boundary_related_lsgi_arrows.R` 是 Cottrazm 边界定义与 LSGI 局部梯度的耦合层。其输入不是原始表达矩阵，而是 `04` 输出的 LSGI arrow table、公共 grid membership 表和 `03` 输出的 spot 反卷积矩阵。该模块保留 LSGI 对局部方向场的定义，同时引入 Cottrazm 的 `Bdy` 标签，将全切片范围内的梯度箭头筛选为真正发生在肿瘤边界附近的方向性变化。

具体而言，脚本先分别在两种 grid-spot 关系上计算每个 grid 的边界负荷。`partition` 关系把每个 spot 分配给最近 grid，用于描述该 grid 管辖空间区域内 `Bdy` spot 的数量和比例；`local` 关系使用 LSGI 局部线性回归实际纳入的最近 spot，用于判断该箭头的回归邻域是否直接覆盖边界。基于这两类指标，模块内置 `partition_any`、`partition_relaxed`、`partition_primary`、`partition_strict`、`local_broad`、`local_relaxed`、`local_primary`、`local_strict`、`primary_consensus` 和 `primary_union` 等策略。`partition_*` 强调空间管辖区域中的边界富集，`local_*` 强调方向估计本身由边界邻域 spot 支持，`consensus/union` 分别提供更保守或更宽松的组合定义。

这种筛选不是简单地在边界图上叠加所有箭头，而是把“边界是否存在”和“该方向场是否由边界局部邻域支持”拆成可审计的 grid-level 指标。每个被选中的 arrow 会保留其 component、grid 坐标、回归方向、R²、梯度强度、箭头长度和头部张角等 LSGI 字段，并合并对应 grid 的 `partition_n_Bdy`、`partition_frac_Bdy`、`local_n_Bdy` 和 `local_frac_Bdy`。模块同时输出 selected grids、selected arrows 和 Bdy spot-arrow assignment 表；后者把每个边界 spot 与其来源 grid 和 component arrow 连接起来，并可附加反卷积细胞组成列。可视化上，模块在 Mal/Bdy/nMal 背景和 H&E 背景上突出显示被筛选的边界相关箭头，未入选箭头可作为灰色背景保留，从而直观看到哪些细胞组分、转录程序或基因的空间变化真正聚焦于边界区域。

### 4. `06_tumor_arrow_guided_boundary_profile.R`：以肿瘤上皮方向场为局部轴的跨边界 profile

`06_tumor_arrow_guided_boundary_profile.R` 是在 Cottrazm-LSGI 耦合结果上的新增下游分析。该模块默认读取 `05` 中 `cell_component` 的边界相关 arrows，并选取 `Cancer.Epithelial` 或 `Tumor` 等肿瘤上皮 component 作为局部方向轴。其核心假设是：肿瘤上皮比例的 LSGI arrow 近似指向局部恶性上皮成分增加的方向，因此可作为每个边界 patch 的“肿瘤侧方向”。与全局 signed-distance 方法相比，该模块允许不同边界片段拥有各自的局部轴，从而适应真实切片中弯曲、分叉、破碎或多灶性的肿瘤边界。

对每支肿瘤上皮边界 arrow，脚本先从 `grid_local_spot_membership.csv` 取出该 grid 的局部 spot，并将 arrow 的 `vx.u/vy.u` 或 `vx/vy` 归一化为单位方向向量。随后把每个局部 spot 投影到 arrow 方向和其垂直方向上，寻找落在轴附近的 `Bdy` spot 作为 0 点；若多个候选存在，则优先选择横向距离更小、纵向距离更近且 local rank 更高的 spot。以该 0 点为中心，脚本在 arrow 方向上建立 `profile_x_px` 坐标，正向代表肿瘤上皮梯度增加方向，负向代表相反方向；在垂直方向上按 `tube_half_width_mode` 定义窄带，只保留局部轴附近的 spot。若同一局部邻域内存在多个分离的 Bdy 簇，脚本会根据 `bdy_cluster_gap_mode` 截断 profile 范围，避免把相邻但不属于同一边界段的 spot 混入同一个 patch。

每个 patch 会被分为 `cross_boundary`、`one_sided`、`no_zero`、`no_direction` 或 `no_local_spots`。主分析默认只纳入 `cross_boundary` patch，并要求窄带内达到最小 spot 数且包含至少一个指定正常侧标签（默认 `nMal`）的 spot；其他 patch 不会被删除，而是保留 exclusion reason，作为 QC 和敏感性分析依据。随后，脚本从 `cell_component`、`marker_module`、`single_gene` 和可选 `pathway` embedding 中抽取用户指定特征，默认包括 redox 相关 module 和基因，并强制加入 `cell_component::Cancer.Epithelial` 作为阳性对照。所有 feature 值按以 0 为中心的距离 bin 汇总为 patch-level bin mean，再汇总为全样本的 `primary_cross_boundary` profile，输出均值、中位数、标准误和 95% CI。

该模块的统计输出包括 `tumor_arrow_profile_patches.csv`、`tumor_arrow_profile_spots.csv`、`tumor_arrow_profile_rectangles.csv`、`selected_features.csv`、`feature_profile_by_patch_bin.csv`、`feature_profile_summary.csv` 和 `feature_boundary_tests.csv`。其中 `feature_boundary_tests.csv` 对主分析 patch 执行 boundary-near 与 tumor-side、boundary-near 与 stroma/opposite-side 的成对差值检验，并估计 tumor-side 和 stroma-side 的 profile slope。对应图形包括每个 feature 的距离曲线、跨 feature 的距离 profile heatmap、patch QC 柱状图、空间 patch 图和窄带矩形覆盖图。因而，该模块不仅回答“某个 feature 是否在 Bdy 高”，还进一步回答其峰值是否位于边界、是否沿肿瘤上皮方向升高，以及是否存在向非恶性侧扩散的局部模式。

### 5. 模块组合的生物学解释框架

这三个关键模块形成了从区域定义到方向解释再到局部定量的证据链。`02_boundary_definition.R` 给每个 spot 建立可复现的肿瘤边界标签，决定了分析对象和空间参照；`05_boundary_related_lsgi_arrows.R` 将 LSGI 的连续方向场限制到 Cottrazm 定义的 Bdy 生态位，识别哪些细胞组分、基因程序或通路在边界处具有方向性；`06_tumor_arrow_guided_boundary_profile.R` 则进一步把肿瘤上皮方向场转化为 patch-level 的局部坐标轴，对 redox、EMT、CAF/ECM、免疫效应、髓系支持、血管和单基因 marker 进行跨边界剖面统计。

这种组合的优势在于避免了三类常见误读。第一，仅靠空间聚类或 CNV 标签容易把肿瘤边界视为静态区域，而 `02` 的迭代邻接扩展强调边界是恶性核心与周围组织之间的连续过渡。第二，仅画全局 LSGI arrows 容易把远离边界的组织梯度误认为肿瘤交界效应，而 `05` 通过 `partition` 和 `local` 两套 grid-aware 策略限定了边界相关性。第三，仅比较 Mal/Bdy/nMal 的平均表达会丢失方向信息，而 `06` 把每个边界 patch 按肿瘤上皮梯度校准到统一的局部轴上，使边界富集、肿瘤侧富集和非恶性侧扩散可以被分开量化。基于这些输出，研究者可以构建 immune-active boundary、immune-excluded boundary、fibrotic/immunosuppressive boundary、invasive-residual boundary 或 redox-protected boundary 等可检验的空间生态位假说，并在样本分组、治疗前后或临床结局之间进行比较。

默认主模块如下：

- `01_spatial_preprocess_cluster.R`：读取 Space Ranger 表达矩阵、spot 坐标和 H&E 图像，构建 Seurat 空间对象，输出 QC，并使用 stLearn/SME 生成 `Morph` assay、PCA、UMAP、聚类和 `NormalScore`。
- `02_boundary_definition.R`：运行 inferCNV，生成统一 `03_cnv_calls.tsv`，把 `CNVLabel` 与 `cnv_score` 写回 Seurat 对象，并结合 CNV、形态校正 UMAP 和空间邻接关系定义 `Mal`、`Bdy`、`nMal`。核心结果是 `intermediate/<样本名>/05_TumorST_boundary_defined.rds.gz`。
- `03_spatial_deconvolution.R`：准备单细胞 signature matrix 和 marker list，估计每个 spot 的细胞类型比例，并输出 LSGI 使用的 spot composition 表。核心结果是 `intermediate/<样本名>/07_DeconData.rds.gz` 和 `output/<样本名>/07_spatial_deconvolution/spot_matrix_pre_lsgi.csv`。
- `04_lsgi_gradient.R`：在组织空间上构建统一 grid，并对 `cell_component`、`nmf`、`marker_module`、`pathway`、`single_gene` 等 embedding/component 来源运行 LSGI 局部线性梯度分析。输出目录保持为 `output/<样本名>/11_lsgi_gradient/`。
- `05_boundary_related_lsgi_arrows.R`：将全局 LSGI arrows 与边界 spot 对齐，按 partition/local/consensus/union 策略筛选边界相关梯度箭头，并输出 arrow、grid 和 Bdy spot-arrow assignment 三层结果。输出目录保持为 `output/<样本名>/12_boundary_related_lsgi_arrows/`。
- `06_tumor_arrow_guided_boundary_profile.R`：使用 `05` 中筛选出的肿瘤上皮 LSGI arrow 作为局部方向轴，沿边界附近窄带统计 redox 或自定义 embedding 的跨边界 profile，并默认加入 `cell_component::Cancer.Epithelial` 作为阳性对照。输出目录为 `output/<样本名>/13_tumor_arrow_guided_boundary_profile/`。


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

`input/CRC1/` 是样本输入目录。`intermediate/CRC1/` 是脚本之间传递的 RDS、TSV 等中间文件目录。`output/CRC1/` 是图、表、inferCNV 输出和可直接用于下游分析的结构化结果目录。

默认完整运行会保留 `intermediate/<样本名>/`，因为其中的 `05_TumorST_boundary_defined.rds.gz` 和 `07_DeconData.rds.gz` 是主线稳定契约的一部分。

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

这些文件会被 `01_spatial_preprocess_cluster.R` 读取并转换为 Seurat 对象。样本名前缀来自 `sample_name`，默认是 `CRC1`。

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

其中 `single_cell_seurat.rds` 应为 Seurat 对象；`define_types.txt` 用于指定或映射细胞类型。`03_spatial_deconvolution.R` 会调用旧子步骤生成：

```text
intermediate/CRC1/06_sig_exp.rds.gz
intermediate/CRC1/06_clustermarkers_list.rds.gz
```

## 运行方式

Linux 上处理样本时，默认方式是通过 `run_all_entry/runall_<sample>.sh` 以 `nohup` 后台提交。这样即使关闭终端，任务也会继续执行，日志统一写到 `run_all_entry/nohup_logs/`。

当前推荐入口：

```bash
bash run_all_entry/runall_CRC1.sh
bash run_all_entry/runall_BRCA1.sh
bash run_all_entry/runall_BRCA2.sh
```

这些 wrapper 默认后台运行，并默认保留 intermediate；如果需要前台调试，再显式加 `--foreground`。

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

也可以分模块运行，例如只跑到反卷积：

```bash
bash run_all_entry/run_all.sh --sample-name BRCA2 --from-step 1 --to-step 3
```

续跑 LSGI 及其下游模块：

```bash
bash run_all_entry/run_all.sh --sample-name BRCA2 --from-step 4 --to-step 6 --resume
```

Linux 下如果不使用 wrapper，而是直接调用通用入口，也建议显式放到后台：

```bash
nohup bash --noprofile --norc run_all_entry/run_all.sh --sample-name BRCA2 > run_all_entry/nohup_logs/BRCA2_manual.nohup.log 2>&1 < /dev/null &
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

CNV 与边界定义会生成：

```text
intermediate/<样本名>/03_cnv_calls.tsv
intermediate/<样本名>/04_TumorST_cnv_scored.rds.gz
intermediate/<样本名>/05_TumorST_boundary_defined.rds.gz
```

反卷积与 LSGI 主线会生成：

```text
intermediate/<样本名>/07_DeconData.rds.gz
output/<样本名>/07_spatial_deconvolution/spot_matrix_pre_lsgi.csv
output/<样本名>/11_lsgi_gradient/component_method_summary.csv
output/<样本名>/11_lsgi_gradient/<component_method>/arrow_tables/arrows_by_grid.csv
output/<样本名>/12_boundary_related_lsgi_arrows/strategy_summary.csv
output/<样本名>/12_boundary_related_lsgi_arrows/<component_method>/arrow_tables/<strategy_id>_arrows.csv
output/<样本名>/12_boundary_related_lsgi_arrows/<component_method>/spot_tables/<strategy_id>_bdy_spot_arrow_assignments.csv
output/<样本名>/13_tumor_arrow_guided_boundary_profile/tumor_arrow_profile_patches.csv
output/<样本名>/13_tumor_arrow_guided_boundary_profile/tumor_arrow_profile_rectangles.csv
output/<样本名>/13_tumor_arrow_guided_boundary_profile/selected_features.csv
output/<样本名>/13_tumor_arrow_guided_boundary_profile/missing_positive_control_features.csv
output/<样本名>/13_tumor_arrow_guided_boundary_profile/feature_profile_summary.csv
output/<样本名>/13_tumor_arrow_guided_boundary_profile/feature_boundary_tests.csv
```

`13_tumor_arrow_guided_boundary_profile` 的主分析只使用通过 QC 的 `cross_boundary` patches：默认要求窄带内至少有 3 个 profile spots，且至少纳入 1 个 `Location == "nMal"` spot。`one_sided` 表示已找到 0 点 Bdy spot，但该 arrow 窄带内保留的 profile spots 只落在 0 点一侧，因此只进入 secondary summary；其他被筛掉的 patch 会保留在 patch/spot/rectangle 表中，并通过 `profile_exclusion_reason` 记录原因。

## 常见检查点

- 输入目录名必须和运行样本名一致，例如 `CRC1` 对应 `input/CRC1/`。
- 新样本请新建 `input/<新样本名>/`，不要覆盖已有样本。
- `output/<样本名>/` 和 `intermediate/<样本名>/` 可以删除后重跑。
- `input/<样本名>/` 是原始输入目录，不应被流程脚本写入或清空。
- 若只需要重新做边界相关箭头筛选，可以在已有 `11_lsgi_gradient/` 结果基础上运行第 5 模块。
- 若只需要重新做 tumor-arrow-guided boundary profile，可以在已有 `11_lsgi_gradient/` 和 `12_boundary_related_lsgi_arrows/` 结果基础上运行第 6 模块。
