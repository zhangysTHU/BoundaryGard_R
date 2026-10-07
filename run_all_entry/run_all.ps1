param(
  [int]$FromStep = 1,
  [int]$ToStep = 6,
  [string]$SampleName = "CRC1",
  [string]$Rscript = ""
)

$ErrorActionPreference = "Stop"
# 通用 R 流水线入口。
# 默认运行整合后的 6 个主模块，并沿用旧关键中间文件和输出目录。
$entryRoot = $PSScriptRoot
$root = Split-Path -Parent $entryRoot
$boundarySelectionFile = Join-Path $root "config/boundary_sample_selection.tsv"
$previousSampleName = [Environment]::GetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", "Process")
$boundaryEnvironmentNames = @(
  "COTTRAZM_BOUNDARY_RUN_ID",
  "COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS",
  "COTTRAZM_BOUNDARY_MAL_CLUSTER_FRACTION",
  "COTTRAZM_BOUNDARY_UMAP_MAL_RATIO",
  "COTTRAZM_BOUNDARY_EXPAND_MAL_RADIUS",
  "COTTRAZM_BOUNDARY_MAX_ROUNDS"
)
$previousBoundaryEnvironment = @{}
foreach ($name in $boundaryEnvironmentNames) {
  $previousBoundaryEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, "Process")
}

if ([string]::IsNullOrWhiteSpace($Rscript)) {
  $Rscript = [Environment]::GetEnvironmentVariable("COTTRAZM_RSCRIPT", "Process")
}
if ([string]::IsNullOrWhiteSpace($Rscript)) {
  $Rscript = "Rscript"
}
$rscriptCommand = Get-Command $Rscript -ErrorAction SilentlyContinue
if (-not $rscriptCommand) {
  throw "Rscript command not found: $Rscript"
}
$Rscript = $rscriptCommand.Source

[Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $SampleName, "Process")

if (-not (Test-Path -LiteralPath $boundarySelectionFile -PathType Leaf)) {
  throw "Boundary sample selection table not found: $boundarySelectionFile"
}
$boundaryRows = @(Import-Csv -LiteralPath $boundarySelectionFile -Delimiter "`t" | Where-Object { $_.sample_name -eq $SampleName })
if ($boundaryRows.Count -gt 1) {
  throw "Duplicate boundary selections for sample $SampleName in $boundarySelectionFile"
}

$boundarySelectionFound = $boundaryRows.Count -eq 1
$boundarySelectedRunId = ""
if ($boundarySelectionFound) {
  $boundaryRow = $boundaryRows[0]
  $requiredColumns = @(
    "sample_name", "boundary_run_id", "boundary_malignant_cnv_labels",
    "boundary_mal_cluster_fraction", "boundary_umap_mal_ratio",
    "boundary_expand_mal_radius", "boundary_max_rounds"
  )
  foreach ($column in $requiredColumns) {
    if ([string]::IsNullOrWhiteSpace($boundaryRow.$column)) {
      throw "Incomplete boundary selection for sample $SampleName: missing $column"
    }
  }
  if ($boundaryRow.boundary_run_id -notmatch '^[A-Za-z0-9._-]+$') {
    throw "Unsafe boundary_run_id for sample ${SampleName}: $($boundaryRow.boundary_run_id)"
  }

  $boundarySelectedRunId = $boundaryRow.boundary_run_id
  $runIdForR = if ($boundarySelectedRunId -eq "legacy") { $null } else { $boundarySelectedRunId }
  [Environment]::SetEnvironmentVariable("COTTRAZM_BOUNDARY_RUN_ID", $runIdForR, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS", $boundaryRow.boundary_malignant_cnv_labels, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_BOUNDARY_MAL_CLUSTER_FRACTION", $boundaryRow.boundary_mal_cluster_fraction, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_BOUNDARY_UMAP_MAL_RATIO", $boundaryRow.boundary_umap_mal_ratio, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_BOUNDARY_EXPAND_MAL_RADIUS", $boundaryRow.boundary_expand_mal_radius, "Process")
  [Environment]::SetEnvironmentVariable("COTTRAZM_BOUNDARY_MAX_ROUNDS", $boundaryRow.boundary_max_rounds, "Process")
  Write-Host ("Boundary selection: sample={0} run={1} labels={2} mal_fraction={3} umap_ratio={4} radius={5} rounds={6}" -f $SampleName, $boundarySelectedRunId, $boundaryRow.boundary_malignant_cnv_labels, $boundaryRow.boundary_mal_cluster_fraction, $boundaryRow.boundary_umap_mal_ratio, $boundaryRow.boundary_expand_mal_radius, $boundaryRow.boundary_max_rounds)
}
else {
  foreach ($name in $boundaryEnvironmentNames) {
    [Environment]::SetEnvironmentVariable($name, $null, "Process")
  }
  Write-Host "Boundary selection: $SampleName is not listed; using 00_config.R defaults."
}

function Activate-SelectedBoundary {
  if (-not $boundarySelectionFound -or $boundarySelectedRunId -eq "legacy") {
    return
  }

  $intermediateRoot = Join-Path $root "intermediate/$SampleName"
  $outputRoot = Join-Path $root "output/$SampleName/05_boundary"
  $selectedIntermediate = Join-Path $intermediateRoot $boundarySelectedRunId
  $selectedOutput = Join-Path $outputRoot $boundarySelectedRunId

  foreach ($name in @("05_TumorST_boundary_subset.rds.gz", "05_TumorST_boundary_defined.rds.gz", "params.tsv", "location_counts.tsv")) {
    $source = Join-Path $selectedIntermediate $name
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
      throw "Selected boundary intermediate is missing: $source"
    }
    Copy-Item -LiteralPath $source -Destination (Join-Path $intermediateRoot $name) -Force
  }

  foreach ($name in @("${SampleName}_BoundaryDefine.pdf", "params.tsv", "location_counts.tsv")) {
    $source = Join-Path $selectedOutput $name
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
      throw "Selected boundary output is missing: $source"
    }
    Copy-Item -LiteralPath $source -Destination (Join-Path $outputRoot $name) -Force
  }

  # Historical *_out_*.pdf diagnostics stay inside their versioned run
  # directories and are not activated as canonical downstream outputs.
  Get-ChildItem -LiteralPath $outputRoot -File -Filter "${SampleName}_out_*.pdf" -ErrorAction SilentlyContinue | Remove-Item -Force
  Write-Host "Activated boundary selection $SampleName/$boundarySelectedRunId for downstream steps."
}

$defaults = @{
  COTTRAZM_FAST_CNV = ""
  COTTRAZM_INFERCNV_THREADS = "2"
  COTTRAZM_INFERCNV_PARTITION = "qnorm"
  COTTRAZM_INFERCNV_ANALYSIS_MODE = "subclusters"
  COTTRAZM_INFERCNV_USE_CHECKPOINT = "1"
  COTTRAZM_INFERCNV_POSTPROCESS_ONLY = ""
  COTTRAZM_CONDA_ROOT = "/lulabdata3/huangkeyun/zhangys/tools/miniforge3"
  COTTRAZM_CONDA_ENV = "BoundaryGrad"
  COTTRAZM_PYTHON = "/lulabdata3/huangkeyun/zhangys/tools/miniforge3/envs/BoundaryGrad/bin/python"
  COTTRAZM_RSCRIPT = "/lulabdata3/huangkeyun/zhangys/tools/miniforge3/envs/BoundaryGrad/bin/Rscript"
  COTTRAZM_USE_RENV = ""
}

foreach ($key in $defaults.Keys) {
  if ([string]::IsNullOrEmpty([Environment]::GetEnvironmentVariable($key, "Process"))) {
    [Environment]::SetEnvironmentVariable($key, $defaults[$key], "Process")
  }
}

$steps = @(
  @{ Step = 1; Script = "01_spatial_preprocess_cluster.R" },
  @{ Step = 2; Script = "02_boundary_definition.R" },
  @{ Step = 3; Script = "03_spatial_deconvolution.R" },
  @{ Step = 4; Script = "04_lsgi_gradient.R" },
  @{ Step = 5; Script = "05_boundary_related_lsgi_arrows.R" },
  @{ Step = 6; Script = "06_tumor_arrow_guided_boundary_profile.R" }
)

Push-Location $root
try {
  if ($FromStep -gt 2) {
    Activate-SelectedBoundary
  }
  foreach ($item in $steps) {
    if ($item.Step -lt $FromStep -or $item.Step -gt $ToStep) {
      continue
    }
    Write-Host ("[{0}] Running {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $item.Script)
    & $Rscript --vanilla $item.Script
    if ($LASTEXITCODE -ne 0) {
      throw "$($item.Script) failed with exit code $LASTEXITCODE"
    }
    if ($item.Step -eq 2) {
      Activate-SelectedBoundary
    }
  }
  Write-Host ("[{0}] R workflow completed." -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"))
}
finally {
  [Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $previousSampleName, "Process")
  foreach ($name in $boundaryEnvironmentNames) {
    [Environment]::SetEnvironmentVariable($name, $previousBoundaryEnvironment[$name], "Process")
  }
  Pop-Location
}
