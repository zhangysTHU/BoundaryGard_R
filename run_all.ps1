param(
  [int]$FromStep = 1,
  [int]$ToStep = 11,
  [string]$SampleName = "CRC1",
  [string]$Rscript = "Rscript"
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$previousSampleName = [Environment]::GetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", "Process")

$rscriptCommand = Get-Command $Rscript -ErrorAction SilentlyContinue
if (-not $rscriptCommand) {
  throw "Rscript not found: $Rscript"
}
$Rscript = $rscriptCommand.Source

[Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $SampleName, "Process")

$defaults = @{
  COTTRAZM_FAST_CNV = ""
  COTTRAZM_INFERCNV_THREADS = "2"
  COTTRAZM_INFERCNV_PARTITION = "qnorm"
  COTTRAZM_INFERCNV_ANALYSIS_MODE = "subclusters"
  COTTRAZM_INFERCNV_USE_CHECKPOINT = "1"
  COTTRAZM_INFERCNV_POSTPROCESS_ONLY = ""
}

foreach ($key in $defaults.Keys) {
  if ([string]::IsNullOrEmpty([Environment]::GetEnvironmentVariable($key, "Process"))) {
    [Environment]::SetEnvironmentVariable($key, $defaults[$key], "Process")
  }
}

$steps = @(
  @{ Step = 1; Script = "01_preprocess_st.R" },
  @{ Step = 2; Script = "02_morphology_adjusted_cluster.R" },
  @{ Step = 3; Script = "03_run_infercnv.R" },
  @{ Step = 4; Script = "04_score_cnv.R" },
  @{ Step = 5; Script = "05_define_boundary.R" },
  @{ Step = 6; Script = "06_prepare_single_cell_reference.R" },
  @{ Step = 7; Script = "07_spatial_deconvolution.R" },
  @{ Step = 8; Script = "08_spatial_reconstruction.R" },
  @{ Step = 9; Script = "09_diff_and_enrichment.R" },
  @{ Step = 10; Script = "10_plot_results.R" },
  @{ Step = 11; Script = "11_lsgi_gradient.R" }
)

Push-Location $root
try {
  foreach ($item in $steps) {
    if ($item.Step -lt $FromStep -or $item.Step -gt $ToStep) {
      continue
    }
    Write-Host ("[{0}] Running {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $item.Script)
    & $Rscript $item.Script
    if ($LASTEXITCODE -ne 0) {
      throw "$($item.Script) failed with exit code $LASTEXITCODE"
    }
  }
  Write-Host ("[{0}] R workflow completed." -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"))
}
finally {
  [Environment]::SetEnvironmentVariable("COTTRAZM_SAMPLE_NAME", $previousSampleName, "Process")
  Pop-Location
}
