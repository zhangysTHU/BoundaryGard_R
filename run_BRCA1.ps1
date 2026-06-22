param(
  [int]$FromStep = 1,
  [int]$ToStep = 11,
  [string]$Rscript = "Rscript",
  [string]$TempDir = "",
  [switch]$Resume,
  [switch]$KeepIntermediate
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$sampleName = "BRCA1"

if ($FromStep -lt 1 -or $ToStep -gt 11 -or $FromStep -gt $ToStep) {
  throw "Invalid step range: FromStep=$FromStep, ToStep=$ToStep"
}

$rscriptCommand = Get-Command $Rscript -ErrorAction SilentlyContinue
if (-not $rscriptCommand) {
  throw "Rscript not found: $Rscript"
}
$Rscript = $rscriptCommand.Source

$inputDir = Join-Path $root "input\$sampleName"
$intermediateDir = Join-Path $root "intermediate\$sampleName"
$outputDir = Join-Path $root "output\$sampleName"
$renvLibrary = Join-Path $root "renv\library\R-4.3\x86_64-w64-mingw32"

if (-not (Test-Path -LiteralPath $inputDir -PathType Container)) {
  throw "BRCA1 input directory not found: $inputDir"
}

if (-not (Test-Path -LiteralPath $renvLibrary -PathType Container)) {
  throw "Project renv library not found: $renvLibrary"
}

if ([string]::IsNullOrWhiteSpace($TempDir)) {
  $driveRoot = [System.IO.Path]::GetPathRoot($root)
  $TempDir = Join-Path $driveRoot "cottrazm_tmp\$sampleName"
}
$TempDir = [System.IO.Path]::GetFullPath($TempDir)

if ($TempDir -match "\s" -or $TempDir -match "[^\x00-\x7F]") {
  throw "TempDir must contain only ASCII characters and no spaces: $TempDir"
}

function Remove-SampleDirectory {
  param([Parameter(Mandatory = $true)][string]$Path)

  $fullPath = [System.IO.Path]::GetFullPath($Path)
  $fullRoot = [System.IO.Path]::GetFullPath($root)
  if (-not $fullPath.StartsWith(
      $fullRoot + [System.IO.Path]::DirectorySeparatorChar,
      [System.StringComparison]::OrdinalIgnoreCase
    )) {
    throw "Refusing to delete a path outside scripts_format_R: $fullPath"
  }

  if (Test-Path -LiteralPath $fullPath) {
    Remove-Item -LiteralPath $fullPath -Recurse -Force
  }
}

if (-not $Resume) {
  if ($FromStep -ne 1) {
    throw "A clean run must start at step 1. Use -Resume for a partial continuation."
  }
  Write-Host "Removing previous BRCA1 intermediate files and outputs..."
  Remove-SampleDirectory -Path $intermediateDir
  Remove-SampleDirectory -Path $outputDir
}

New-Item -ItemType Directory -Path $TempDir -Force | Out-Null

$previousEnvironment = @{}
$environment = @{
  TEMP                              = $TempDir
  TMP                               = $TempDir
  TMPDIR                            = $TempDir
  R_LIBS_USER                       = $renvLibrary
  COTTRAZM_SAMPLE_NAME              = $sampleName
  COTTRAZM_FAST_CNV                 = ""
  COTTRAZM_INFERCNV_THREADS         = "2"
  COTTRAZM_INFERCNV_PARTITION       = "qnorm"
  COTTRAZM_INFERCNV_ANALYSIS_MODE   = "subclusters"
  COTTRAZM_INFERCNV_USE_CHECKPOINT  = "1"
  COTTRAZM_INFERCNV_POSTPROCESS_ONLY = ""
}

foreach ($key in $environment.Keys) {
  $previousEnvironment[$key] = [Environment]::GetEnvironmentVariable($key, "Process")
  [Environment]::SetEnvironmentVariable($key, $environment[$key], "Process")
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
    & $Rscript --vanilla $item.Script
    if ($LASTEXITCODE -ne 0) {
      throw "$($item.Script) failed with exit code $LASTEXITCODE"
    }
  }

  Write-Host ("[{0}] BRCA1 R workflow completed." -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"))

  if (-not $KeepIntermediate -and $ToStep -eq 11) {
    Write-Host "Removing BRCA1 intermediate files..."
    Remove-SampleDirectory -Path $intermediateDir
  }
}
finally {
  foreach ($key in $previousEnvironment.Keys) {
    [Environment]::SetEnvironmentVariable($key, $previousEnvironment[$key], "Process")
  }
  Pop-Location
}
