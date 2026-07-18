#!/usr/bin/env bash
set -euo pipefail

# 通用 R 流水线入口。
# - 默认按当前 scripts_format_R 的新 03/04 inferCNV 逻辑运行。
# - 负责样本目录检查、环境变量注入、日志落盘和按步执行。

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

FROM_STEP="${FROM_STEP:-1}"
TO_STEP="${TO_STEP:-11}"
SAMPLE_NAME="${SAMPLE_NAME:-CRC1}"
MINIFORGE_ROOT="/lulabdata3/huangkeyun/zhangys/tools/miniforge3"
CONDA_ENV_NAME="BoundaryGrad"
ENV_PREFIX="${MINIFORGE_ROOT}/envs/${CONDA_ENV_NAME}"
R_SCRIPT_BIN="${ENV_PREFIX}/bin/Rscript"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
TEMP_DIR="${TEMP_DIR:-/tmp/BoundaryGrad_${SAMPLE_NAME}}"
LOW_LOAD="${LOW_LOAD:-1}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
RESUME="${RESUME:-0}"
STEP11_SCRIPT="${COTTRAZM_STEP11_SCRIPT:-11_lsgi_gradient.R}"

usage() {
  cat <<'EOF'
Usage: run_all.sh [options]

Options:
  --from-step N
  --to-step N
  --sample-name NAME
  --temp-dir PATH
  --resume
  --keep-intermediate
  --normal-load
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from-step) FROM_STEP="$2"; shift 2 ;;
    --to-step) TO_STEP="$2"; shift 2 ;;
    --sample-name) SAMPLE_NAME="$2"; shift 2 ;;
    --temp-dir) TEMP_DIR="$2"; shift 2 ;;
    --resume) RESUME=1; shift ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
    --normal-load) LOW_LOAD=0; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ ! "$FROM_STEP" =~ ^[0-9]+$ ]] || [[ ! "$TO_STEP" =~ ^[0-9]+$ ]]; then
  echo "From/to step must be integers" >&2
  exit 1
fi
if (( FROM_STEP < 1 || TO_STEP > 11 || FROM_STEP > TO_STEP )); then
  echo "Invalid step range: ${FROM_STEP}-${TO_STEP}" >&2
  exit 1
fi

INPUT_DIR="${ROOT_DIR}/input/${SAMPLE_NAME}"
INTERMEDIATE_DIR="${ROOT_DIR}/intermediate/${SAMPLE_NAME}"
OUTPUT_DIR="${ROOT_DIR}/output/${SAMPLE_NAME}"
LOG_DIR="${OUTPUT_DIR}/run_logs"

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Sample input directory not found: ${INPUT_DIR}" >&2
  exit 1
fi

mkdir -p "${TEMP_DIR}" "${LOG_DIR}"

export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH
export LANG=C
export LC_ALL=C
export TEMP="${TEMP_DIR}"
export TMP="${TEMP_DIR}"
export TMPDIR="${TEMP_DIR}"
export COTTRAZM_SAMPLE_NAME="${SAMPLE_NAME}"
export COTTRAZM_FAST_CNV="${COTTRAZM_FAST_CNV:-}"
export COTTRAZM_INFERCNV_THREADS="${COTTRAZM_INFERCNV_THREADS:-2}"
export COTTRAZM_INFERCNV_PARTITION="${COTTRAZM_INFERCNV_PARTITION:-qnorm}"
export COTTRAZM_INFERCNV_ANALYSIS_MODE="${COTTRAZM_INFERCNV_ANALYSIS_MODE:-subclusters}"
export COTTRAZM_INFERCNV_USE_CHECKPOINT="${COTTRAZM_INFERCNV_USE_CHECKPOINT:-1}"
export COTTRAZM_INFERCNV_POSTPROCESS_ONLY="${COTTRAZM_INFERCNV_POSTPROCESS_ONLY:-}"
export COTTRAZM_CONDA_ROOT="${MINIFORGE_ROOT}"
export COTTRAZM_CONDA_ENV="${CONDA_ENV_NAME}"
export COTTRAZM_PYTHON="${PYTHON_BIN}"
export COTTRAZM_RSCRIPT="${R_SCRIPT_BIN}"
export COTTRAZM_USE_RENV="${COTTRAZM_USE_RENV:-}"

if [[ ! -x "${R_SCRIPT_BIN}" ]]; then
  echo "Rscript not executable: ${R_SCRIPT_BIN}" >&2
  exit 1
fi
if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Python not executable: ${PYTHON_BIN}" >&2
  exit 1
fi

if (( RESUME == 0 )); then
  if (( FROM_STEP != 1 )); then
    echo "A clean run must start at step 1. Use --resume for continuation." >&2
    exit 1
  fi
  rm -rf "${INTERMEDIATE_DIR}" "${OUTPUT_DIR}"
  mkdir -p "${LOG_DIR}"
fi

declare -a STEP_SCRIPTS=(
  "01_preprocess_st.R"
  "02_morphology_adjusted_cluster.R"
  "03_run_infercnv.R"
  "04_score_cnv.R"
  "05_define_boundary.R"
  "06_prepare_single_cell_reference.R"
  "07_spatial_deconvolution.R"
  "08_spatial_reconstruction.R"
  "09_diff_and_enrichment.R"
  "10_plot_results.R"
  "${STEP11_SCRIPT}"
)

run_one_step() {
  local step_num="$1"
  local script_name="$2"
  local log_file="${LOG_DIR}/step_$(printf '%02d' "${step_num}")_${script_name%.R}.log"
  local stdout_file="${LOG_DIR}/step_$(printf '%02d' "${step_num}")_${script_name%.R}.stdout.log"
  local stderr_file="${LOG_DIR}/step_$(printf '%02d' "${step_num}")_${script_name%.R}.stderr.log"
  echo "[$(date '+%F %T')] Running ${script_name}" | tee -a "${log_file}"
  if (( LOW_LOAD == 1 )); then
    nice -n 10 "${R_SCRIPT_BIN}" --vanilla "${ROOT_DIR}/${script_name}" \
      > >(tee -a "${stdout_file}" | tee -a "${log_file}") \
      2> >(tee -a "${stderr_file}" | tee -a "${log_file}" >&2)
  else
    "${R_SCRIPT_BIN}" --vanilla "${ROOT_DIR}/${script_name}" \
      > >(tee -a "${stdout_file}" | tee -a "${log_file}") \
      2> >(tee -a "${stderr_file}" | tee -a "${log_file}" >&2)
  fi
}

for i in "${!STEP_SCRIPTS[@]}"; do
  step_num=$((i + 1))
  if (( step_num < FROM_STEP || step_num > TO_STEP )); then
    continue
  fi
  run_one_step "${step_num}" "${STEP_SCRIPTS[$i]}"
done

if (( KEEP_INTERMEDIATE == 0 && TO_STEP == 11 )); then
  rm -rf "${INTERMEDIATE_DIR}"
fi

echo "[$(date '+%F %T')] R workflow completed for ${SAMPLE_NAME}."
