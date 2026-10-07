#!/usr/bin/env bash
set -euo pipefail

# 通用 R 流水线入口。
# - 默认运行整合后的 6 个主模块。
# - 负责样本目录检查、环境变量注入、日志落盘和按模块执行。

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BOUNDARY_SELECTION_FILE="${ROOT_DIR}/config/boundary_sample_selection.tsv"

FROM_STEP="${FROM_STEP:-1}"
TO_STEP="${TO_STEP:-6}"
SAMPLE_NAME="${SAMPLE_NAME:-CRC1}"
MINIFORGE_ROOT="/lulabdata3/huangkeyun/zhangys/tools/miniforge3"
CONDA_ENV_NAME="BoundaryGrad"
ENV_PREFIX="${MINIFORGE_ROOT}/envs/${CONDA_ENV_NAME}"
R_SCRIPT_BIN="${ENV_PREFIX}/bin/Rscript"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
TEMP_DIR="${TEMP_DIR:-/tmp/BoundaryGrad_${SAMPLE_NAME}}"
LOW_LOAD="${LOW_LOAD:-1}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-1}"
RESUME="${RESUME:-0}"

usage() {
  cat <<'EOF'
Usage: run_all.sh [options]

Options:
  --from-step N
  --to-step N
  --sample-name NAME
  --temp-dir PATH
  --resume
  --keep-intermediate  Keep intermediate/<sample>/ after a full run; this is the default.
  --normal-load
  --help

Default module steps:
  1  01_spatial_preprocess_cluster.R
  2  02_boundary_definition.R
  3  03_spatial_deconvolution.R
  4  04_lsgi_gradient.R
  5  05_boundary_related_lsgi_arrows.R
  6  06_tumor_arrow_guided_boundary_profile.R

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
if (( FROM_STEP < 1 || TO_STEP > 6 || FROM_STEP > TO_STEP )); then
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

BOUNDARY_SELECTION_FOUND=0
BOUNDARY_SELECTED_RUN_ID=""

load_boundary_selection() {
  local -a matching_rows=()
  local expected_header=$'sample_name\tboundary_run_id\tboundary_malignant_cnv_labels\tboundary_mal_cluster_fraction\tboundary_umap_mal_ratio\tboundary_expand_mal_radius\tboundary_max_rounds'
  local actual_header=""
  local row=""
  local table_sample=""
  local table_run_id=""
  local table_labels=""
  local table_fraction=""
  local table_umap_ratio=""
  local table_radius=""
  local table_rounds=""

  if [[ ! -f "${BOUNDARY_SELECTION_FILE}" ]]; then
    echo "Boundary sample selection table not found: ${BOUNDARY_SELECTION_FILE}" >&2
    exit 1
  fi

  IFS= read -r actual_header < "${BOUNDARY_SELECTION_FILE}"
  actual_header="${actual_header%$'\r'}"
  if [[ "${actual_header}" != "${expected_header}" ]]; then
    echo "Unexpected header in boundary sample selection table: ${BOUNDARY_SELECTION_FILE}" >&2
    echo "Expected: ${expected_header}" >&2
    echo "Actual:   ${actual_header}" >&2
    exit 1
  fi

  mapfile -t matching_rows < <(
    awk -F '\t' -v sample="${SAMPLE_NAME}" '
      NR == 1 { next }
      $1 == sample { sub(/\r$/, ""); print }
    ' "${BOUNDARY_SELECTION_FILE}"
  )

  if (( ${#matching_rows[@]} > 1 )); then
    echo "Duplicate boundary selections for sample ${SAMPLE_NAME} in ${BOUNDARY_SELECTION_FILE}" >&2
    exit 1
  fi

  if (( ${#matching_rows[@]} == 0 )); then
    unset COTTRAZM_BOUNDARY_RUN_ID
    unset COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS
    unset COTTRAZM_BOUNDARY_MAL_CLUSTER_FRACTION
    unset COTTRAZM_BOUNDARY_UMAP_MAL_RATIO
    unset COTTRAZM_BOUNDARY_EXPAND_MAL_RADIUS
    unset COTTRAZM_BOUNDARY_MAX_ROUNDS
    echo "Boundary selection: ${SAMPLE_NAME} is not listed; using 00_config.R defaults."
    return
  fi

  row="${matching_rows[0]}"
  IFS=$'\t' read -r table_sample table_run_id table_labels table_fraction table_umap_ratio table_radius table_rounds <<< "${row}"
  table_rounds="${table_rounds%$'\r'}"

  if [[ -z "${table_sample}" || -z "${table_run_id}" || -z "${table_labels}" || \
        -z "${table_fraction}" || -z "${table_umap_ratio}" || -z "${table_radius}" || \
        -z "${table_rounds}" ]]; then
    echo "Incomplete boundary selection row for sample ${SAMPLE_NAME}: ${row}" >&2
    exit 1
  fi
  if [[ ! "${table_run_id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Unsafe boundary_run_id for sample ${SAMPLE_NAME}: ${table_run_id}" >&2
    exit 1
  fi

  BOUNDARY_SELECTION_FOUND=1
  BOUNDARY_SELECTED_RUN_ID="${table_run_id}"
  if [[ "${table_run_id}" == "legacy" ]]; then
    export COTTRAZM_BOUNDARY_RUN_ID=""
  else
    export COTTRAZM_BOUNDARY_RUN_ID="${table_run_id}"
  fi
  export COTTRAZM_BOUNDARY_MALIGNANT_CNV_LABELS="${table_labels}"
  export COTTRAZM_BOUNDARY_MAL_CLUSTER_FRACTION="${table_fraction}"
  export COTTRAZM_BOUNDARY_UMAP_MAL_RATIO="${table_umap_ratio}"
  export COTTRAZM_BOUNDARY_EXPAND_MAL_RADIUS="${table_radius}"
  export COTTRAZM_BOUNDARY_MAX_ROUNDS="${table_rounds}"

  echo "Boundary selection: sample=${SAMPLE_NAME} run=${table_run_id} labels=${table_labels} mal_fraction=${table_fraction} umap_ratio=${table_umap_ratio} radius=${table_radius} rounds=${table_rounds}"
}

activate_selected_boundary() {
  local selected_intermediate_dir=""
  local selected_output_dir=""
  local name=""
  local src=""

  if (( BOUNDARY_SELECTION_FOUND == 0 )) || [[ "${BOUNDARY_SELECTED_RUN_ID}" == "legacy" ]]; then
    return
  fi

  selected_intermediate_dir="${INTERMEDIATE_DIR}/${BOUNDARY_SELECTED_RUN_ID}"
  selected_output_dir="${OUTPUT_DIR}/05_boundary/${BOUNDARY_SELECTED_RUN_ID}"

  for name in 05_TumorST_boundary_subset.rds.gz 05_TumorST_boundary_defined.rds.gz params.tsv location_counts.tsv; do
    src="${selected_intermediate_dir}/${name}"
    if [[ ! -f "${src}" ]]; then
      echo "Selected boundary intermediate is missing: ${src}" >&2
      exit 1
    fi
    cp -f -- "${src}" "${INTERMEDIATE_DIR}/${name}"
  done

  for name in "${SAMPLE_NAME}_BoundaryDefine.pdf" params.tsv location_counts.tsv; do
    src="${selected_output_dir}/${name}"
    if [[ ! -f "${src}" ]]; then
      echo "Selected boundary output is missing: ${src}" >&2
      exit 1
    fi
    cp -f -- "${src}" "${OUTPUT_DIR}/05_boundary/${name}"
  done

  # Historical *_out_*.pdf diagnostics stay inside their versioned run
  # directories.  They are not activated as canonical downstream outputs.
  find "${OUTPUT_DIR}/05_boundary" -maxdepth 1 -type f -name "${SAMPLE_NAME}_out_*.pdf" -delete

  echo "Activated boundary selection ${SAMPLE_NAME}/${BOUNDARY_SELECTED_RUN_ID} for downstream steps."
}

load_boundary_selection

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

if (( FROM_STEP > 2 )); then
  activate_selected_boundary
fi

declare -a STEP_SCRIPTS=(
  "01_spatial_preprocess_cluster.R"
  "02_boundary_definition.R"
  "03_spatial_deconvolution.R"
  "04_lsgi_gradient.R"
  "05_boundary_related_lsgi_arrows.R"
  "06_tumor_arrow_guided_boundary_profile.R"
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
  if (( step_num == 2 )); then
    activate_selected_boundary
  fi
done

if (( KEEP_INTERMEDIATE == 0 && TO_STEP == ${#STEP_SCRIPTS[@]} )); then
  rm -rf "${INTERMEDIATE_DIR}"
fi

echo "[$(date '+%F %T')] R workflow completed for ${SAMPLE_NAME}."
