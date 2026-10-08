#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MANIFEST="${MANIFEST:-${SCRIPT_DIR}/new_samples_2026-10-07.tsv}"
STATUS_DIR="${PIPELINE_DIR}/output/run_logs"
STATUS_FILE="${STATUS_DIR}/new_samples_batch_status.tsv"
RUNNER="${SCRIPT_DIR}/run_all.sh"

mkdir -p "${STATUS_DIR}"
printf 'sample_name\tstatus\tstarted_at\tfinished_at\tdetail\n' > "${STATUS_FILE}"

export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH
export LANG=C
export LC_ALL=C
export COTTRAZM_INFERCNV_THREADS="${COTTRAZM_INFERCNV_THREADS:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"

record_status() {
  local sample="$1"
  local state="$2"
  local started="$3"
  local finished="$4"
  local detail="$5"
  detail="${detail//$'\t'/ }"
  detail="${detail//$'\n'/ }"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${sample}" "${state}" "${started}" "${finished}" "${detail}" >> "${STATUS_FILE}"
}

failure_count=0
completed_count=0
target_count=0
while IFS=$'\t' read -r sample_name; do
  sample_name="${sample_name%$'\r'}"
  [[ -z "${sample_name}" || "${sample_name}" == "sample_name" ]] && continue
  target_count=$((target_count + 1))
  input_dir="${PIPELINE_DIR}/input/${sample_name}"
  output_dir="${PIPELINE_DIR}/output/${sample_name}"
  completion_marker="${output_dir}/run_logs/pipeline_complete.tsv"
  started_at="$(date '+%F %T')"

  if [[ -s "${completion_marker}" ]]; then
    record_status "${sample_name}" "complete" "${started_at}" "${started_at}" "existing completion marker"
    completed_count=$((completed_count + 1))
    continue
  fi
  if [[ ! -d "${input_dir}" ]]; then
    record_status "${sample_name}" "failed" "${started_at}" "$(date '+%F %T')" "input directory missing"
    failure_count=$((failure_count + 1))
    continue
  fi

  printf '[%s] START %s\n' "${started_at}" "${sample_name}"
  if "${RUNNER}" --sample-name "${sample_name}" --keep-intermediate; then
    finished_at="$(date '+%F %T')"
    printf 'sample_name\tcompleted_at\tsteps\n%s\t%s\t1-6\n' \
      "${sample_name}" "${finished_at}" > "${completion_marker}"
    record_status "${sample_name}" "complete" "${started_at}" "${finished_at}" "steps 1-6 completed"
    completed_count=$((completed_count + 1))
    printf '[%s] COMPLETE %s\n' "${finished_at}" "${sample_name}"
  else
    finished_at="$(date '+%F %T')"
    record_status "${sample_name}" "failed" "${started_at}" "${finished_at}" \
      "see ${output_dir}/run_logs"
    failure_count=$((failure_count + 1))
    printf '[%s] FAILED %s; continuing with next sample\n' "${finished_at}" "${sample_name}" >&2
  fi
done < "${MANIFEST}"

printf 'target_samples\t%s\ncompleted\t%s\nfailed\t%s\nfinished_at\t%s\n' \
  "${target_count}" "${completed_count}" "${failure_count}" "$(date '+%F %T')" \
  > "${STATUS_DIR}/new_samples_batch_summary.tsv"

exit 0
