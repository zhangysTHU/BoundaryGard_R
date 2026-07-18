#!/usr/bin/env bash
set -euo pipefail

# BRCA1 的当前推荐入口。
# 目的：
# - 固定 sample name 为 BRCA1；
# - 默认走当前 scripts_format_R 新版 03/04 inferCNV 流水线；
# - 支持前台直跑，也支持带 nohup/setsid 的后台提交；
# - 把运行日志和启动脚本集中写到 run_all_entry/nohup_logs。

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_RUNNER="${SCRIPT_DIR}/run_all.sh"
SAMPLE_NAME="BRCA1"
MINIFORGE_ROOT="/lulabdata3/huangkeyun/zhangys/tools/miniforge3"
CONDA_ENV_NAME="BoundaryGrad"
ENV_PREFIX="${MINIFORGE_ROOT}/envs/${CONDA_ENV_NAME}"
TEMP_DIR="${TEMP_DIR:-/tmp/BoundaryGrad_${SAMPLE_NAME}}"
RUN_LOG_DIR="${SCRIPT_DIR}/nohup_logs"
mkdir -p "${RUN_LOG_DIR}"

FROM_STEP="${FROM_STEP:-1}"
TO_STEP="${TO_STEP:-11}"
BACKGROUND="${BACKGROUND:-1}"
RESUME="${RESUME:-0}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
LOW_LOAD="${LOW_LOAD:-1}"

usage() {
  cat <<'EOF'
Usage: runall_BRCA1.sh [options]

Options:
  --from-step N
  --to-step N
  --temp-dir PATH
  --resume
  --keep-intermediate
  --foreground
  --background
  --normal-load
  --help

Examples:
  bash runall_BRCA1.sh --foreground
  bash runall_BRCA1.sh --from-step 3 --to-step 5 --resume
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from-step) FROM_STEP="$2"; shift 2 ;;
    --to-step) TO_STEP="$2"; shift 2 ;;
    --temp-dir) TEMP_DIR="$2"; shift 2 ;;
    --resume) RESUME=1; shift ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
    --foreground) BACKGROUND=0; shift ;;
    --background) BACKGROUND=1; shift ;;
    --normal-load) LOW_LOAD=0; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

cmd=(
  "${ROOT_RUNNER}"
  --sample-name "${SAMPLE_NAME}"
  --from-step "${FROM_STEP}"
  --to-step "${TO_STEP}"
  --temp-dir "${TEMP_DIR}"
)
if (( RESUME == 1 )); then
  cmd+=(--resume)
fi
if (( KEEP_INTERMEDIATE == 1 )); then
  cmd+=(--keep-intermediate)
fi
if (( LOW_LOAD == 0 )); then
  cmd+=(--normal-load)
fi

# 固定到本地 BoundaryGrad 环境，避免依赖用户当前 shell 的 activate 状态。
export COTTRAZM_CONDA_ROOT="${MINIFORGE_ROOT}"
export COTTRAZM_CONDA_ENV="${CONDA_ENV_NAME}"
export COTTRAZM_PYTHON="${ENV_PREFIX}/bin/python"
export COTTRAZM_RSCRIPT="${ENV_PREFIX}/bin/Rscript"
export LOW_LOAD

if (( BACKGROUND == 1 )); then
  timestamp="$(date '+%Y%m%d_%H%M%S')"
  nohup_log="${RUN_LOG_DIR}/BRCA1_${timestamp}.nohup.log"
  launcher_script="${RUN_LOG_DIR}/BRCA1_${timestamp}.launcher.sh"
  runtime_pid_file="${RUN_LOG_DIR}/BRCA1_${timestamp}.runtime.pid"

  cat > "${launcher_script}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
echo "\$\$" > "${runtime_pid_file}"
export COTTRAZM_CONDA_ROOT="${MINIFORGE_ROOT}"
export COTTRAZM_CONDA_ENV="${CONDA_ENV_NAME}"
export COTTRAZM_PYTHON="${ENV_PREFIX}/bin/python"
export COTTRAZM_RSCRIPT="${ENV_PREFIX}/bin/Rscript"
export LOW_LOAD="${LOW_LOAD}"
exec "${ROOT_RUNNER}" --sample-name "${SAMPLE_NAME}" --from-step "${FROM_STEP}" --to-step "${TO_STEP}" --temp-dir "${TEMP_DIR}"$( (( RESUME == 1 )) && printf ' --resume' )$( (( KEEP_INTERMEDIATE == 1 )) && printf ' --keep-intermediate' )$( (( LOW_LOAD == 0 )) && printf ' --normal-load' )
EOF
  chmod +x "${launcher_script}"
  : > "${runtime_pid_file}"

  nohup setsid bash "${launcher_script}" > "${nohup_log}" 2>&1 < /dev/null &
  wrapper_pid="$!"
  echo "${wrapper_pid}" > "${RUN_LOG_DIR}/BRCA1_latest.pid"
  printf '%s\n' "${nohup_log}" > "${RUN_LOG_DIR}/BRCA1_latest.logpath"
  printf '%s\n' "${launcher_script}" > "${RUN_LOG_DIR}/BRCA1_latest.launcher"
  printf '%s\n' "${runtime_pid_file}" > "${RUN_LOG_DIR}/BRCA1_latest.runtime_pid_path"

  echo "Started BRCA1 workflow in background"
  echo "Wrapper PID: ${wrapper_pid}"
  echo "Log: ${nohup_log}"
  echo "Launcher: ${launcher_script}"
else
  exec "${cmd[@]}"
fi
