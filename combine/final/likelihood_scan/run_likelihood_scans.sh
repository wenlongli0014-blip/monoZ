#!/usr/bin/env bash
set -euo pipefail

BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/combine/final"
OUT="${BASE}/likelihood_scan"
COMBINE_DIR="/afs/cern.ch/user/l/liwe/Combine_v9_2_1_root626"

COUNTING_WS="${BASE}/counting_sr_asimov_cr_data_sim_v921_root626_asimov.root"
SHAPE_WS="${BASE}/shape_sr_asimov_cr_data_sim_v921_root626_asimov.root"

mkdir -p "${OUT}/logs" "${OUT}/root" "${OUT}/plots" "${OUT}/tables"

cd "${COMBINE_DIR}"
set +u
source env_lcg.sh
set -u

cd "${OUT}/root"

COMMON_ARGS=(
  -M MultiDimFit
  -m 125
  --dataset asimovData_1
  --redefineSignalPOIs r
  -P r
  --floatOtherPOIs 1
  --cminDefaultMinimizerStrategy 0
  --robustFit 1
  --setParameterRanges r=-1,2
)

run_one() {
  local label="$1"
  local ws="$2"
  local points="$3"

  combine "${COMMON_ARGS[@]}" "${ws}" \
    -n ".Expected_${label}_singles" \
    --algo singles \
    > "${OUT}/logs/${label}_singles.log" 2>&1

  combine "${COMMON_ARGS[@]}" "${ws}" \
    -n ".Expected_${label}_grid" \
    --algo grid \
    --points "${points}" \
    > "${OUT}/logs/${label}_grid.log" 2>&1
}

run_one "counting" "${COUNTING_WS}" 181
run_one "shape" "${SHAPE_WS}" 181

python3 "${OUT}/plot_likelihood_scan.py"
