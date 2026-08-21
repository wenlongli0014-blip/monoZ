#!/usr/bin/env bash
set -euo pipefail

BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/combine/final"
OUT="${BASE}/cls_limit_ranking"
COMBINE_DIR="/afs/cern.ch/user/l/liwe/Combine_v9_2_1_root626"

COUNTING_WS="${BASE}/counting_sr_asimov_cr_data_sim_v921_root626_asimov.root"
SHAPE_WS="${BASE}/shape_sr_asimov_cr_data_sim_v921_root626_asimov.root"

mkdir -p "${OUT}/logs" "${OUT}/root" "${OUT}/plots" "${OUT}/tables" "${OUT}/json"

cd "${COMBINE_DIR}"
set +u
source env_lcg.sh
set -u

cd "${OUT}/root"
rm -f higgsCombine.Expected_*.AsymptoticLimits.mH125.root
rm -f higgsCombine_initialFit_.*_impact.MultiDimFit.mH125.root
rm -f higgsCombine_paramFit_.*_impact_*.MultiDimFit.mH125.root

RANGES="r=-5,5:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5"
COMMON_FIT_ARGS=(
  -m 125
  --dataset asimovData_1
  --redefineSignalPOIs r
  --cminDefaultMinimizerStrategy 0
)

run_limit() {
  local label="$1"
  local ws="$2"

  combine -M AsymptoticLimits "${ws}" \
    "${COMMON_FIT_ARGS[@]}" \
    -n ".Expected_${label}" \
    --rMin 0 --rMax 5 \
    > "${OUT}/logs/${label}_asymptotic_limits.log" 2>&1
}

run_ranking() {
  local label="$1"
  local ws="$2"

  combine -M MultiDimFit "${ws}" \
    "${COMMON_FIT_ARGS[@]}" \
    -n "_initialFit_.${label}_impact" \
    --algo singles \
    --robustFit 1 \
    --setParameterRanges "${RANGES}" \
    > "${OUT}/logs/${label}_impacts_initial.log" 2>&1

  for param in k_Zjet k_WZ k_emu; do
    combine -M MultiDimFit "${ws}" \
      "${COMMON_FIT_ARGS[@]}" \
      -n "_paramFit_.${label}_impact_${param}" \
      --algo impact \
      -P "${param}" \
      --floatOtherPOIs 1 \
      --saveInactivePOI 1 \
      --robustFit 1 \
      --setParameterRanges "${RANGES}" \
      > "${OUT}/logs/${label}_impact_${param}.log" 2>&1
  done
}

run_limit "counting" "${COUNTING_WS}"
run_limit "shape" "${SHAPE_WS}"

run_ranking "counting" "${COUNTING_WS}"
run_ranking "shape" "${SHAPE_WS}"

python3 "${OUT}/plot_cls_limit_ranking.py"
