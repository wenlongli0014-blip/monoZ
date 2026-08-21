#!/usr/bin/env bash
set -euo pipefail

DNN_BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"
OUT="${DNN_BASE}/three_panel_results"

mkdir -p "${OUT}/logs" "${OUT}/root" "${OUT}/plots" "${OUT}/tables"

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

run_one() {
  local label="$1"
  local workspace="$2"
  local ranges="$3"
  shift 3
  local params=("$@")

  echo "[$(date)] ${label}: MultiDimFit singles"
  combine -M MultiDimFit "${workspace}" \
    -m 125 \
    --dataset asimovData_1 \
    --redefineSignalPOIs r \
    -P r \
    --floatOtherPOIs 1 \
    --cminDefaultMinimizerStrategy 0 \
    --robustFit 1 \
    --setParameterRanges "r=-1,2" \
    -n ".${label}_singles" \
    --algo singles \
    > "${OUT}/logs/${label}_singles.log" 2>&1
  mv "higgsCombine.${label}_singles.MultiDimFit.mH125.root" "${OUT}/root/" 2>/dev/null || true

  echo "[$(date)] ${label}: MultiDimFit grid"
  combine -M MultiDimFit "${workspace}" \
    -m 125 \
    --dataset asimovData_1 \
    --redefineSignalPOIs r \
    -P r \
    --floatOtherPOIs 1 \
    --cminDefaultMinimizerStrategy 0 \
    --robustFit 1 \
    --setParameterRanges "r=-1,2" \
    -n ".${label}_grid" \
    --algo grid \
    --points 181 \
    > "${OUT}/logs/${label}_grid.log" 2>&1
  mv "higgsCombine.${label}_grid.MultiDimFit.mH125.root" "${OUT}/root/" 2>/dev/null || true

  echo "[$(date)] ${label}: ranking initial"
  combine -M MultiDimFit "${workspace}" \
    -m 125 \
    --dataset asimovData_1 \
    --redefineSignalPOIs r \
    --cminDefaultMinimizerStrategy 0 \
    -n "_initialFit_.${label}_impact" \
    --algo singles \
    --robustFit 1 \
    --setParameterRanges "${ranges}" \
    > "${OUT}/logs/${label}_impacts_initial.log" 2>&1
  mv "higgsCombine_initialFit_.${label}_impact.MultiDimFit.mH125.root" "${OUT}/root/" 2>/dev/null || true

  for param in "${params[@]}"; do
    echo "[$(date)] ${label}: impact ${param}"
    combine -M MultiDimFit "${workspace}" \
      -m 125 \
      --dataset asimovData_1 \
      --redefineSignalPOIs r \
      --cminDefaultMinimizerStrategy 0 \
      -n "_paramFit_.${label}_impact_${param}" \
      --algo impact \
      -P "${param}" \
      --floatOtherPOIs 1 \
      --saveInactivePOI 1 \
      --robustFit 1 \
      --setParameterRanges "${ranges}" \
      > "${OUT}/logs/${label}_impact_${param}.log" 2>&1
    mv "higgsCombine_paramFit_.${label}_impact_${param}.MultiDimFit.mH125.root" "${OUT}/root/" 2>/dev/null || true
  done
}

cd "${OUT}/root"
rm -f ./*.root
rm -f "${OUT}/logs"/*.log

run_one \
  "dnn2017" \
  "${DNN_BASE}/2017/workspaces/workspace_dnn_score_2017_asimov.root" \
  "r=-5,5:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5" \
  k_Zjet k_WZ k_emu

run_one \
  "dnn2018" \
  "${DNN_BASE}/2018/workspaces/workspace_dnn_score_2018_asimov.root" \
  "r=-5,5:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5" \
  k_Zjet k_WZ k_emu

run_one \
  "dnn2017_2018" \
  "${DNN_BASE}/combined_2017_2018/workspaces/workspace_dnn_score_2017_2018_asimov.root" \
  "r=-5,5:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5" \
  k_Zjet k_WZ k_emu

python3 "${DNN_BASE}/plot_three_panel_dnn_results.py"

echo "Done."
echo "Plots: ${OUT}/plots"
