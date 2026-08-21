#!/usr/bin/env bash
set -euo pipefail

BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"
WORKSPACE_INPUT="${BASE}/workspaces/workspace_dnn_score.root"
WORKSPACE="${BASE}/workspaces/workspace_dnn_score_asimov.root"

CLS_OUT="${BASE}/cls_limit_ranking"
SCAN_OUT="${BASE}/likelihood_scan"

mkdir -p "${CLS_OUT}/logs" "${CLS_OUT}/root" "${CLS_OUT}/plots" "${CLS_OUT}/tables" "${CLS_OUT}/json"
mkdir -p "${SCAN_OUT}/logs" "${SCAN_OUT}/root" "${SCAN_OUT}/plots" "${SCAN_OUT}/tables"

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

RANGES="r=-5,5:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5"
COMMON_ARGS=(
  -m 125
  --dataset asimovData_1
  --redefineSignalPOIs r
  --cminDefaultMinimizerStrategy 0
)

cd "${BASE}/workspaces"
rm -f higgsCombine.DnnScoreAsimov.FitDiagnostics.mH125*.root
rm -f fitDiagnostics.DnnScoreAsimov.root
rm -f "${WORKSPACE}"

combine -M FitDiagnostics "${WORKSPACE_INPUT}" \
  -m 125 \
  -n ".DnnScoreAsimov" \
  --redefineSignalPOIs r \
  -t -1 \
  --expectSignal 0 \
  --toysFrequentist \
  --saveToys \
  --saveWorkspace \
  --cminDefaultMinimizerStrategy 0 \
  --rMin -5 --rMax 10 \
  --setParameterRanges "r=-5,10:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5" \
  > "${BASE}/logs/dnn_score_make_asimov_fitdiagnostics.log" 2>&1

ASIMOV_TOY_FILE="$(ls -1 higgsCombine.DnnScoreAsimov.FitDiagnostics.mH125.*.root | head -n 1)"
python3 "${BASE}/make_dnn_asimov_workspace.py" \
  --input-workspace "${WORKSPACE_INPUT}" \
  --toy-file "${BASE}/workspaces/${ASIMOV_TOY_FILE}" \
  --output-workspace "${WORKSPACE}" \
  > "${BASE}/logs/dnn_score_make_asimov_workspace.log" 2>&1

cd "${CLS_OUT}/root"
rm -f higgsCombine.Expected_dnn_score.AsymptoticLimits.mH125.root
rm -f higgsCombine_initialFit_.dnn_score_impact.MultiDimFit.mH125.root
rm -f higgsCombine_paramFit_.dnn_score_impact_*.MultiDimFit.mH125.root

combine -M AsymptoticLimits "${WORKSPACE}" \
  "${COMMON_ARGS[@]}" \
  -n ".Expected_dnn_score" \
  --rMin 0 --rMax 5 \
  > "${CLS_OUT}/logs/dnn_score_asymptotic_limits.log" 2>&1

combine -M MultiDimFit "${WORKSPACE}" \
  "${COMMON_ARGS[@]}" \
  -n "_initialFit_.dnn_score_impact" \
  --algo singles \
  --robustFit 1 \
  --setParameterRanges "${RANGES}" \
  > "${CLS_OUT}/logs/dnn_score_impacts_initial.log" 2>&1

for param in k_Zjet k_WZ k_emu; do
  combine -M MultiDimFit "${WORKSPACE}" \
    "${COMMON_ARGS[@]}" \
    -n "_paramFit_.dnn_score_impact_${param}" \
    --algo impact \
    -P "${param}" \
    --floatOtherPOIs 1 \
    --saveInactivePOI 1 \
    --robustFit 1 \
    --setParameterRanges "${RANGES}" \
    > "${CLS_OUT}/logs/dnn_score_impact_${param}.log" 2>&1
done

python3 "${BASE}/plot_dnn_cls_limit_ranking.py"

cd "${SCAN_OUT}/root"
rm -f higgsCombine.Expected_dnn_score_singles.MultiDimFit.mH125.root
rm -f higgsCombine.Expected_dnn_score_grid.MultiDimFit.mH125.root

SCAN_ARGS=(
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

combine "${SCAN_ARGS[@]}" "${WORKSPACE}" \
  -n ".Expected_dnn_score_singles" \
  --algo singles \
  > "${SCAN_OUT}/logs/dnn_score_singles.log" 2>&1

combine "${SCAN_ARGS[@]}" "${WORKSPACE}" \
  -n ".Expected_dnn_score_grid" \
  --algo grid \
  --points 181 \
  > "${SCAN_OUT}/logs/dnn_score_grid.log" 2>&1

python3 "${BASE}/plot_dnn_likelihood_scan.py"

echo "Done."
echo "CLs/ranking plots: ${CLS_OUT}/plots"
echo "Likelihood scan plots: ${SCAN_OUT}/plots"
