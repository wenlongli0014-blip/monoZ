#!/usr/bin/env bash
set -euo pipefail

BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"
TAG="dnn_score_real_data"
CARD="${BASE}/cards/combined_${TAG}.txt"
WORKSPACE="${BASE}/workspaces/workspace_${TAG}.root"
SHAPES="${BASE}/shapes_dnn_score_data.root"

OUT="${BASE}/real_data_fit"
mkdir -p "${OUT}/logs" "${OUT}/root" "${OUT}/plots" "${OUT}/tables" "${OUT}/json" "${BASE}/cards" "${BASE}/workspaces"

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

cd "${BASE}"
python3 "${BASE}/make_dnn_datacard.py" \
  --shapes "${SHAPES}" \
  --tag "${TAG}" \
  --outdir "${BASE}/cards" \
  > "${OUT}/logs/make_datacard.log" 2>&1

text2workspace.py "${CARD}" \
  -o "${WORKSPACE}" \
  > "${OUT}/logs/text2workspace.log" 2>&1

RANGES="r=-1,2:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5"
COMMON_ARGS=(
  -m 125
  --redefineSignalPOIs r
  --cminDefaultMinimizerStrategy 0
)

cd "${OUT}/root"
rm -f higgsCombine.RealData_dnn_score.AsymptoticLimits.mH125.root
rm -f higgsCombine_initialFit_.real_data_dnn_score_impact.MultiDimFit.mH125.root
rm -f higgsCombine_paramFit_.real_data_dnn_score_impact_*.MultiDimFit.mH125.root
rm -f higgsCombine.RealData_dnn_score_singles.MultiDimFit.mH125.root
rm -f higgsCombine.RealData_dnn_score_grid.MultiDimFit.mH125.root

combine -M AsymptoticLimits "${WORKSPACE}" \
  "${COMMON_ARGS[@]}" \
  -n ".RealData_dnn_score" \
  --rMin 0 --rMax 5 \
  > "${OUT}/logs/dnn_score_asymptotic_limits.log" 2>&1

combine -M MultiDimFit "${WORKSPACE}" \
  "${COMMON_ARGS[@]}" \
  -n "_initialFit_.real_data_dnn_score_impact" \
  --algo singles \
  --robustFit 1 \
  --setParameterRanges "${RANGES}" \
  > "${OUT}/logs/dnn_score_impacts_initial.log" 2>&1

for param in k_Zjet k_WZ k_emu; do
  combine -M MultiDimFit "${WORKSPACE}" \
    "${COMMON_ARGS[@]}" \
    -n "_paramFit_.real_data_dnn_score_impact_${param}" \
    --algo impact \
    -P "${param}" \
    --floatOtherPOIs 1 \
    --saveInactivePOI 1 \
    --robustFit 1 \
    --setParameterRanges "${RANGES}" \
    > "${OUT}/logs/dnn_score_impact_${param}.log" 2>&1
done

SCAN_ARGS=(
  -M MultiDimFit
  -m 125
  --redefineSignalPOIs r
  -P r
  --floatOtherPOIs 1
  --cminDefaultMinimizerStrategy 0
  --robustFit 1
  --setParameterRanges r=-1,2
)

combine "${SCAN_ARGS[@]}" "${WORKSPACE}" \
  -n ".RealData_dnn_score_singles" \
  --algo singles \
  > "${OUT}/logs/dnn_score_singles.log" 2>&1

combine "${SCAN_ARGS[@]}" "${WORKSPACE}" \
  -n ".RealData_dnn_score_grid" \
  --algo grid \
  --points 181 \
  > "${OUT}/logs/dnn_score_grid.log" 2>&1

python3 "${BASE}/plot_dnn_real_data_results.py"

echo "Done."
echo "Real-data plots: ${OUT}/plots"
