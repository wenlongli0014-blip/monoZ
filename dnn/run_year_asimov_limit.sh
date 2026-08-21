#!/usr/bin/env bash
set -euo pipefail

YEAR="${1:?Usage: run_year_asimov_limit.sh YEAR}"

DNN_BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
YEAR_DIR="${DNN_BASE}/${YEAR}"
CONFIG="${YEAR_DIR}/dnn_config_${YEAR}.json"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"
LCG_VIEW="/cvmfs/sft.cern.ch/lcg/views/LCG_109_cuda/x86_64-el9-gcc13-opt/setup.sh"
TAG="dnn_score_${YEAR}"

INPUT_BASE="/eos/user/l/liwe/monoZ_combine_17_18_dnn/${YEAR}_dnn"
SHAPES_DATA="${YEAR_DIR}/shapes_${TAG}_data.root"
SHAPES_ASIMOV="${YEAR_DIR}/shapes_${TAG}_sr_asimov.root"
CARD="${YEAR_DIR}/cards/combined_${TAG}.txt"
WORKSPACE_INPUT="${YEAR_DIR}/workspaces/workspace_${TAG}.root"
WORKSPACE_ASIMOV="${YEAR_DIR}/workspaces/workspace_${TAG}_asimov.root"
LIMIT_DIR="${YEAR_DIR}/limits"

mkdir -p \
  "${YEAR_DIR}/logs" \
  "${YEAR_DIR}/cards" \
  "${YEAR_DIR}/workspaces" \
  "${LIMIT_DIR}"

cd "${DNN_BASE}"

echo "[$(date)] Build DNN-score shapes for ${YEAR}"
set +u
source "${LCG_VIEW}"
set -u
python3 "${DNN_BASE}/make_dnn_shapes.py" \
  --config "${CONFIG}" \
  --input-base "${INPUT_BASE}" \
  --signal-file signal.root \
  --out "${SHAPES_DATA}" \
  > "${YEAR_DIR}/logs/03_make_shapes.log" 2>&1

echo "[$(date)] Replace SR data_obs with background-only Asimov; keep CR data_obs real"
cp -f "${SHAPES_DATA}" "${SHAPES_ASIMOV}"
python3 /afs/cern.ch/user/l/liwe/hzz2l2nu/combine/make_pseudo_asimov_data.py \
  --shapes "${SHAPES_ASIMOV}" \
  --asimov-channels SR \
  > "${YEAR_DIR}/logs/03b_sr_asimov.log" 2>&1

echo "[$(date)] Build datacard/workspace for ${YEAR}"
cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

cd "${DNN_BASE}"
python3 "${DNN_BASE}/make_dnn_datacard.py" \
  --shapes "${SHAPES_ASIMOV}" \
  --tag "${TAG}" \
  --outdir "${YEAR_DIR}/cards" \
  > "${YEAR_DIR}/logs/04_make_datacard.log" 2>&1

text2workspace.py "${CARD}" \
  -o "${WORKSPACE_INPUT}" \
  > "${YEAR_DIR}/logs/05_text2workspace.log" 2>&1

echo "[$(date)] Generate fitted Asimov dataset as asimovData_1 for ${YEAR}"
cd "${YEAR_DIR}/workspaces"
rm -f higgsCombine.MakeAsimov_${TAG}.FitDiagnostics.mH125*.root fitDiagnostics.MakeAsimov_${TAG}.root "${WORKSPACE_ASIMOV}"

combine -M FitDiagnostics "${WORKSPACE_INPUT}" \
  -m 125 \
  -n ".MakeAsimov_${TAG}" \
  --redefineSignalPOIs r \
  -t -1 \
  --expectSignal 0 \
  --toysFrequentist \
  --saveToys \
  --saveWorkspace \
  --cminDefaultMinimizerStrategy 0 \
  --rMin -5 --rMax 10 \
  --setParameterRanges "r=-5,10:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5" \
  > "${YEAR_DIR}/logs/06_make_asimov_fitdiagnostics.log" 2>&1

ASIMOV_TOY_FILE="$(ls -1 higgsCombine.MakeAsimov_${TAG}.FitDiagnostics.mH125.*.root | head -n 1)"
python3 "${DNN_BASE}/make_dnn_asimov_workspace.py" \
  --input-workspace "${WORKSPACE_INPUT}" \
  --toy-file "${YEAR_DIR}/workspaces/${ASIMOV_TOY_FILE}" \
  --output-workspace "${WORKSPACE_ASIMOV}" \
  > "${YEAR_DIR}/logs/07_import_asimov_workspace.log" 2>&1

echo "[$(date)] Run expected CLs limit for ${YEAR}"
cd "${LIMIT_DIR}"
rm -f "higgsCombine.Expected_${TAG}.AsymptoticLimits.mH125.root"
combine -M AsymptoticLimits "${WORKSPACE_ASIMOV}" \
  -m 125 \
  --dataset asimovData_1 \
  --redefineSignalPOIs r \
  --cminDefaultMinimizerStrategy 0 \
  -n ".Expected_${TAG}" \
  --rMin 0 --rMax 5 \
  > "${YEAR_DIR}/logs/08_expected_limit.log" 2>&1

echo "[$(date)] Done ${YEAR}"
echo "Shapes data: ${SHAPES_DATA}"
echo "Shapes SR-Asimov/CR-data: ${SHAPES_ASIMOV}"
echo "Workspace: ${WORKSPACE_ASIMOV}"
echo "Limit log: ${YEAR_DIR}/logs/08_expected_limit.log"
