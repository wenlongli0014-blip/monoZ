#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
INPUT_BASE="/eos/user/l/liwe/monoZ_combine"
OUTPUT_BASE="/eos/user/l/liwe/monoZ_combine_dnn"
LCG_VIEW="/cvmfs/sft.cern.ch/lcg/views/LCG_109_cuda/x86_64-el9-gcc13-opt/setup.sh"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"
TAG="dnn_score"
CONFIG="${WORKDIR}/dnn_config.json"

mkdir -p "${WORKDIR}/logs" "${WORKDIR}/output" "${WORKDIR}/cards" "${WORKDIR}/workspaces" "${WORKDIR}/limits"

cd "${WORKDIR}"

echo "[1/6] Train weighted DNN"
set +u
source "${LCG_VIEW}"
set -u
python3 train_dnn_score.py \
  --config "${CONFIG}" \
  --input-base "${INPUT_BASE}" \
  --signal-file signal.root \
  --outdir "${WORKDIR}/output" \
  --epochs 50 \
  --patience 10 \
  --max-events-per-process 0 \
  | tee "${WORKDIR}/logs/01_train.log"

echo "[2/6] Add dnn_score branch to copied ROOT files"
python3 add_dnn_score_to_root.py \
  --input-base "${INPUT_BASE}" \
  --output-base "${OUTPUT_BASE}" \
  --model "${WORKDIR}/output/dnn_model.pt" \
  --metadata "${WORKDIR}/output/dnn_metadata.json" \
  --signal-file signal.root \
  --overwrite \
  | tee "${WORKDIR}/logs/02_add_scores.log"

echo "[3/6] Build DNN score shapes"
python3 make_dnn_shapes.py \
  --config "${CONFIG}" \
  --input-base "${OUTPUT_BASE}" \
  --signal-file signal.root \
  --out "${WORKDIR}/shapes_${TAG}_data.root" \
  | tee "${WORKDIR}/logs/03_make_shapes.log"

echo "[3b/6] Replace SR data_obs with background-only Asimov; keep CR data_obs as real data"
cp -f "${WORKDIR}/shapes_${TAG}_data.root" "${WORKDIR}/shapes_${TAG}_sr_asimov.root"
python3 /afs/cern.ch/user/l/liwe/hzz2l2nu/combine/make_pseudo_asimov_data.py \
  --shapes "${WORKDIR}/shapes_${TAG}_sr_asimov.root" \
  --asimov-channels SR \
  | tee "${WORKDIR}/logs/03b_sr_asimov.log"

echo "[4/6] Enter CMSSW/Combine environment and build datacard"
cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

cd "${WORKDIR}"
python3 make_dnn_datacard.py \
  --shapes "${WORKDIR}/shapes_${TAG}_sr_asimov.root" \
  --tag "${TAG}" \
  --outdir "${WORKDIR}/cards" \
  | tee "${WORKDIR}/logs/04_make_datacard.log"

echo "[5/6] text2workspace"
text2workspace.py "${WORKDIR}/cards/combined_${TAG}.txt" \
  -o "${WORKDIR}/workspaces/workspace_${TAG}.root" \
  | tee "${WORKDIR}/logs/05_text2workspace.log"

echo "[6/6] Expected Asimov limit"
combine -M AsymptoticLimits \
  -d "${WORKDIR}/workspaces/workspace_${TAG}.root" \
  -m 125 \
  -n ".${TAG}_exp" \
  | tee "${WORKDIR}/logs/06_limit_expected.log"

mv "higgsCombine.${TAG}_exp.AsymptoticLimits.mH125.root" "${WORKDIR}/limits/" 2>/dev/null || true

echo "Done."
echo "Scored ROOT files: ${OUTPUT_BASE}"
echo "Shapes: ${WORKDIR}/shapes_${TAG}.root"
echo "Card: ${WORKDIR}/cards/combined_${TAG}.txt"
echo "Workspace: ${WORKDIR}/workspaces/workspace_${TAG}.root"
