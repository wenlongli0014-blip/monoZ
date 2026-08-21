#!/usr/bin/env bash
set -euo pipefail

DNN_BASE="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"
OUT="${DNN_BASE}/combined_2017_2018"
TAG="dnn_score_2017_2018"

mkdir -p "${OUT}/cards" "${OUT}/logs" "${OUT}/workspaces" "${OUT}/limits"

CARD_2017_IN="${DNN_BASE}/2017/cards/combined_dnn_score_2017.txt"
CARD_2018_IN="${DNN_BASE}/2018/cards/combined_dnn_score_2018.txt"
CARD_2017="${OUT}/cards/combined_dnn_score_2017_sharedparams.txt"
CARD_2018="${OUT}/cards/combined_dnn_score_2018_sharedparams.txt"
CARD_COMBINED="${OUT}/cards/combined_${TAG}.txt"
WORKSPACE_INPUT="${OUT}/workspaces/workspace_${TAG}.root"
WORKSPACE_ASIMOV="${OUT}/workspaces/workspace_${TAG}_asimov.root"

if [[ ! -f "${CARD_2017_IN}" || ! -f "${CARD_2018_IN}" ]]; then
  echo "Missing per-year cards. Run dnn/run_year_asimov_limit.sh 2017 and 2018 first." >&2
  exit 1
fi

python3 - "${CARD_2017_IN}" "${CARD_2017}" <<'PY'
from pathlib import Path
import sys

src, dst = sys.argv[1:]
Path(dst).write_text(Path(src).read_text())
PY

python3 - "${CARD_2018_IN}" "${CARD_2018}" <<'PY'
from pathlib import Path
import sys

src, dst = sys.argv[1:]
Path(dst).write_text(Path(src).read_text())
PY

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

echo "[$(date)] Combine 2017 and 2018 cards"
combineCards.py y2017="${CARD_2017}" y2018="${CARD_2018}" > "${CARD_COMBINED}"
python3 - "${CARD_COMBINED}" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
lines = path.read_text().splitlines()
for idx, line in enumerate(lines):
    if line.startswith("kmax "):
        lines[idx] = "kmax * number of nuisance parameters"
        break
path.write_text("\n".join(lines) + "\n")
PY

echo "[$(date)] text2workspace"
text2workspace.py "${CARD_COMBINED}" \
  -o "${WORKSPACE_INPUT}" \
  > "${OUT}/logs/01_text2workspace.log" 2>&1

PARAM_RANGES="r=-5,10:k_Zjet=0,5:k_WZ=0,5:k_emu=0,5"

echo "[$(date)] Generate combined fitted Asimov dataset"
cd "${OUT}/workspaces"
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
  --setParameterRanges "${PARAM_RANGES}" \
  > "${OUT}/logs/02_make_asimov_fitdiagnostics.log" 2>&1

ASIMOV_TOY_FILE="$(ls -1 higgsCombine.MakeAsimov_${TAG}.FitDiagnostics.mH125.*.root | head -n 1)"
python3 "${DNN_BASE}/make_dnn_asimov_workspace.py" \
  --input-workspace "${WORKSPACE_INPUT}" \
  --toy-file "${OUT}/workspaces/${ASIMOV_TOY_FILE}" \
  --output-workspace "${WORKSPACE_ASIMOV}" \
  > "${OUT}/logs/03_import_asimov_workspace.log" 2>&1

echo "[$(date)] Run combined expected CLs limit"
cd "${OUT}/limits"
rm -f "higgsCombine.Expected_${TAG}.AsymptoticLimits.mH125.root"
combine -M AsymptoticLimits "${WORKSPACE_ASIMOV}" \
  -m 125 \
  --dataset asimovData_1 \
  --redefineSignalPOIs r \
  --cminDefaultMinimizerStrategy 0 \
  -n ".Expected_${TAG}" \
  --rMin 0 --rMax 5 \
  > "${OUT}/logs/04_expected_limit.log" 2>&1

echo "[$(date)] Done combined 2017+2018"
echo "Card: ${CARD_COMBINED}"
echo "Workspace: ${WORKSPACE_ASIMOV}"
echo "Limit log: ${OUT}/logs/04_expected_limit.log"
