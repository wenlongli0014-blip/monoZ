#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <variable> [tag]"
  echo "  variable: ptmiss | mT | ptmiss_significance_corrected"
  echo "  tag     : optional output tag (default: variable)"
  exit 1
fi

VAR="$1"
RAW_TAG="${2:-$VAR}"
TAG="$(echo "$RAW_TAG" | sed -E 's/[^A-Za-z0-9._-]/_/g')"

WORKDIR="/afs/cern.ch/user/l/liwe/hzz2l2nu/combine"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"

LOGDIR="${WORKDIR}/logs/${TAG}"
SHAPES="shapes_${TAG}_statonly.root"
CARD="cards/combined_${TAG}_statonly.txt"
WORKSPACE="workspace_${TAG}_statonly.root"

mkdir -p "${LOGDIR}"

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

cd "${WORKDIR}"

echo "[1/6] make shapes (${VAR})"
python3 make_shapes_var_statonly.py --var "${VAR}" --out "${SHAPES}" | tee "${LOGDIR}/01_make_shapes.log"

echo "[2/6] make datacards (${TAG})"
python3 make_datacards_var_statonly.py --shapes "${SHAPES}" --tag "${TAG}" | tee "${LOGDIR}/02_make_cards.log"

echo "[3/6] text2workspace (${TAG})"
text2workspace.py "${CARD}" -o "${WORKSPACE}" | tee "${LOGDIR}/03_text2workspace.log"

echo "[4/6] FitDiagnostics (${TAG})"
combine -M FitDiagnostics \
  -d "${WORKSPACE}" \
  -m 125 \
  -n ".${TAG}_fit" \
  --saveShapes --saveWithUncertainties --saveNormalizations \
  | tee "${LOGDIR}/04_fitdiagnostics.log"

echo "[5/6] Asymptotic observed (${TAG})"
combine -M AsymptoticLimits \
  -d "${WORKSPACE}" \
  -m 125 \
  -n ".${TAG}_obs" \
  | tee "${LOGDIR}/05_limit_observed.log"

echo "[6/6] Asymptotic expected (Asimov) (${TAG})"
combine -M AsymptoticLimits \
  -d "${WORKSPACE}" \
  -m 125 \
  -t -1 \
  -n ".${TAG}_exp" \
  | tee "${LOGDIR}/06_limit_expected.log"

echo "Done. Main outputs:"
echo "  ${WORKDIR}/${WORKSPACE}"
echo "  ${WORKDIR}/higgsCombine.${TAG}_obs.AsymptoticLimits.mH125.root"
echo "  ${WORKDIR}/higgsCombine.${TAG}_exp.AsymptoticLimits.mH125.root"