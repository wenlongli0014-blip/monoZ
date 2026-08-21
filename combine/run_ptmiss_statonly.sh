#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/afs/cern.ch/user/l/liwe/hzz2l2nu/combine"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"

mkdir -p "${WORKDIR}/logs"

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

cd "${WORKDIR}"

echo "[1/6] make shapes"
python3 make_shapes_ptmiss_statonly.py | tee logs/01_make_shapes.log

echo "[2/6] make datacards"
python3 make_datacards_ptmiss_statonly.py | tee logs/02_make_cards.log

echo "[3/6] text2workspace"
text2workspace.py cards/combined_statonly.txt -o workspace_ptmiss_statonly.root | tee logs/03_text2workspace.log

echo "[4/6] FitDiagnostics"
combine -M FitDiagnostics \
  -d workspace_ptmiss_statonly.root \
  -m 125 \
  -n .ptmiss_fit \
  --saveShapes --saveWithUncertainties --saveNormalizations \
  | tee logs/04_fitdiagnostics.log

echo "[5/6] Asymptotic observed"
combine -M AsymptoticLimits \
  -d workspace_ptmiss_statonly.root \
  -m 125 \
  -n .ptmiss_obs \
  | tee logs/05_limit_observed.log

echo "[6/6] Asymptotic expected (Asimov)"
combine -M AsymptoticLimits \
  -d workspace_ptmiss_statonly.root \
  -m 125 \
  -t -1 \
  -n .ptmiss_exp \
  | tee logs/06_limit_expected.log

echo "Done. Main outputs:"
echo "  ${WORKDIR}/workspace_ptmiss_statonly.root"
echo "  ${WORKDIR}/higgsCombine.ptmiss_obs.AsymptoticLimits.mH125.root"
echo "  ${WORKDIR}/higgsCombine.ptmiss_exp.AsymptoticLimits.mH125.root"
