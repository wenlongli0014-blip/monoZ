#!/bin/bash

set -euo pipefail

REPO="/afs/cern.ch/user/l/liwe/hzz2l2nu"
STUDY_BASE="/eos/user/l/liwe/hadronic_tau_test"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"

if [[ $# -eq 0 ]]; then
  variants=(tau_medium tau_vloose tau_vvvloose tau_none)
else
  variants=("$@")
fi

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

for variant in "${variants[@]}"; do
  input_base="${STUDY_BASE}/${variant}/combine_input"
  output_base="${STUDY_BASE}/${variant}/combine"
  shapes="${output_base}/shapes_ptmiss_2017_2018_statonly.root"
  crfit_workspace="${output_base}/workspaces/ptmiss_crfit.root"
  final_workspace="${output_base}/workspaces/ptmiss_final.root"
  crfit_result="${output_base}/fitDiagnostics.${variant}_ptmiss_crfit.root"
  closure_result="${output_base}/fitDiagnostics.${variant}_ptmiss_expected_closure.root"
  mkdir -p \
    "${output_base}/cards" \
    "${output_base}/fits" \
    "${output_base}/logs" \
    "${output_base}/workspaces"
  cd "${output_base}"

  echo "[${variant}] Build blinded fixed-binning ptmiss shapes"
  python3 "${REPO}/combine/make_shapes_ptmiss_2017_2018_blinded.py" \
    --input-base "${input_base}" \
    --mode shape \
    --out "${shapes}" \
    2>&1 | tee "${output_base}/logs/01_make_shapes.log"

  echo "[${variant}] Build three-control-region fit card"
  python3 "${REPO}/combine/fit_check/scripts/make_fit_check_datacard.py" \
    --shapes "${shapes}" \
    --tag ptmiss_crfit \
    --channels DYCR CR3L EMUCR \
    --k-factors k_Zjet k_WZ k_emu \
    --dummy-signal \
    --out-dir "${output_base}/cards" \
    2>&1 | tee "${output_base}/logs/02_make_crfit_card.log"

  echo "[${variant}] Build control-region workspace"
  text2workspace.py \
    "${output_base}/cards/ptmiss_crfit.txt" \
    -o "${crfit_workspace}" \
    2>&1 | tee "${output_base}/logs/03_text2workspace_crfit.log"

  echo "[${variant}] Fit actual DYCR, CR3L, and EMUCR data"
  combine -M FitDiagnostics \
    -d "${crfit_workspace}" \
    -m 125 \
    --redefineSignalPOIs k_Zjet \
    --setParameters r=0 \
    --freezeParameters r \
    --saveShapes \
    --saveWithUncertainties \
    --saveNormalizations \
    --cminDefaultMinimizerStrategy 0 \
    -n ".${variant}_ptmiss_crfit" \
    2>&1 | tee "${output_base}/logs/04_crfit.log"

  python3 "${REPO}/bin/check_combine_fit_result.py" \
    --fit "${crfit_result}" \
    --params k_Zjet k_WZ k_emu \
    --out "${output_base}/fits/cr_fit_summary.json" \
    2>&1 | tee "${output_base}/logs/05_check_crfit.log"

  k_zjet="$(
    python3 "${REPO}/combine/fit_check/scripts/read_fit_param.py" \
      --fit "${crfit_result}" --param k_Zjet
  )"
  k_wz="$(
    python3 "${REPO}/combine/fit_check/scripts/read_fit_param.py" \
      --fit "${crfit_result}" --param k_WZ
  )"
  k_emu="$(
    python3 "${REPO}/combine/fit_check/scripts/read_fit_param.py" \
      --fit "${crfit_result}" --param k_emu
  )"
  fitted_parameters="k_Zjet=${k_zjet},k_WZ=${k_wz},k_emu=${k_emu}"
  printf '%s\n' "${fitted_parameters}" \
    | tee "${output_base}/fits/fitted_parameters.txt"

  echo "[${variant}] Build full SR+CR stat-only card"
  python3 "${REPO}/combine/fit_check/scripts/make_fit_check_datacard.py" \
    --shapes "${shapes}" \
    --tag ptmiss_final \
    --channels SR DYCR CR3L EMUCR \
    --k-factors k_Zjet k_WZ k_emu \
    --out-dir "${output_base}/cards" \
    2>&1 | tee "${output_base}/logs/06_make_final_card.log"

  if grep -Eq '^[[:space:]]*[^#].*autoMCStats' \
    "${output_base}/cards/ptmiss_final.txt"; then
    echo "ERROR: autoMCStats is enabled in the final card" >&2
    exit 1
  fi

  echo "[${variant}] Build full SR+CR workspace"
  text2workspace.py \
    "${output_base}/cards/ptmiss_final.txt" \
    -o "${final_workspace}" \
    2>&1 | tee "${output_base}/logs/07_text2workspace_final.log"

  echo "[${variant}] Check fitted-k background-only Asimov closure"
  combine -M FitDiagnostics \
    -d "${final_workspace}" \
    -m 125 \
    -t -1 \
    --expectSignal 0 \
    --setParameters "${fitted_parameters}" \
    --redefineSignalPOIs r \
    --rMin -5 \
    --rMax 5 \
    --saveShapes \
    --saveWithUncertainties \
    --saveNormalizations \
    --cminDefaultMinimizerStrategy 0 \
    -n ".${variant}_ptmiss_expected_closure" \
    2>&1 | tee "${output_base}/logs/08_expected_closure.log"

  python3 "${REPO}/bin/check_combine_fit_result.py" \
    --fit "${closure_result}" \
    --params r k_Zjet k_WZ k_emu \
    --out "${output_base}/fits/expected_closure_summary.json" \
    2>&1 | tee "${output_base}/logs/09_check_expected_closure.log"

  echo "[${variant}] Run fitted-k background-only Asimov expected limit"
  combine -M AsymptoticLimits \
    -d "${final_workspace}" \
    -m 125 \
    -t -1 \
    --expectSignal 0 \
    --setParameters "${fitted_parameters}" \
    --redefineSignalPOIs r \
    --rMin 0 \
    --rMax 5 \
    --cminDefaultMinimizerStrategy 0 \
    -n ".${variant}_ptmiss_crconditioned_expected" \
    2>&1 | tee "${output_base}/logs/10_limit_expected.log"
done
