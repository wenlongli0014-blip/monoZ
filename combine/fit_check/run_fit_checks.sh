#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/afs/cern.ch/user/l/liwe/hzz2l2nu/combine/fit_check"
CMSSW_SRC="/afs/cern.ch/user/l/liwe/CMSSW_14_1_0_pre4/src"

cd "${CMSSW_SRC}"
unset PYTHONHOME
unset PYTHONPATH
export SCRAM_ARCH=el9_amd64_gcc12
eval "$(scramv1 runtime -sh)"

cd "${WORKDIR}"
mkdir -p shapes cards workspaces fits logs tables plots

rm -f tables/fit_parameters.csv tables/fit_status.csv tables/yield_summary.csv \
      tables/k_factors_used.csv tables/manifest.tsv

printf "tag\tcategory\tdescription\tshapes\tcard\tworkspace\tfit_root\tlog\n" > tables/manifest.tsv

run_fit() {
  local tag="$1"
  local category="$2"
  local description="$3"
  local shapes="$4"
  local channels="$5"
  local k_factors="$6"
  local pois="$7"
  local allow_r="$8"

  local card="cards/${tag}.txt"
  local workspace="workspaces/${tag}.root"
  local fit_dir="fits/${tag}"
  local log="logs/${tag}_fitdiagnostics.log"
  local dummy_signal_args=()
  mkdir -p "${fit_dir}"

  if [[ " ${channels} " != *" SR "* ]]; then
    dummy_signal_args+=(--dummy-signal)
  fi

  echo "[CARD] ${tag}"
  if [[ -n "${k_factors}" ]]; then
    python3 scripts/make_fit_check_datacard.py \
      --shapes "${shapes}" \
      --tag "${tag}" \
      --channels ${channels} \
      --k-factors ${k_factors} \
      "${dummy_signal_args[@]}" \
      > "logs/${tag}_make_card.log"
  else
    python3 scripts/make_fit_check_datacard.py \
      --shapes "${shapes}" \
      --tag "${tag}" \
      --channels ${channels} \
      "${dummy_signal_args[@]}" \
      > "logs/${tag}_make_card.log"
  fi

  echo "[WORKSPACE] ${tag}"
  text2workspace.py "${card}" -o "${workspace}" > "logs/${tag}_text2workspace.log" 2>&1

  local workspace_abs
  workspace_abs="$(readlink -f "${workspace}")"
  local fit_args=(
    -M FitDiagnostics
    -d "${workspace_abs}"
    -m 125
    -n ".${tag}"
    --saveShapes
    --saveWithUncertainties
    --saveNormalizations
    --cminDefaultMinimizerStrategy 0
  )
  if [[ -n "${pois}" ]]; then
    fit_args+=(--redefineSignalPOIs "${pois}")
  fi
  if [[ " ${channels} " != *" SR "* ]]; then
    fit_args+=(--setParameters r=0 --freezeParameters r)
  fi
  if [[ "${allow_r}" == "yes" ]]; then
    fit_args+=(--rMin -5 --rMax 10)
  fi

  echo "[FIT] ${tag}"
  (
    cd "${fit_dir}"
    combine "${fit_args[@]}"
  ) > "${log}" 2>&1

  local fit_root="${fit_dir}/fitDiagnostics.${tag}.root"
  python3 scripts/collect_fit_check_results.py \
    --tag "${tag}" \
    --fit "${fit_root}" \
    --shapes "${shapes}" \
    --channels ${channels}

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${tag}" "${category}" "${description}" "${shapes}" "${card}" "${workspace}" "${fit_root}" "${log}" \
    >> tables/manifest.tsv
}

echo "[1/7] Build real-data ptmiss shapes"
python3 scripts/make_fit_check_shapes.py \
  --mode shape \
  --out shapes/shape_ptmiss.root \
  > logs/shape_ptmiss_make_shapes.log 2>&1

echo "[2/7] Build SR-Asimov ptmiss shapes"
python3 scripts/make_fit_check_shapes.py \
  --mode shape \
  --out shapes/shape_ptmiss_sr_asimov.root \
  --asimov-channels SR \
  > logs/shape_ptmiss_sr_asimov_make_shapes.log 2>&1

echo "[3/7] Build real-data counting shapes"
python3 scripts/make_fit_check_shapes.py \
  --mode counting \
  --out shapes/counting_onebin.root \
  > logs/counting_onebin_make_shapes.log 2>&1

echo "[4/7] Build SR-Asimov counting shapes"
python3 scripts/make_fit_check_shapes.py \
  --mode counting \
  --out shapes/counting_onebin_sr_asimov.root \
  --asimov-channels SR \
  > logs/counting_onebin_sr_asimov_make_shapes.log 2>&1

echo "[5/7] Shape fits with ptmiss"
run_fit "shape_cronly_dycr" "shape" "CR only data fit: Z+jet CR -> k_Zjet" \
  "shapes/shape_ptmiss.root" "DYCR" "k_Zjet" "k_Zjet" "no"
run_fit "shape_cronly_cr3l" "shape" "CR only data fit: 3l CR -> k_WZ" \
  "shapes/shape_ptmiss.root" "CR3L" "k_WZ" "k_WZ" "no"
run_fit "shape_cronly_emucr" "shape" "CR only data fit: emu CR -> k_emu" \
  "shapes/shape_ptmiss.root" "EMUCR" "k_emu" "k_emu" "no"

K_SHAPE_ZJET="$(python3 scripts/read_fit_param.py --fit fits/shape_cronly_dycr/fitDiagnostics.shape_cronly_dycr.root --param k_Zjet)"
K_SHAPE_WZ="$(python3 scripts/read_fit_param.py --fit fits/shape_cronly_cr3l/fitDiagnostics.shape_cronly_cr3l.root --param k_WZ)"
K_SHAPE_EMU="$(python3 scripts/read_fit_param.py --fit fits/shape_cronly_emucr/fitDiagnostics.shape_cronly_emucr.root --param k_emu)"

python3 scripts/scale_fit_check_templates.py \
  --input shapes/shape_ptmiss.root \
  --output shapes/shape_ptmiss_sr_scaled_by_shape_cr.root \
  --channels SR \
  --k-zjet "${K_SHAPE_ZJET}" \
  --k-wz "${K_SHAPE_WZ}" \
  --k-emu "${K_SHAPE_EMU}" \
  > logs/shape_sr_scale_templates.log 2>&1

{
  echo "source,mode,k_Zjet,k_WZ,k_emu"
  echo "CR-only ptmiss shape,shape,${K_SHAPE_ZJET},${K_SHAPE_WZ},${K_SHAPE_EMU}"
} > tables/k_factors_used.csv

run_fit "shape_sr_asimov_nok" "shape" "SR only Asimov fit, no k-factors" \
  "shapes/shape_ptmiss_sr_asimov.root" "SR" "" "" "yes"
run_fit "shape_sr_data_withk" "shape" "SR only data fit, floating k_Zjet/k_WZ/k_emu" \
  "shapes/shape_ptmiss.root" "SR" "k_Zjet k_WZ k_emu" "" "yes"
run_fit "shape_sr_data_scaledk" "shape" "SR only data fit, templates scaled by CR-only k-factors, no k-factors in card" \
  "shapes/shape_ptmiss_sr_scaled_by_shape_cr.root" "SR" "" "" "yes"

run_fit "shape_crsim_cr3l_emucr" "shape" "3l CR + emu CR simultaneous data fit" \
  "shapes/shape_ptmiss.root" "CR3L EMUCR" "k_WZ k_emu" "k_WZ,k_emu" "no"
run_fit "shape_crsim_cr3l_dycr" "shape" "3l CR + Z+jet CR simultaneous data fit" \
  "shapes/shape_ptmiss.root" "CR3L DYCR" "k_WZ k_Zjet" "k_WZ,k_Zjet" "no"
run_fit "shape_crsim_emucr_dycr" "shape" "emu CR + Z+jet CR simultaneous data fit" \
  "shapes/shape_ptmiss.root" "EMUCR DYCR" "k_emu k_Zjet" "k_emu,k_Zjet" "no"
run_fit "shape_crsim_allcr" "shape" "3l CR + Z+jet CR + emu CR simultaneous data fit" \
  "shapes/shape_ptmiss.root" "CR3L DYCR EMUCR" "k_WZ k_Zjet k_emu" "k_WZ,k_Zjet,k_emu" "no"

echo "[6/7] Counting fits"
run_fit "counting_cronly_dycr" "counting" "CR only data counting fit: Z+jet CR -> k_Zjet" \
  "shapes/counting_onebin.root" "DYCR" "k_Zjet" "k_Zjet" "no"
run_fit "counting_cronly_cr3l" "counting" "CR only data counting fit: 3l CR -> k_WZ" \
  "shapes/counting_onebin.root" "CR3L" "k_WZ" "k_WZ" "no"
run_fit "counting_cronly_emucr" "counting" "CR only data counting fit: emu CR -> k_emu" \
  "shapes/counting_onebin.root" "EMUCR" "k_emu" "k_emu" "no"

K_COUNT_ZJET="$(python3 scripts/read_fit_param.py --fit fits/counting_cronly_dycr/fitDiagnostics.counting_cronly_dycr.root --param k_Zjet)"
K_COUNT_WZ="$(python3 scripts/read_fit_param.py --fit fits/counting_cronly_cr3l/fitDiagnostics.counting_cronly_cr3l.root --param k_WZ)"
K_COUNT_EMU="$(python3 scripts/read_fit_param.py --fit fits/counting_cronly_emucr/fitDiagnostics.counting_cronly_emucr.root --param k_emu)"

python3 scripts/scale_fit_check_templates.py \
  --input shapes/counting_onebin.root \
  --output shapes/counting_onebin_sr_scaled_by_counting_cr.root \
  --channels SR \
  --k-zjet "${K_COUNT_ZJET}" \
  --k-wz "${K_COUNT_WZ}" \
  --k-emu "${K_COUNT_EMU}" \
  > logs/counting_sr_scale_templates.log 2>&1

echo "CR-only one-bin counting,counting,${K_COUNT_ZJET},${K_COUNT_WZ},${K_COUNT_EMU}" >> tables/k_factors_used.csv

run_fit "counting_sr_data_scaledk" "counting" "SR only data counting fit, templates scaled by CR-only k-factors, no k-factors in card" \
  "shapes/counting_onebin_sr_scaled_by_counting_cr.root" "SR" "" "" "yes"
run_fit "counting_sr_asimov_cr_data_sim" "counting" "SR Asimov + CR data simultaneous counting fit" \
  "shapes/counting_onebin_sr_asimov.root" "SR DYCR CR3L EMUCR" "k_Zjet k_WZ k_emu" "" "yes"
run_fit "counting_sr_cr_data_sim" "counting" "SR+CR data simultaneous counting fit" \
  "shapes/counting_onebin.root" "SR DYCR CR3L EMUCR" "k_Zjet k_WZ k_emu" "" "yes"

echo "[7/7] Done"
python3 scripts/generate_counting_yield_tables.py

echo "Tables:"
echo "  ${WORKDIR}/tables/fit_parameters.csv"
echo "  ${WORKDIR}/tables/fit_status.csv"
echo "  ${WORKDIR}/tables/yield_summary.csv"
echo "  ${WORKDIR}/tables/k_factors_used.csv"
echo "  ${WORKDIR}/tables/manifest.tsv"
echo "  ${WORKDIR}/tables/counting_yields_before_after.md"
