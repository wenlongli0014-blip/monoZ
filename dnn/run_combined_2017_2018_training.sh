#!/usr/bin/env bash
set -euo pipefail

DNN_DIR="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn"
CONFIG="${DNN_DIR}/combined_2017_2018/dnn_config_combined_2017_2018.json"
OUTDIR="${DNN_DIR}/combined_2017_2018/output"
LOGDIR="${DNN_DIR}/combined_2017_2018/logs"
LCG_VIEW="/cvmfs/sft.cern.ch/lcg/views/LCG_109_cuda/x86_64-el9-gcc13-opt/setup.sh"

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${DNN_DIR}"

set +u
source "${LCG_VIEW}"
set -u

python3 train_dnn_score.py \
  --config "${CONFIG}" \
  --outdir "${OUTDIR}" \
  --epochs 50 \
  --patience 10 \
  --max-events-per-process 0 \
  --torch-threads 4 \
  | tee "${LOGDIR}/01_train_combined.log"
