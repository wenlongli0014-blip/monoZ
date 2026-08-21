#!/usr/bin/env bash
set -euo pipefail

YEAR="${1:?Usage: run_year_train_score.sh YEAR INPUT_BASE OUTPUT_BASE}"
INPUT_BASE="${2:?Missing input base}"
OUTPUT_BASE="${3:?Missing output base}"

REPO="/afs/cern.ch/user/l/liwe/hzz2l2nu"
DNN_DIR="${REPO}/dnn"
YEAR_DIR="${DNN_DIR}/${YEAR}"
CONFIG="${YEAR_DIR}/dnn_config_${YEAR}.json"
LCG_VIEW="/cvmfs/sft.cern.ch/lcg/views/LCG_109_cuda/x86_64-el9-gcc13-opt/setup.sh"

mkdir -p "${YEAR_DIR}/output" "${YEAR_DIR}/logs" "${OUTPUT_BASE}"

cd "${DNN_DIR}"

echo "[$(date)] year=${YEAR}"
echo "Input base: ${INPUT_BASE}"
echo "Output base: ${OUTPUT_BASE}"
echo "Config: ${CONFIG}"

set +u
source "${LCG_VIEW}"
set -u

echo "[$(date)] Train DNN"
python3 train_dnn_score.py \
  --config "${CONFIG}" \
  --input-base "${INPUT_BASE}" \
  --signal-file signal.root \
  --outdir "${YEAR_DIR}/output" \
  --epochs 50 \
  --patience 10 \
  --max-events-per-process 0 \
  --torch-threads 4 \
  | tee "${YEAR_DIR}/logs/01_train.log"

echo "[$(date)] Add dnn_score to ROOT files"
python3 add_dnn_score_to_root.py \
  --input-base "${INPUT_BASE}" \
  --output-base "${OUTPUT_BASE}" \
  --model "${YEAR_DIR}/output/dnn_model.pt" \
  --metadata "${YEAR_DIR}/output/dnn_metadata.json" \
  --signal-file signal.root \
  --overwrite \
  | tee "${YEAR_DIR}/logs/02_add_scores.log"

echo "[$(date)] Basic scored ROOT sanity check"
python3 - "${INPUT_BASE}" "${OUTPUT_BASE}" "${YEAR_DIR}/output/dnn_metadata.json" "${YEAR_DIR}/output/basic_sanity.json" <<'PY'
import json
import os
import sys

import numpy as np
import uproot

input_base, output_base, metadata_path, out_path = sys.argv[1:]
channels = ["SR", "DY_CR", "emu_CR", "3l_CR"]
processes = ["Data", "DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
rows = []

with open(metadata_path) as handle:
    metadata = json.load(handle)
features = metadata["feature_names"]
for forbidden in ["weight", "ptmiss_significance", "num_pv_good", "ll_phi", "ptmiss_phi"]:
    if forbidden in features:
        raise RuntimeError(f"{forbidden} is present in DNN features")

for channel in channels:
    filenames = [f"{p}.root" for p in processes]
    if channel == "SR":
        filenames.append("signal.root")
    for filename in filenames:
        source = os.path.join(input_base, channel, filename)
        output = os.path.join(output_base, channel, filename)
        with uproot.open(source)["Vars"] as src, uproot.open(output)["Vars"] as out:
            src_branches = set(src.keys())
            out_branches = set(out.keys())
            if src.num_entries != out.num_entries:
                raise RuntimeError(f"Entry mismatch: {output}")
            missing = sorted(src_branches - out_branches)
            added = sorted(out_branches - src_branches)
            if missing:
                raise RuntimeError(f"Missing branches in {output}: {missing}")
            if added != ["dnn_score"]:
                raise RuntimeError(f"Unexpected added branches in {output}: {added}")
            scores = out["dnn_score"].array(library="np")
            if not np.all(np.isfinite(scores)):
                raise RuntimeError(f"Nonfinite dnn_score in {output}")
            if len(scores) and (scores.min() < -1e-6 or scores.max() > 1.0 + 1e-6):
                raise RuntimeError(f"dnn_score outside [0,1] in {output}")
            rows.append({
                "file": output,
                "entries": int(out.num_entries),
                "score_min": float(scores.min()) if len(scores) else None,
                "score_mean": float(scores.mean()) if len(scores) else None,
                "score_max": float(scores.max()) if len(scores) else None,
            })

payload = {
    "status": "ok",
    "num_features": len(features),
    "features": features,
    "checked_files": rows,
}
with open(out_path, "w") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
print(json.dumps({"status": "ok", "num_features": len(features), "checked_files": len(rows)}, indent=2))
PY

echo "[$(date)] Done year=${YEAR}"
echo "Model: ${YEAR_DIR}/output/dnn_model.pt"
echo "Metadata: ${YEAR_DIR}/output/dnn_metadata.json"
echo "Scored ROOTs: ${OUTPUT_BASE}"
