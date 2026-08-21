#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import ROOT
import uproot

BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]

CHANNELS = [
    ("SR", "SR", "lepton_cat!=2 && ptmiss>=100", True),
    ("DYCR", "DY_CR", "lepton_cat!=2 && ptmiss>=30 && ptmiss<=90", False),
    ("EMUCR", "emu_CR", "lepton_cat==2 && ptmiss>=100", False),
    ("CR3L", "3l_CR", "ptmiss>=30", False),
]

DEFAULT_CONFIG_PATH = Path("/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/dnn_config.json")


def load_config(path):
    return json.loads(Path(path).read_text())


def load_metadata(path):
    return json.loads(Path(path).read_text())


def feature_names(feature_specs):
    return [spec["name"] for spec in feature_specs]


def tree_branch_names(path):
    with uproot.open(path)["Vars"] as tree:
        return set(tree.keys()), tree.num_entries


def check_scored_files(input_base, output_base, signal_file):
    rows = []
    for _, source_dir, _, has_signal in CHANNELS:
        filenames = ["Data.root"] + [f"{p}.root" for p in BKG_PROCESSES]
        if has_signal:
            filenames.append(signal_file)
        for filename in filenames:
            source = os.path.join(input_base, source_dir, filename)
            output = os.path.join(output_base, source_dir, "signal.root" if filename == signal_file else filename)
            source_branches, source_entries = tree_branch_names(source)
            output_branches, output_entries = tree_branch_names(output)
            missing = sorted(source_branches - output_branches)
            added = sorted(output_branches - source_branches)
            if missing:
                raise RuntimeError(f"Missing branches in {output}: {missing}")
            if added != ["dnn_score"]:
                raise RuntimeError(f"Unexpected added branches in {output}: {added}")
            if source_entries != output_entries:
                raise RuntimeError(f"Entry mismatch in {output}: {source_entries} vs {output_entries}")

            with uproot.open(output)["Vars"] as tree:
                scores = tree["dnn_score"].array(library="np")
            if not np.all(np.isfinite(scores)):
                raise RuntimeError(f"Non-finite dnn_score in {output}")
            if len(scores) and (scores.min() < -1e-6 or scores.max() > 1.0 + 1e-6):
                raise RuntimeError(f"dnn_score outside [0,1] in {output}: {scores.min()}..{scores.max()}")
            rows.append({
                "file": output,
                "entries": int(output_entries),
                "score_min": float(scores.min()) if len(scores) else None,
                "score_mean": float(scores.mean()) if len(scores) else None,
                "score_max": float(scores.max()) if len(scores) else None,
            })
    return rows


def check_datacard(card_path):
    text = Path(card_path).read_text()
    required = [
        "k_Zjet  rateParam *     DY",
        "k_WZ    rateParam *     WZ",
        "k_emu   rateParam *     ttbar",
        "k_emu   rateParam *     WW",
        "k_emu   rateParam *     ST",
    ]
    missing = [line for line in required if line not in text]
    if missing:
        raise RuntimeError(f"Datacard is missing expected rateParam lines: {missing}")
    if "ptmiss_significance" in text:
        raise RuntimeError("Datacard unexpectedly mentions ptmiss_significance")


def check_workspace_dataset(workspace_path, dataset_name):
    root_file = ROOT.TFile.Open(str(workspace_path))
    if not root_file or root_file.IsZombie():
        raise OSError(workspace_path)
    workspace = root_file.Get("w")
    if not workspace:
        raise RuntimeError(f"Missing workspace w in {workspace_path}")
    dataset = workspace.data(dataset_name)
    if not dataset:
        raise RuntimeError(f"Missing dataset {dataset_name} in {workspace_path}")
    root_file.Close()


def main():
    parser = argparse.ArgumentParser(description="Run DNN workflow sanity checks.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--metadata", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/output/dnn_metadata.json")
    parser.add_argument("--card", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/cards/combined_dnn_score.txt")
    parser.add_argument("--asimov-workspace", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/workspaces/workspace_dnn_score_asimov.root")
    parser.add_argument("--out", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/output/sanity_checks.json")
    args = parser.parse_args()

    config = load_config(args.config)
    metadata = load_metadata(args.metadata)
    names = feature_names(metadata["feature_specs"])
    if "ptmiss_significance" in names:
        raise RuntimeError("ptmiss_significance is still present in DNN input features")
    if "weight" in names:
        raise RuntimeError("weight is present in DNN input features")

    file_rows = check_scored_files(config["input_base"], config["scored_output_base"], config["signal_file"])
    check_datacard(args.card)
    check_workspace_dataset(args.asimov_workspace, "asimovData_1")

    payload = {
        "status": "ok",
        "num_features": len(names),
        "features": names,
        "skipped_branches": metadata.get("skipped_branches", []),
        "checked_files": file_rows,
        "card": args.card,
        "asimov_workspace": args.asimov_workspace,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(json.dumps({"status": "ok", "num_features": len(names), "checked_files": len(file_rows)}, indent=2))


if __name__ == "__main__":
    main()
