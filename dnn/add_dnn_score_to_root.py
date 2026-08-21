#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
from array import array
from pathlib import Path

import numpy as np
import ROOT
import torch
import uproot

from dnn_common import (
    BKG_PROCESSES,
    CHANNELS,
    MISSING_SENTINEL,
    available_branches,
    branches_for_specs,
    extract_features_from_arrays,
    feature_names,
    load_metadata,
    load_model,
    predict_scores,
    transform_features,
)


def score_file(input_path, model, mean, scale, feature_specs, metadata, device, batch_size):
    requested = branches_for_specs(feature_specs)
    with uproot.open(input_path)["Vars"] as tree:
        arrays = tree.arrays(available_branches(tree, requested), library="ak")
        n_entries = int(tree.num_entries)
    raw_x = extract_features_from_arrays(arrays, feature_specs, length=n_entries).astype(np.float32)
    x = transform_features(
        raw_x,
        mean,
        scale,
        missing_sentinel=float(metadata.get("missing_sentinel", MISSING_SENTINEL)),
        preprocessing=str(metadata.get("preprocessing", "standardize_nonmissing_then_map_missing_to_zero")),
    )
    scores = predict_scores(model, x, device, batch_size=batch_size)
    return scores.astype(np.float32)


def clone_with_score(input_path, output_path, scores):
    in_file = ROOT.TFile.Open(input_path, "READ")
    if not in_file or in_file.IsZombie():
        raise OSError(f"Could not open {input_path}")
    in_tree = in_file.Get("Vars")
    if not in_tree:
        raise KeyError(f"Missing Vars tree in {input_path}")
    if in_tree.GetEntries() != len(scores):
        raise RuntimeError(f"Score length mismatch for {input_path}: {len(scores)} vs {in_tree.GetEntries()}")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    out_file = ROOT.TFile(output_path, "RECREATE")
    out_tree = in_tree.CloneTree(0)
    score_buf = array("f", [0.0])
    out_tree.Branch("dnn_score", score_buf, "dnn_score/F")

    for i_entry, score in enumerate(scores):
        in_tree.GetEntry(i_entry)
        score_buf[0] = float(score)
        out_tree.Fill()

    out_tree.Write("Vars")
    out_file.Close()
    in_file.Close()


def iter_input_files(input_base, signal_file):
    for _, source_dir, _, has_signal in CHANNELS:
        region_dir = os.path.join(input_base, source_dir)
        yield source_dir, "Data.root", os.path.join(region_dir, "Data.root")
        for process in BKG_PROCESSES:
            yield source_dir, f"{process}.root", os.path.join(region_dir, f"{process}.root")
        if has_signal:
            yield source_dir, signal_file, os.path.join(region_dir, signal_file)


def parse_args():
    parser = argparse.ArgumentParser(description="Copy monoZ ROOT trees and add a dnn_score branch.")
    parser.add_argument("--input-base", default="/eos/user/l/liwe/monoZ_combine")
    parser.add_argument("--output-base", default="/eos/user/l/liwe/monoZ_combine_dnn")
    parser.add_argument("--model", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/output/dnn_model.pt")
    parser.add_argument("--metadata", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/output/dnn_metadata.json")
    parser.add_argument("--signal-file", default="signal.root")
    parser.add_argument("--batch-size", type=int, default=65536)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    metadata, mean, scale = load_metadata(args.metadata)
    feature_specs = metadata["feature_specs"]
    print(f"Loaded {len(feature_specs)} DNN input features from metadata.")
    print("Features:", ", ".join(feature_names(feature_specs)))

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = load_model(args.model, device)
    print(f"Using device: {device}")
    print(f"Input base : {args.input_base}")
    print(f"Output base: {args.output_base}")

    for source_dir, filename, input_path in iter_input_files(args.input_base, args.signal_file):
        output_name = "signal.root" if filename == args.signal_file else filename
        output_path = os.path.join(args.output_base, source_dir, output_name)
        if os.path.exists(output_path) and not args.overwrite:
            print(f"[SKIP] {output_path} exists")
            continue
        if not os.path.exists(input_path):
            raise FileNotFoundError(input_path)

        scores = score_file(input_path, model, mean, scale, feature_specs, metadata, device, args.batch_size)
        clone_with_score(input_path, output_path, scores)
        if len(scores):
            score_summary = f"({scores.min():.4f},{scores.mean():.4f},{scores.max():.4f})"
        else:
            score_summary = "(empty,empty,empty)"
        print(
            f"[OK] {source_dir}/{filename:12s} -> {output_path} "
            f"entries={len(scores)} score[min,mean,max]={score_summary}"
        )


if __name__ == "__main__":
    main()
