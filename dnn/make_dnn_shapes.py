#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
from array import array

import ROOT

from dnn_common import BKG_PROCESSES, CHANNELS, DEFAULT_CONFIG_PATH, load_config


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.ROOT.EnableImplicitMT()

DEFAULT_BINS = {
    "SR": [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.875, 1.0],
    "DYCR": [0.0, 0.15, 0.30, 1.0],
    "EMUCR": [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 1.0],
    "CR3L": [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 1.0],
}


def build_hist(file_path, selection, bins, hist_name, weight_expr=None):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    rdf = ROOT.RDataFrame("Vars", file_path)
    if selection and selection != "1":
        rdf = rdf.Filter(selection)
    rdf = rdf.Filter("std::isfinite(dnn_score) && dnn_score >= 0 && dnn_score <= 1")
    model = ROOT.RDF.TH1DModel(hist_name, "", len(bins) - 1, array("d", bins))
    if weight_expr:
        rdf = rdf.Filter(f"std::isfinite({weight_expr})")
        rdf = rdf.Define("w_for_hist", weight_expr)
        hist = rdf.Histo1D(model, "dnn_score", "w_for_hist").GetValue().Clone(hist_name)
    else:
        hist = rdf.Histo1D(model, "dnn_score").GetValue().Clone(hist_name)
    hist.SetDirectory(ROOT.nullptr)
    return hist


def floor_mc_hist(hist, context):
    for i_bin in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(i_bin) <= 0:
            print(f"  [WARN] {context}: bin {i_bin} <= 0; flooring to 1e-6")
            hist.SetBinContent(i_bin, 1e-6)
        if hist.GetBinError(i_bin) <= 0:
            hist.SetBinError(i_bin, 1e-6)


def parse_edges(text):
    edges = [float(item) for item in text.split(",") if item.strip()]
    if len(edges) < 2 or any(high <= low for low, high in zip(edges[:-1], edges[1:])):
        raise ValueError(f"Invalid bin edges: {text}")
    return edges


def parse_args():
    parser = argparse.ArgumentParser(description="Build DNN score shape templates for monoZ Combine.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--input-base", default="/eos/user/l/liwe/monoZ_combine_dnn")
    parser.add_argument("--signal-file", default="signal.root")
    parser.add_argument("--out", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/shapes_dnn_score.root")
    parser.add_argument(
        "--bin-edges",
        default=None,
        help="Optional comma-separated DNN score bin edges used for all channels.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    common_bins = parse_edges(args.bin_edges) if args.bin_edges else None
    out = ROOT.TFile(args.out, "RECREATE")
    print(f"Writing {args.out}")
    print(f"Input base = {args.input_base}")
    configured_bins = config.get("binning", DEFAULT_BINS)
    print(f"DNN score bins = {common_bins if common_bins is not None else configured_bins}")

    for output_name, source_dir, selection, has_signal in CHANNELS:
        bins = common_bins if common_bins is not None else configured_bins[output_name]
        print(f"\n[{output_name}] source={source_dir} selection=({selection})")
        print(f"  bins={bins}")
        out.mkdir(output_name)
        out.cd(output_name)

        h_data = build_hist(os.path.join(args.input_base, source_dir, "Data.root"), selection, bins, "data_obs")
        h_data.Write("data_obs")
        print(f"  data_obs integral = {h_data.Integral(1, h_data.GetNbinsX()):.3f}")

        for process in BKG_PROCESSES:
            hist = build_hist(os.path.join(args.input_base, source_dir, f"{process}.root"), selection, bins, process, "weight")
            floor_mc_hist(hist, f"{output_name}/{process}")
            hist.Write(process)
            print(f"  {process:6s}: {hist.Integral(1, hist.GetNbinsX()):10.3f}")

        if has_signal:
            hist = build_hist(os.path.join(args.input_base, source_dir, args.signal_file), selection, bins, "signal", "weight")
            floor_mc_hist(hist, f"{output_name}/signal")
            hist.Write("signal")
            print(f"  {'signal':6s}: {hist.Integral(1, hist.GetNbinsX()):10.3f}")

    out.Close()
    print("\nDone.")


if __name__ == "__main__":
    main()
