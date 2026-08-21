#!/usr/bin/env python3

import argparse
import os
from array import array

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.ROOT.EnableImplicitMT()

INPUT_BASE = "/eos/user/l/liwe/monoZ_combine"
SIGNAL_FILE = "signal.root"

BACKGROUNDS = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]

CHANNELS = [
    {
        "name": "SR",
        "source_dir": "SR",
        "selection": "lepton_cat!=2 && ptmiss>=100",
        "has_signal": True,
        "ptmiss_bins": [100, 140, 200, 320, 550, 3000],
    },
    {
        "name": "DYCR",
        "source_dir": "DY_CR",
        "selection": "lepton_cat!=2 && ptmiss>=30 && ptmiss<=90",
        "has_signal": False,
        "ptmiss_bins": [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90],
    },
    {
        "name": "EMUCR",
        "source_dir": "emu_CR",
        "selection": "lepton_cat==2 && ptmiss>=100",
        "has_signal": False,
        "ptmiss_bins": [100, 140, 180, 1100],
    },
    {
        "name": "CR3L",
        "source_dir": "3l_CR",
        "selection": "1",
        "has_signal": False,
        "ptmiss_bins": [30, 70, 120, 220, 420, 2000],
    },
]


def build_hist(file_path, selection, mode, bins, hist_name, weight_expr=None):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Missing input file: {file_path}")

    rdf = ROOT.RDataFrame("Vars", file_path)
    if selection and selection != "1":
        rdf = rdf.Filter(selection)

    if mode == "shape":
        variable = "ptmiss"
        rdf = rdf.Filter("std::isfinite(ptmiss)")
    elif mode == "counting":
        variable = "fit_check_onebin"
        rdf = rdf.Define(variable, "0.5")
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    model = ROOT.RDF.TH1DModel(hist_name, "", len(bins) - 1, array("d", bins))
    if weight_expr:
        rdf = rdf.Filter(f"std::isfinite({weight_expr})")
        rdf = rdf.Define("fit_check_weight", weight_expr)
        hist = rdf.Histo1D(model, variable, "fit_check_weight").GetValue().Clone(hist_name)
    else:
        hist = rdf.Histo1D(model, variable).GetValue().Clone(hist_name)

    hist.SetDirectory(ROOT.nullptr)
    return hist


def floor_mc_hist(hist, context):
    for i_bin in range(1, hist.GetNbinsX() + 1):
        value = hist.GetBinContent(i_bin)
        error = hist.GetBinError(i_bin)
        if value <= 0:
            print(f"  [WARN] {context}: bin {i_bin} yield {value:.6g}; flooring to 1e-6")
            hist.SetBinContent(i_bin, 1e-6)
        if error <= 0:
            hist.SetBinError(i_bin, 1e-6)


def make_asimov(directory):
    total = None
    for process in BACKGROUNDS:
        hist = directory.Get(process)
        if not hist:
            continue
        if total is None:
            total = hist.Clone("data_obs")
            total.Reset("ICES")
            total.SetDirectory(ROOT.nullptr)
        total.Add(hist)
    if total is None:
        raise RuntimeError(f"No background histograms found in {directory.GetPath()}")
    return total


def make_dummy_signal(template):
    dummy = template.Clone("signal")
    dummy.Reset("ICES")
    dummy.SetDirectory(ROOT.nullptr)
    for i_bin in range(1, dummy.GetNbinsX() + 1):
        dummy.SetBinContent(i_bin, 1e-9)
        dummy.SetBinError(i_bin, 1e-12)
    return dummy


def parse_args():
    parser = argparse.ArgumentParser(description="Build fit-check shapes from monoZ ROOT ntuples.")
    parser.add_argument("--mode", choices=["shape", "counting"], required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--input-base", default=INPUT_BASE)
    parser.add_argument("--signal-file", default=SIGNAL_FILE)
    parser.add_argument(
        "--asimov-channels",
        nargs="*",
        default=[],
        choices=[channel["name"] for channel in CHANNELS],
        help="Channels whose data_obs should be replaced by sum(backgrounds).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    asimov_channels = set(args.asimov_channels)

    output = ROOT.TFile(args.out, "RECREATE")
    if not output or output.IsZombie():
        raise OSError(f"Could not create {args.out}")

    print(f"Writing {args.out}")
    print(f"Input base = {args.input_base}")
    print(f"Mode = {args.mode}")
    print(f"Signal file = {args.signal_file}")
    if asimov_channels:
        print(f"Asimov channels = {', '.join(sorted(asimov_channels))}")

    for channel in CHANNELS:
        name = channel["name"]
        source_dir = channel["source_dir"]
        selection = channel["selection"]
        bins = channel["ptmiss_bins"] if args.mode == "shape" else [0.0, 1.0]

        print(f"\n[{name}] source={source_dir} selection=({selection}) bins={bins}")
        output.mkdir(name)
        directory = output.GetDirectory(name)
        directory.cd()

        data_path = os.path.join(args.input_base, source_dir, "Data.root")
        data_hist = build_hist(data_path, selection, args.mode, bins, "data_obs")
        data_hist.Write("data_obs")
        print(f"  data_obs: {data_hist.Integral(1, data_hist.GetNbinsX()):10.3f}")

        for process in BACKGROUNDS:
            process_path = os.path.join(args.input_base, source_dir, f"{process}.root")
            hist = build_hist(process_path, selection, args.mode, bins, process, "weight")
            floor_mc_hist(hist, f"{name}/{process}")
            hist.Write(process)
            print(f"  {process:6s}: {hist.Integral(1, hist.GetNbinsX()):10.3f}")

        if channel["has_signal"]:
            signal_path = os.path.join(args.input_base, source_dir, args.signal_file)
            signal_hist = build_hist(signal_path, selection, args.mode, bins, "signal", "weight")
            floor_mc_hist(signal_hist, f"{name}/signal")
            signal_hist.Write("signal")
            print(f"  {'signal':6s}: {signal_hist.Integral(1, signal_hist.GetNbinsX()):10.3f}")
        else:
            # Combine requires at least one signal process in a datacard. CR-only
            # cards include this negligible template only as a technical placeholder.
            dummy_signal = make_dummy_signal(data_hist)
            dummy_signal.Write("signal")
            print(f"  {'signal':6s}: {dummy_signal.Integral(1, dummy_signal.GetNbinsX()):10.3e} (dummy)")

        if name in asimov_channels:
            asimov = make_asimov(directory)
            asimov.Write("data_obs", ROOT.TObject.kOverwrite)
            print(f"  [ASIMOV] data_obs -> sum(backgrounds) = {asimov.Integral(1, asimov.GetNbinsX()):.6f}")

    output.Close()
    print("\nDone.")


if __name__ == "__main__":
    main()
