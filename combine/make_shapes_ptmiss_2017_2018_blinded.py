#!/usr/bin/env python3

"""Build combined multi-year ptmiss templates with a blinded SR."""

import argparse
import os
from array import array

import numpy as np
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

BACKGROUNDS = [
    "DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"
]

CHANNELS = [
    {
        "name": "SR",
        "source_dir": "SR",
        "bins": [100, 140, 200, 320, 550, 3000],
        "selection": "sr",
        "has_signal": True,
        "asimov": True,
    },
    {
        "name": "DYCR",
        "source_dir": "DY_CR",
        "bins": [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90],
        "selection": "dy",
        "has_signal": False,
        "asimov": False,
    },
    {
        "name": "EMUCR",
        "source_dir": "emu_CR",
        "bins": [100, 140, 180, 1100],
        "selection": "emu",
        "has_signal": False,
        "asimov": False,
    },
    {
        "name": "CR3L",
        "source_dir": "3l_CR",
        "bins": [30, 70, 120, 220, 420, 2000],
        "selection": "all",
        "has_signal": False,
        "asimov": False,
    },
]


def fill_arrays(paths, bins, selection, weighted, mode):
    chain = ROOT.TChain("Vars")
    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        if chain.Add(path) != 1:
            raise OSError(f"Could not add {path} to TChain")

    selections = {
        "sr": "lepton_cat != 2 && ptmiss >= 100",
        "dy": "lepton_cat != 2 && ptmiss >= 30 && ptmiss <= 90",
        "emu": "lepton_cat == 2 && ptmiss >= 100",
        "all": "true",
    }
    expression = f"std::isfinite(ptmiss) && ({selections[selection]})"
    if weighted:
        expression += " && std::isfinite(weight)"

    frame = ROOT.RDataFrame(chain).Filter(expression)
    observable = "ptmiss"
    if mode == "counting":
        frame = frame.Define("tau_study_counting_observable", "0.5")
        observable = "tau_study_counting_observable"

    model = ROOT.RDF.TH1DModel(
        "", "", len(bins) - 1, array("d", bins)
    )
    if weighted:
        result = frame.Histo1D(model, observable, "weight")
    else:
        result = frame.Histo1D(model, observable)
    hist = result.GetValue()

    contents = np.asarray([
        hist.GetBinContent(index)
        for index in range(1, hist.GetNbinsX() + 1)
    ], dtype=np.float64)
    errors = np.asarray([
        hist.GetBinError(index)
        for index in range(1, hist.GetNbinsX() + 1)
    ], dtype=np.float64)
    return contents, errors


def make_hist(name, bins, contents, errors):
    hist = ROOT.TH1D(name, name, len(bins) - 1, array("d", bins))
    hist.Sumw2()
    hist.SetDirectory(ROOT.nullptr)
    for index, (value, error) in enumerate(zip(contents, errors), start=1):
        hist.SetBinContent(index, float(value))
        hist.SetBinError(index, float(error))
    return hist


def floor_mc(contents, errors, context):
    for index in range(len(contents)):
        if contents[index] <= 0:
            print(
                f"  [WARN] {context}: bin {index + 1} yield "
                f"{contents[index]:.6g}; flooring to 1e-6"
            )
            contents[index] = 1e-6
        if errors[index] <= 0:
            errors[index] = 1e-6


def paths_for(input_base, years, source_dir, sample):
    return [
        os.path.join(input_base, year, source_dir, f"{sample}.root")
        for year in years
    ]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-base",
        default="/eos/user/l/liwe/monoZ_combine_17_18_dnn",
    )
    parser.add_argument(
        "--years",
        nargs="+",
        default=["2017", "2018"],
        choices=["2016", "2017", "2018"],
        help="Years to combine (default: 2017 2018).",
    )
    parser.add_argument(
        "--sr-bins",
        nargs="+",
        type=float,
        default=None,
        help=(
            "Optional SR ptmiss bin edges in GeV. The default is "
            "100 140 200 320 550 3000."
        ),
    )
    parser.add_argument("--mode", choices=["shape", "counting"], default="shape")
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.sr_bins is not None:
        if len(args.sr_bins) < 2:
            raise ValueError("--sr-bins needs at least two bin edges")
        if any(right <= left for left, right in zip(args.sr_bins, args.sr_bins[1:])):
            raise ValueError("--sr-bins must be strictly increasing")
        if args.sr_bins[0] < 100:
            raise ValueError("The SR selection starts at 100 GeV")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    output = ROOT.TFile(args.out, "RECREATE")
    if not output or output.IsZombie():
        raise OSError(f"Could not create {args.out}")

    print(f"Input base = {args.input_base}")
    print(f"Years = {' + '.join(args.years)}")
    print(f"Mode = {args.mode}")
    print("SR data_obs = background-only Asimov (SR/Data.root is not read)")
    print("EMUCR selection = lepton_cat == 2 && ptmiss >= 100 GeV")

    for channel in CHANNELS:
        name = channel["name"]
        channel_bins = (
            args.sr_bins
            if name == "SR" and args.sr_bins is not None
            else channel["bins"]
        )
        bins = np.asarray(
            channel_bins if args.mode == "shape" else [0.0, 1.0],
            dtype=np.float64,
        )
        output.mkdir(name)
        directory = output.GetDirectory(name)
        directory.cd()
        print(f"\n[{name}] bins = {list(bins)}")

        background_hists = []
        for process in BACKGROUNDS:
            paths = paths_for(
                args.input_base, args.years, channel["source_dir"], process)
            contents, errors = fill_arrays(
                paths, bins, channel["selection"], weighted=True,
                mode=args.mode)
            floor_mc(contents, errors, f"{name}/{process}")
            hist = make_hist(process, bins, contents, errors)
            hist.Write(process)
            background_hists.append(hist)
            print(f"  {process:6s}: {contents.sum():12.4f}")

        if channel["asimov"]:
            data_hist = background_hists[0].Clone("data_obs")
            data_hist.Reset("ICES")
            data_hist.SetDirectory(ROOT.nullptr)
            for hist in background_hists:
                data_hist.Add(hist)
            print(
                "  data_obs: "
                f"{data_hist.Integral(1, data_hist.GetNbinsX()):12.4f} "
                "(background-only Asimov)"
            )
        else:
            data_paths = paths_for(
                args.input_base, args.years, channel["source_dir"], "Data")
            contents, errors = fill_arrays(
                data_paths, bins, channel["selection"], weighted=False,
                mode=args.mode)
            data_hist = make_hist("data_obs", bins, contents, errors)
            print(f"  data_obs: {contents.sum():12.4f}")
        data_hist.Write("data_obs")

        if channel["has_signal"]:
            signal_paths = paths_for(
                args.input_base, args.years, channel["source_dir"], "signal")
            contents, errors = fill_arrays(
                signal_paths, bins, channel["selection"], weighted=True,
                mode=args.mode)
            floor_mc(contents, errors, f"{name}/signal")
            signal_hist = make_hist("signal", bins, contents, errors)
            signal_hist.Write("signal")
            print(f"  signal: {contents.sum():12.4f}")
        else:
            dummy_contents = np.full(len(bins) - 1, 1e-9)
            dummy_errors = np.full(len(bins) - 1, 1e-12)
            dummy_signal = make_hist(
                "signal", bins, dummy_contents, dummy_errors)
            dummy_signal.Write("signal")

    output.Close()
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
