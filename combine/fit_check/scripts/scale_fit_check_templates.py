#!/usr/bin/env python3

import argparse
import os
import shutil

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

SCALE_MAP = {
    "DY": "k_zjet",
    "WZ": "k_wz",
    "ttbar": "k_emu",
    "WW": "k_emu",
    "ST": "k_emu",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Copy a shapes file and scale SR templates by fixed k-factors.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--channels", nargs="+", default=["SR"])
    parser.add_argument("--k-zjet", type=float, required=True)
    parser.add_argument("--k-wz", type=float, required=True)
    parser.add_argument("--k-emu", type=float, required=True)
    return parser.parse_args()


def scale_hist(directory, process, factor):
    hist = directory.Get(process)
    if not hist:
        raise KeyError(f"Missing template {directory.GetPath()}/{process}")
    before = hist.Integral(1, hist.GetNbinsX())
    hist.Scale(factor)
    hist.Write(process, ROOT.TObject.kOverwrite)
    after = hist.Integral(1, hist.GetNbinsX())
    return before, after


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    shutil.copyfile(args.input, args.output)

    factors = {
        "k_zjet": args.k_zjet,
        "k_wz": args.k_wz,
        "k_emu": args.k_emu,
    }

    root_file = ROOT.TFile(args.output, "UPDATE")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {args.output}")

    print(f"Input  = {args.input}")
    print(f"Output = {args.output}")
    for channel in args.channels:
        directory = root_file.GetDirectory(channel)
        if not directory:
            raise KeyError(f"Missing channel {channel} in {args.output}")
        directory.cd()
        for process, factor_name in SCALE_MAP.items():
            before, after = scale_hist(directory, process, factors[factor_name])
            print(f"{channel:5s} {process:6s} {factor_name:7s}: {before:.6f} -> {after:.6f}")

    root_file.Close()


if __name__ == "__main__":
    main()
