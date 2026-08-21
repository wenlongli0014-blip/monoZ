#!/usr/bin/env python3

import argparse

import ROOT


BACKGROUNDS = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]


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


def parse_args():
    parser = argparse.ArgumentParser(description="Replace data_obs with background-only pseudo-Asimov data.")
    parser.add_argument("--shapes", required=True, help="Shapes ROOT file to update in place")
    parser.add_argument(
        "--asimov-channels",
        nargs="+",
        required=True,
        help="Channels whose data_obs should be replaced by sum(backgrounds), e.g. SR DYCR EMUCR CR3L",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    root_file = ROOT.TFile(args.shapes, "UPDATE")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {args.shapes}")

    for channel in args.asimov_channels:
        directory = root_file.GetDirectory(channel)
        if not directory:
            raise KeyError(f"Missing channel directory {channel} in {args.shapes}")
        directory.cd()
        asimov = make_asimov(directory)
        asimov.Write("data_obs", ROOT.TObject.kOverwrite)
        print(f"{channel:5s}: replaced data_obs with b-only Asimov, yield={asimov.Integral():.6f}")

    root_file.Close()
    print(f"Updated {args.shapes}")


if __name__ == "__main__":
    main()
