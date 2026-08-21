#!/usr/bin/env python3

import argparse

import ROOT


def parse_args():
    parser = argparse.ArgumentParser(description="Scale process templates in a shapes ROOT file in place.")
    parser.add_argument("--shapes", required=True)
    parser.add_argument("--process", required=True)
    parser.add_argument("--scale", type=float, required=True)
    parser.add_argument("--channels", nargs="+", default=["SR", "DYCR", "EMUCR", "CR3L"])
    parser.add_argument(
        "--include-variations",
        action="store_true",
        help="Also scale histograms whose names start with '<process>_'.",
    )
    return parser.parse_args()


def should_scale(name, process, include_variations):
    if name == process:
        return True
    return include_variations and name.startswith(f"{process}_")


def main():
    args = parse_args()
    root_file = ROOT.TFile(args.shapes, "UPDATE")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {args.shapes}")

    for channel in args.channels:
        directory = root_file.GetDirectory(channel)
        if not directory:
            raise KeyError(f"Missing channel {channel} in {args.shapes}")
        directory.cd()
        keys = [key.GetName() for key in directory.GetListOfKeys()]
        for name in sorted(set(keys)):
            if not should_scale(name, args.process, args.include_variations):
                continue
            hist = directory.Get(name)
            if not hist:
                continue
            before = hist.Integral(1, hist.GetNbinsX())
            hist.Scale(args.scale)
            hist.Write(name, ROOT.TObject.kOverwrite)
            after = hist.Integral(1, hist.GetNbinsX())
            print(f"{channel:5s} {name:28s}: {before:.6f} -> {after:.6f}")

    root_file.Close()
    print(f"Updated {args.shapes}")


if __name__ == "__main__":
    main()
