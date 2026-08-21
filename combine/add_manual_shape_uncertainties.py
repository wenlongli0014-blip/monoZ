#!/usr/bin/env python3

import argparse

import ROOT


UNCERTAINTIES = {
    "manual_DY_40": (["DY"], 0.40),
    "manual_ZZ_15": (["ggZZ", "qqZZ"], 0.15),
    "manual_WZ_10": (["WZ"], 0.10),
    "manual_ttbar_15": (["ttbar"], 0.15),
    "manual_WW_15": (["WW"], 0.15),
    "manual_ST_15": (["ST"], 0.15),
}


def clone_scaled(hist, name, scale):
    out = hist.Clone(name)
    out.Scale(scale)
    out.SetDirectory(ROOT.nullptr)
    return out


def parse_args():
    parser = argparse.ArgumentParser(description="Add manual Combine shape Up/Down templates to a shapes ROOT file.")
    parser.add_argument("--shapes", required=True, help="Shapes ROOT file to update in place")
    parser.add_argument("--channels", nargs="+", default=["SR", "DYCR", "EMUCR", "CR3L"])
    return parser.parse_args()


def main():
    args = parse_args()
    root_file = ROOT.TFile(args.shapes, "UPDATE")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {args.shapes}")

    for channel in args.channels:
        directory = root_file.GetDirectory(channel)
        if not directory:
            raise KeyError(f"Missing channel directory {channel} in {args.shapes}")
        directory.cd()
        for nuisance, (processes, fraction) in UNCERTAINTIES.items():
            for process in processes:
                hist = directory.Get(process)
                if not hist:
                    continue
                up = clone_scaled(hist, f"{process}_{nuisance}Up", 1.0 + fraction)
                down = clone_scaled(hist, f"{process}_{nuisance}Down", 1.0 - fraction)
                up.Write(up.GetName(), ROOT.TObject.kOverwrite)
                down.Write(down.GetName(), ROOT.TObject.kOverwrite)
                print(f"{channel:5s} {process:6s} {nuisance:16s}: down={1.0-fraction:.3f}, up={1.0+fraction:.3f}")

    root_file.Close()
    print(f"Updated {args.shapes}")


if __name__ == "__main__":
    main()
