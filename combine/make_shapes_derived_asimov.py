#!/usr/bin/env python3

import argparse
import os
from array import array

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare("#include <TVector2.h>")


INPUT_BASE = "/eos/user/l/liwe/monoZ_combine"
SIGNAL_FILE = "signal.root"
BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]

VARIABLES = {
    "deltaR_ll": "std::sqrt(std::pow(lepton_eta[0]-lepton_eta[1],2)+std::pow(TVector2::Phi_mpi_pi(lepton_phi[0]-lepton_phi[1]),2))",
    "dphi_met_ll": "std::abs(TVector2::Phi_mpi_pi(ll_phi-ptmiss_phi))",
}

CHANNELS = [
    {
        "output_name": "SR",
        "source_dir": "SR",
        "selection": "lepton_cat!=2 && ptmiss>=100",
        "has_signal": True,
        "bins": {
            "deltaR_ll": [0.0, 0.4, 0.7, 1.0, 1.2, 1.4, 1.6, 1.8],
            "dphi_met_ll": [2.4, 2.55, 2.7, 2.85, 3.0, 3.1416],
        },
    },
    {
        "output_name": "DYCR",
        "source_dir": "DY_CR",
        "selection": "lepton_cat!=2 && ptmiss>=30 && ptmiss<=90",
        "has_signal": False,
        "bins": {
            "deltaR_ll": [0.0, 0.4, 0.7, 1.0, 1.2, 1.4, 1.6, 1.8],
            "dphi_met_ll": [2.4, 2.55, 2.7, 2.85, 3.0, 3.1416],
        },
    },
    {
        "output_name": "EMUCR",
        "source_dir": "emu_CR",
        "selection": "lepton_cat==2",
        "has_signal": False,
        "bins": {
            "deltaR_ll": [0.0, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8],
            "dphi_met_ll": [2.4, 2.55, 2.7, 2.85, 3.0, 3.1416],
        },
    },
    {
        "output_name": "CR3L",
        "source_dir": "3l_CR",
        "selection": "1",
        "has_signal": False,
        "bins": {
            "deltaR_ll": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0],
            "dphi_met_ll": [0.0, 0.6, 1.2, 1.8, 2.4, 2.7, 3.1416],
        },
    },
]


def build_hist(file_path, selection, variable, expression, bins, hist_name, weight_expr=None):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    rdf = ROOT.RDataFrame("Vars", file_path)
    if selection and selection != "1":
        rdf = rdf.Filter(selection)
    rdf = rdf.Define(variable, expression).Filter(f"std::isfinite({variable})")

    model = ROOT.RDF.TH1DModel(hist_name, "", len(bins) - 1, array("d", bins))
    if weight_expr:
        rdf = rdf.Filter(f"std::isfinite({weight_expr})")
        rdf = rdf.Define("w_for_hist", weight_expr)
        hist = rdf.Histo1D(model, variable, "w_for_hist").GetValue().Clone(hist_name)
    else:
        hist = rdf.Histo1D(model, variable).GetValue().Clone(hist_name)
    hist.SetDirectory(ROOT.nullptr)
    return hist


def floor_mc_hist(hist, context):
    for i_bin in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(i_bin) <= 0:
            print(f"  [WARN] {context}: bin {i_bin} non-positive yield; flooring to 1e-6")
            hist.SetBinContent(i_bin, 1e-6)
        if hist.GetBinError(i_bin) <= 0:
            hist.SetBinError(i_bin, 1e-6)


def integral(hist):
    return float(hist.Integral(1, hist.GetNbinsX()))


def parse_args():
    parser = argparse.ArgumentParser(description="Build shapes for derived dilepton angular variables.")
    parser.add_argument("--var", required=True, choices=sorted(VARIABLES))
    parser.add_argument("--out", required=True)
    parser.add_argument("--input-base", default=INPUT_BASE)
    parser.add_argument("--signal-file", default=SIGNAL_FILE)
    return parser.parse_args()


def main():
    args = parse_args()
    expression = VARIABLES[args.var]
    out = ROOT.TFile(args.out, "RECREATE")
    print(f"Writing {args.out}")
    print(f"Variable = {args.var}")
    print(f"Expression = {expression}")

    for channel in CHANNELS:
        out_name = channel["output_name"]
        src_dir = channel["source_dir"]
        selection = channel["selection"]
        bins = channel["bins"][args.var]

        print(f"\n[{out_name}] source={src_dir} selection=({selection}) bins={bins}")
        out.mkdir(out_name)
        out.cd(out_name)

        data_path = os.path.join(args.input_base, src_dir, "Data.root")
        h_data = build_hist(data_path, selection, args.var, expression, bins, "data_obs")
        h_data.Write("data_obs")
        print(f"  data_obs: {integral(h_data):10.3f}")

        for process in BKG_PROCESSES:
            path = os.path.join(args.input_base, src_dir, f"{process}.root")
            hist = build_hist(path, selection, args.var, expression, bins, process, "weight")
            floor_mc_hist(hist, f"{out_name}/{process}")
            hist.Write(process)
            print(f"  {process:6s}: {integral(hist):10.3f}")

        if channel["has_signal"]:
            signal_path = os.path.join(args.input_base, src_dir, args.signal_file)
            hist = build_hist(signal_path, selection, args.var, expression, bins, "signal", "weight")
            floor_mc_hist(hist, f"{out_name}/signal")
            hist.Write("signal")
            print(f"  {'signal':6s}: {integral(hist):10.3f}")

    out.Close()
    print("\nDone.")


if __name__ == "__main__":
    main()
