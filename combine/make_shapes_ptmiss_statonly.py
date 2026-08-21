#!/usr/bin/env python3

import os
from array import array

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.ROOT.EnableImplicitMT()

INPUT_BASE = "/eos/user/l/liwe/monoZ_combine"
OUTPUT_FILE = "shapes_ptmiss_statonly.root"

BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
SIGNAL_FILE = "ZH_ZToLL_HToInvisible_M125_UL17_local.root"

# output_name: name used in datacard/workspace
# source_dir:  subdirectory in monoZ_combine
CHANNELS = [
    {
        "output_name": "SR",
        "source_dir": "SR",
        "selection": "lepton_cat!=2 && ptmiss>=100",
        # Robust SR binning: merge very-low-yield high-ptmiss tail bins.
        "bins": [100, 140, 200, 320, 550, 3000],
        "has_signal": True,
    },
    {
        "output_name": "DYCR",
        "source_dir": "DY_CR",
        "selection": "lepton_cat!=2 && ptmiss>=30 && ptmiss<=90",
        # Drop the 90-95 empty tail bin since DY_CR is cut at ptmiss<=90.
        "bins": [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90],
        "has_signal": False,
    },
    {
        "output_name": "EMUCR",
        "source_dir": "emu_CR",
        "selection": "lepton_cat==2",
        # Coarse EMUCR binning to avoid unstable sparse tail bins.
        "bins": [100, 140, 180, 1100],
        "has_signal": False,
    },
    {
        "output_name": "CR3L",
        "source_dir": "3l_CR",
        "selection": "1",
        # Coarse 3l control-region binning with merged high-ptmiss tail.
        "bins": [30, 70, 120, 220, 420, 2000],
        "has_signal": False,
    },
]


def build_hist(file_path, selection, bins, hist_name, weight_expr=None):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    rdf = ROOT.RDataFrame("Vars", file_path)
    if selection and selection != "1":
        rdf = rdf.Filter(selection)
    rdf = rdf.Filter("std::isfinite(ptmiss)")

    model = ROOT.RDF.TH1DModel(hist_name, "", len(bins) - 1, array("d", bins))
    if weight_expr:
        rdf = rdf.Filter(f"std::isfinite({weight_expr})")
        rdf = rdf.Define("w_for_hist", weight_expr)
        h = rdf.Histo1D(model, "ptmiss", "w_for_hist").GetValue().Clone(hist_name)
    else:
        h = rdf.Histo1D(model, "ptmiss").GetValue().Clone(hist_name)

    h.SetDirectory(ROOT.nullptr)
    return h


def floor_mc_hist(hist, context):
    for i_bin in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(i_bin) <= 0:
            print(
                f"  [WARN] {context}: bin {i_bin} has non-positive yield "
                f"{hist.GetBinContent(i_bin):.6g}; flooring to 1e-6"
            )
            hist.SetBinContent(i_bin, 1e-6)
        if hist.GetBinError(i_bin) <= 0:
            hist.SetBinError(i_bin, 1e-6)


def integral(hist):
    return float(hist.Integral(1, hist.GetNbinsX()))


def main():
    out = ROOT.TFile(OUTPUT_FILE, "RECREATE")

    print(f"Writing {OUTPUT_FILE}")
    for ch in CHANNELS:
        out_name = ch["output_name"]
        src_dir = ch["source_dir"]
        selection = ch["selection"]
        bins = ch["bins"]

        print(f"\n[{out_name}] source={src_dir} selection=({selection})")
        out.mkdir(out_name)
        out.cd(out_name)

        # Data
        data_path = os.path.join(INPUT_BASE, src_dir, "Data.root")
        h_data = build_hist(data_path, selection, bins, "data_obs", None)
        h_data.Write("data_obs")
        data_sum = h_data.Integral(1, h_data.GetNbinsX())
        print(f"  data_obs integral = {data_sum:.3f}")

        # Backgrounds
        for proc in BKG_PROCESSES:
            p_path = os.path.join(INPUT_BASE, src_dir, f"{proc}.root")
            h = build_hist(p_path, selection, bins, proc, "weight")
            floor_mc_hist(h, f"{out_name}/{proc}")
            h.Write(proc)
            yld = integral(h)
            print(f"  {proc:6s}: {yld:10.3f}")

        # Signal only in SR
        if ch["has_signal"]:
            s_path = os.path.join(INPUT_BASE, src_dir, SIGNAL_FILE)
            h_sig = build_hist(s_path, selection, bins, "signal", "weight")
            floor_mc_hist(h_sig, f"{out_name}/signal")
            h_sig.Write("signal")
            yld = integral(h_sig)
            print(f"  {'signal':6s}: {yld:10.3f}")

    out.Close()
    print("\nDone.")


if __name__ == "__main__":
    main()
