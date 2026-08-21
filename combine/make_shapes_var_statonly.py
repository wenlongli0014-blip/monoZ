#!/usr/bin/env python3

import argparse
import os
from array import array

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.ROOT.EnableImplicitMT()

INPUT_BASE = "/eos/user/l/liwe/monoZ_combine"

BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
SIGNAL_FILE = "ZH_ZToLL_HToInvisible_M125_UL17_local.root"

SUPPORTED_VARIABLES = ["ptmiss", "mT", "ptmiss_significance_corrected"]

# output_name: name used in datacard/workspace
# source_dir:  subdirectory in monoZ_combine
CHANNELS = [
    {
        "output_name": "SR",
        "source_dir": "SR",
        "selection": "lepton_cat!=2 && ptmiss>=100",
        "has_signal": True,
        "bins": {
            "ptmiss": [100, 140, 200, 320, 550, 3000],
            "mT": [240, 280, 320, 380, 460, 600, 900, 5000],
            "ptmiss_significance_corrected": [12, 16, 22, 32, 50, 80, 130, 220, 500, 12000],
        },
    },
    {
        "output_name": "DYCR",
        "source_dir": "DY_CR",
        "selection": "lepton_cat!=2 && ptmiss>=30 && ptmiss<=90",
        "has_signal": False,
        "bins": {
            "ptmiss": [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90],
            "mT": [210, 240, 260, 280, 320, 380, 500, 900],
            "ptmiss_significance_corrected": [12, 12.5, 13, 14, 16, 18, 22, 28, 40, 200],
        },
    },
    {
        "output_name": "EMUCR",
        "source_dir": "emu_CR",
        "selection": "lepton_cat==2",
        "has_signal": False,
        "bins": {
            "ptmiss": [100, 140, 180, 1100],
            "mT": [240, 280, 320, 380, 500, 1600],
            "ptmiss_significance_corrected": [12, 15, 20, 30, 50, 100, 3000],
        },
    },
    {
        "output_name": "CR3L",
        "source_dir": "3l_CR",
        "selection": "1",
        "has_signal": False,
        "bins": {
            "ptmiss": [30, 70, 120, 220, 420, 2000],
            "mT": [160, 190, 220, 260, 320, 420, 700, 3500],
            "ptmiss_significance_corrected": [3, 4, 5.5, 8, 12, 18, 30, 60, 200, 8000],
        },
    },
]


def sanitize_tag(text):
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def build_hist(file_path, selection, variable, bins, hist_name, weight_expr=None):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    rdf = ROOT.RDataFrame("Vars", file_path)
    if selection and selection != "1":
        rdf = rdf.Filter(selection)
    rdf = rdf.Filter(f"std::isfinite({variable})")

    model = ROOT.RDF.TH1DModel(hist_name, "", len(bins) - 1, array("d", bins))
    if weight_expr:
        rdf = rdf.Filter(f"std::isfinite({weight_expr})")
        rdf = rdf.Define("w_for_hist", weight_expr)
        h = rdf.Histo1D(model, variable, "w_for_hist").GetValue().Clone(hist_name)
    else:
        h = rdf.Histo1D(model, variable).GetValue().Clone(hist_name)

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


def parse_args():
    parser = argparse.ArgumentParser(description="Build stat-only shape templates for a chosen variable.")
    parser.add_argument("--var", required=True, choices=SUPPORTED_VARIABLES, help="Shape variable name")
    parser.add_argument(
        "--out",
        default=None,
        help="Output ROOT file path (default: shapes_<var>_statonly.root)",
    )
    parser.add_argument(
        "--min-value",
        type=float,
        default=None,
        help="Optional lower edge override for the chosen variable. Events below this value go to underflow.",
    )
    parser.add_argument(
        "--bin-edges",
        default=None,
        help="Optional comma-separated bin edges to use for all channels, e.g. '16,17,18,36,12000'.",
    )
    parser.add_argument(
        "--sr-bin-edges",
        default=None,
        help="Optional comma-separated bin edges to use only for SR. Other channels keep their default bins.",
    )
    parser.add_argument(
        "--input-base",
        default=INPUT_BASE,
        help=f"Input directory containing SR/DY_CR/emu_CR/3l_CR subdirectories (default: {INPUT_BASE})",
    )
    parser.add_argument(
        "--signal-file",
        default=SIGNAL_FILE,
        help=f"Signal ROOT file name inside the SR directory (default: {SIGNAL_FILE})",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    var = args.var
    default_out = f"shapes_{sanitize_tag(var)}_statonly.root"
    output_file = args.out if args.out else default_out
    custom_bins = None
    if args.bin_edges:
        custom_bins = [float(edge) for edge in args.bin_edges.split(",") if edge.strip()]
        if len(custom_bins) < 2:
            raise ValueError("--bin-edges must contain at least two edges")
        if any(high <= low for low, high in zip(custom_bins[:-1], custom_bins[1:])):
            raise ValueError(f"--bin-edges must be strictly increasing: {custom_bins}")
    sr_custom_bins = None
    if args.sr_bin_edges:
        sr_custom_bins = [float(edge) for edge in args.sr_bin_edges.split(",") if edge.strip()]
        if len(sr_custom_bins) < 2:
            raise ValueError("--sr-bin-edges must contain at least two edges")
        if any(high <= low for low, high in zip(sr_custom_bins[:-1], sr_custom_bins[1:])):
            raise ValueError(f"--sr-bin-edges must be strictly increasing: {sr_custom_bins}")

    out = ROOT.TFile(output_file, "RECREATE")

    print(f"Writing {output_file}")
    print(f"Variable = {var}")
    print(f"Input base = {args.input_base}")
    print(f"Signal file = {args.signal_file}")
    for ch in CHANNELS:
        out_name = ch["output_name"]
        src_dir = ch["source_dir"]
        selection = ch["selection"]
        if sr_custom_bins is not None and out_name == "SR":
            bins = sr_custom_bins
        else:
            bins = custom_bins if custom_bins is not None else ch["bins"][var]
        if args.min_value is not None:
            if args.min_value >= bins[-1]:
                raise ValueError(f"--min-value {args.min_value} is above the last bin edge for {out_name}: {bins[-1]}")
            bins = [args.min_value] + [edge for edge in bins if edge > args.min_value]

        print(f"\n[{out_name}] source={src_dir} selection=({selection})")
        out.mkdir(out_name)
        out.cd(out_name)

        # Data
        data_path = os.path.join(args.input_base, src_dir, "Data.root")
        h_data = build_hist(data_path, selection, var, bins, "data_obs", None)
        h_data.Write("data_obs")
        data_sum = h_data.Integral(1, h_data.GetNbinsX())
        print(f"  data_obs integral = {data_sum:.3f}")

        # Backgrounds
        for proc in BKG_PROCESSES:
            p_path = os.path.join(args.input_base, src_dir, f"{proc}.root")
            h = build_hist(p_path, selection, var, bins, proc, "weight")
            floor_mc_hist(h, f"{out_name}/{proc}")
            h.Write(proc)
            yld = integral(h)
            print(f"  {proc:6s}: {yld:10.3f}")

        # Signal only in SR
        if ch["has_signal"]:
            s_path = os.path.join(args.input_base, src_dir, args.signal_file)
            h_sig = build_hist(s_path, selection, var, bins, "signal", "weight")
            floor_mc_hist(h_sig, f"{out_name}/signal")
            h_sig.Write("signal")
            yld = integral(h_sig)
            print(f"  {'signal':6s}: {yld:10.3f}")

    out.Close()
    print("\nDone.")


if __name__ == "__main__":
    main()
