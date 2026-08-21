#!/usr/bin/env python3

import argparse
import csv
import os

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

PARAMETERS = ["r", "k_Zjet", "k_WZ", "k_emu"]
BACKGROUNDS = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]


def get_fit_result(root_file):
    for name in ["fit_s", "fit_b"]:
        result = root_file.Get(name)
        if result:
            return name, result
    return "", None


def append_row(path, fieldnames, row):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    write_header = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def integral(hist):
    if not hist:
        return ""
    return f"{hist.Integral(1, hist.GetNbinsX()):.10g}"


def sum_processes(directory, processes):
    total = 0.0
    found = False
    for process in processes:
        hist = directory.Get(process) if directory else None
        if hist:
            total += hist.Integral(1, hist.GetNbinsX())
            found = True
    return f"{total:.10g}" if found else ""


def parse_args():
    parser = argparse.ArgumentParser(description="Append parameter, status, and yield summaries for one fit.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--fit", required=True)
    parser.add_argument("--shapes", required=True)
    parser.add_argument("--channels", nargs="+", required=True)
    parser.add_argument("--params-out", default="tables/fit_parameters.csv")
    parser.add_argument("--status-out", default="tables/fit_status.csv")
    parser.add_argument("--yields-out", default="tables/yield_summary.csv")
    return parser.parse_args()


def main():
    args = parse_args()
    fit_file = ROOT.TFile.Open(args.fit)
    if not fit_file or fit_file.IsZombie():
        raise OSError(f"Could not open {args.fit}")
    shapes_file = ROOT.TFile.Open(args.shapes)
    if not shapes_file or shapes_file.IsZombie():
        raise OSError(f"Could not open {args.shapes}")

    fit_name, fit_result = get_fit_result(fit_file)
    status_row = {
        "tag": args.tag,
        "fit_file": args.fit,
        "fit_result": fit_name,
        "status": "",
        "covQual": "",
        "edm": "",
        "minNll": "",
    }
    if fit_result:
        status_row.update(
            {
                "status": fit_result.status(),
                "covQual": fit_result.covQual(),
                "edm": f"{fit_result.edm():.10g}",
                "minNll": f"{fit_result.minNll():.10g}",
            }
        )
        params = fit_result.floatParsFinal()
        for param_name in PARAMETERS:
            param = params.find(param_name)
            if not param:
                continue
            append_row(
                args.params_out,
                ["tag", "param", "value", "error", "error_lo", "error_hi"],
                {
                    "tag": args.tag,
                    "param": param_name,
                    "value": f"{param.getVal():.10g}",
                    "error": f"{param.getError():.10g}",
                    "error_lo": f"{param.getAsymErrorLo():.10g}",
                    "error_hi": f"{param.getAsymErrorHi():.10g}",
                },
            )

    append_row(
        args.status_out,
        ["tag", "fit_file", "fit_result", "status", "covQual", "edm", "minNll"],
        status_row,
    )

    for channel in args.channels:
        input_dir = shapes_file.GetDirectory(channel)
        fit_dir = fit_file.GetDirectory(f"shapes_{fit_name}/{channel}") if fit_name else None
        postfit_bkg = integral(fit_dir.Get("total_background")) if fit_dir else ""
        postfit_sig = integral(fit_dir.Get("total_signal")) if fit_dir else ""
        append_row(
            args.yields_out,
            [
                "tag",
                "channel",
                "data_obs",
                "prefit_background",
                "prefit_signal",
                "postfit_background",
                "postfit_signal",
            ],
            {
                "tag": args.tag,
                "channel": channel,
                "data_obs": integral(input_dir.Get("data_obs")) if input_dir else "",
                "prefit_background": sum_processes(input_dir, BACKGROUNDS),
                "prefit_signal": integral(input_dir.Get("signal")) if input_dir else "",
                "postfit_background": postfit_bkg,
                "postfit_signal": postfit_sig,
            },
        )

    fit_file.Close()
    shapes_file.Close()


if __name__ == "__main__":
    main()
