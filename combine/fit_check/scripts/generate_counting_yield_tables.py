#!/usr/bin/env python3

import csv
import os
from collections import defaultdict

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

PROCESS_ORDER = [
    "signal",
    "DY",
    "WZ",
    "ttbar",
    "WW",
    "ST",
    "qqZZ",
    "ggZZ",
    "VVV",
    "Other",
]


def read_csv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def read_manifest(path):
    with open(path) as handle:
        return {row["tag"]: row for row in csv.DictReader(handle, delimiter="\t")}


def format_value(value):
    return f"{value:.6f}"


def should_keep_signal(channel, prefit, postfit):
    return channel == "SR" or abs(prefit) > 1e-6 or abs(postfit) > 1e-6


def main():
    workdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tables_dir = os.path.join(workdir, "tables")

    manifest = read_manifest(os.path.join(tables_dir, "manifest.tsv"))
    fit_rows = read_csv(os.path.join(tables_dir, "fit_parameters.csv"))
    yield_rows = read_csv(os.path.join(tables_dir, "yield_summary.csv"))

    fit_params = defaultdict(list)
    for row in fit_rows:
        fit_params[row["tag"]].append(row)

    channel_summary = defaultdict(list)
    for row in yield_rows:
        if row["tag"].startswith("counting_"):
            channel_summary[row["tag"]].append(row)

    summary_rows = []
    summary_md = []
    summary_md.append("# Counting Fit Yield Tables\n")
    summary_md.append(
        "Inputs: `/eos/user/l/liwe/monoZ_combine`; `EMUCR` uses `lepton_cat==2 && ptmiss>=100`.\n"
    )
    summary_md.append(
        "These tables correspond to the six one-bin counting fits in `fit_check`. "
        "Each channel is shown once, with `Observed data` written in the channel header rather than repeated on every row.\n"
    )

    for tag in sorted(channel_summary.keys()):
        info = manifest[tag]
        shapes_path = os.path.join(workdir, info["shapes"])
        fit_path = os.path.join(workdir, info["fit_root"])

        shapes_file = ROOT.TFile.Open(shapes_path)
        fit_file = ROOT.TFile.Open(fit_path)
        if not shapes_file or shapes_file.IsZombie():
            raise OSError(f"Could not open shapes file: {shapes_path}")
        if not fit_file or fit_file.IsZombie():
            raise OSError(f"Could not open fit file: {fit_path}")

        lines = [f"# {tag}\n", f"Description: {info['description']}\n"]
        if fit_params.get(tag):
            rendered = "; ".join(
                f"{row['param']} = {float(row['value']):.6g} +/- {float(row['error']):.6g}"
                for row in fit_params[tag]
            )
            lines.append(f"Fit parameters: {rendered}\n")

        summary_md.append(f"## {tag}\n")
        summary_md.append(f"Description: {info['description']}\n")
        if fit_params.get(tag):
            rendered = "; ".join(
                f"{row['param']} = {float(row['value']):.6g} +/- {float(row['error']):.6g}"
                for row in fit_params[tag]
            )
            summary_md.append(f"Fit parameters: {rendered}\n")

        for channel_row in channel_summary[tag]:
            channel = channel_row["channel"]
            data_obs = float(channel_row["data_obs"])
            prefit_bkg = float(channel_row["prefit_background"])
            postfit_bkg = float(channel_row["postfit_background"] or 0.0)

            lines.append(f"## {channel}\n")
            lines.append(f"Observed data: {format_value(data_obs)}\n")
            lines.append("| Process | Prefit | Postfit |")
            lines.append("|---|---:|---:|")

            summary_md.append(f"### {channel}\n")
            summary_md.append(f"Observed data: {format_value(data_obs)}\n")
            summary_md.append("| Process | Prefit | Postfit |")
            summary_md.append("|---|---:|---:|")

            shapes_dir = shapes_file.Get(channel)
            postfit_dir = fit_file.Get(f"shapes_fit_s/{channel}")
            if not shapes_dir:
                raise KeyError(f"Missing channel {channel} in {shapes_path}")
            if not postfit_dir:
                raise KeyError(f"Missing channel {channel} in {fit_path}")

            for process in PROCESS_ORDER:
                prefit_hist = shapes_dir.Get(process)
                postfit_hist = postfit_dir.Get(process)
                prefit = float(prefit_hist.Integral()) if prefit_hist else 0.0
                postfit = float(postfit_hist.Integral()) if postfit_hist else 0.0

                if process == "signal" and not should_keep_signal(channel, prefit, postfit):
                    continue
                if process != "signal" and abs(prefit) < 1e-12 and abs(postfit) < 1e-12:
                    continue

                row = f"| {process} | {format_value(prefit)} | {format_value(postfit)} |"
                lines.append(row)
                summary_md.append(row)
                summary_rows.append(
                    {
                        "tag": tag,
                        "channel": channel,
                        "process": process,
                        "data_obs": data_obs,
                        "prefit_yield": prefit,
                        "postfit_yield": postfit,
                    }
                )

            total_row = (
                f"| Total background | {format_value(prefit_bkg)} | {format_value(postfit_bkg)} |"
            )
            lines.append(total_row)
            summary_md.append(total_row)
            lines.append("")
            summary_md.append("")

            summary_rows.append(
                {
                    "tag": tag,
                    "channel": channel,
                    "process": "Total background",
                    "data_obs": data_obs,
                    "prefit_yield": prefit_bkg,
                    "postfit_yield": postfit_bkg,
                }
            )

        shapes_file.Close()
        fit_file.Close()

        with open(os.path.join(tables_dir, f"{tag}_table.md"), "w") as handle:
            handle.write("\n".join(lines).rstrip() + "\n")

    with open(os.path.join(tables_dir, "counting_yields_before_after.md"), "w") as handle:
        handle.write("\n".join(summary_md).rstrip() + "\n")

    with open(os.path.join(tables_dir, "counting_yields_before_after.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tag", "channel", "process", "data_obs", "prefit_yield", "postfit_yield"],
        )
        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
