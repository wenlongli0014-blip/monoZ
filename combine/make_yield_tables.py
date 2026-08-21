#!/usr/bin/env python3

import argparse
import csv
import os

import numpy as np
import uproot


PROCESSES = ["signal", "DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
SECTIONS = {
    "Before fit": "shapes_prefit",
    "After fit s+b": "shapes_fit_s",
    "After fit b-only": "shapes_fit_b",
}


def hist_sum(directory, name):
    obj = directory.get(name)
    if obj is None:
        return 0.0
    return float(np.sum(obj.values()))


def data_sum(directory):
    graph = directory["data"]
    _, y_values = graph.values()
    return float(np.sum(y_values))


def build_rows(root_file, channels):
    rows = []
    for section, fit_dir in SECTIONS.items():
        for channel in channels:
            directory = root_file[f"{fit_dir}/{channel}"]
            row = {"section": section, "channel": channel, "Data": data_sum(directory)}
            total_bkg = 0.0
            for process in PROCESSES:
                value = hist_sum(directory, process)
                row[process] = value
                if process != "signal":
                    total_bkg += value
            row["Total background"] = total_bkg
            row["Total signal+background"] = total_bkg + row["signal"]
            rows.append(row)
    return rows


def write_csv(rows, path):
    fieldnames = [
        "section",
        "channel",
        "Data",
        "signal",
        "DY",
        "Other",
        "ST",
        "VVV",
        "WW",
        "WZ",
        "ggZZ",
        "qqZZ",
        "ttbar",
        "Total background",
        "Total signal+background",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{row[key]:.6f}" if isinstance(row[key], float) else row[key] for key in fieldnames})


def write_markdown(rows, path):
    columns = [
        "channel",
        "Data",
        "signal",
        "DY",
        "WZ",
        "ttbar",
        "qqZZ",
        "ggZZ",
        "WW",
        "ST",
        "VVV",
        "Other",
        "Total background",
        "Total signal+background",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        for section in SECTIONS:
            handle.write(f"## {section}\n\n")
            handle.write("| " + " | ".join(columns) + " |\n")
            handle.write("| " + " | ".join(["---"] * len(columns)) + " |\n")
            for row in rows:
                if row["section"] != section:
                    continue
                values = []
                for column in columns:
                    value = row[column]
                    values.append(f"{value:.3f}" if isinstance(value, float) else value)
                handle.write("| " + " | ".join(values) + " |\n")
            handle.write("\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Create yield tables from fitDiagnostics output.")
    parser.add_argument("--fitdiagnostics", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--outdir", default="tables")
    parser.add_argument("--channels", nargs="+", default=["SR", "DYCR", "EMUCR", "CR3L"])
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with uproot.open(args.fitdiagnostics) as root_file:
        rows = build_rows(root_file, args.channels)

    csv_path = os.path.join(args.outdir, f"yields_{args.tag}.csv")
    md_path = os.path.join(args.outdir, f"yields_{args.tag}.md")
    write_csv(rows, csv_path)
    write_markdown(rows, md_path)
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
