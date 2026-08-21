#!/usr/bin/env python3

"""Scan predefined Run-2 SR ptmiss binnings with the previous-limit model."""

import argparse
import csv
import json
import os
import re
import subprocess
import sys

import ROOT


REPO_COMBINE = os.path.dirname(os.path.abspath(__file__))
SHAPE_SCRIPT = os.path.join(REPO_COMBINE, "make_shapes_ptmiss_2017_2018_blinded.py")
CARD_SCRIPT = os.path.join(REPO_COMBINE, "fit_check", "scripts", "make_fit_check_datacard.py")

K_FACTORS = {
    "k_Zjet": 1.27457654256,
    "k_WZ": 0.992271991508,
    "k_emu": 1.24039762208,
}

SCALED_PROCESSES = {
    "DY": K_FACTORS["k_Zjet"],
    "WZ": K_FACTORS["k_WZ"],
    "ttbar": K_FACTORS["k_emu"],
    "WW": K_FACTORS["k_emu"],
    "ST": K_FACTORS["k_emu"],
}

BACKGROUNDS = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]

CANDIDATES = {
    "baseline_5": [100, 140, 200, 320, 550, 3000],
    "coarse_4": [100, 160, 240, 400, 3000],
    "balanced_8": [100, 125, 155, 195, 250, 330, 450, 650, 3000],
    "low_fine_9": [100, 120, 140, 170, 200, 250, 320, 420, 550, 3000],
    "high_split_8": [100, 140, 200, 280, 400, 550, 750, 1000, 3000],
    "uniformish_11": [100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 3000],
    "variable_14": [100, 120, 140, 160, 180, 200, 230, 270, 320, 380, 460, 550, 700, 1000, 3000],
    "fine_20": [100, 110, 120, 130, 140, 150, 160, 175, 190, 210, 235, 260, 290, 325, 365, 410, 470, 550, 700, 1000, 3000],
    "fine_low_18": [100, 110, 120, 130, 140, 150, 160, 175, 190, 210, 235, 260, 290, 325, 365, 410, 470, 550, 3000],
    "fine_low_tail2_19": [100, 110, 120, 130, 140, 150, 160, 175, 190, 210, 235, 260, 290, 325, 365, 410, 470, 550, 700, 3000],
    "fine_low5_26": [100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 170, 180, 190, 200, 215, 235, 260, 290, 325, 365, 410, 470, 550, 3000],
}


def run(command, cwd, log_path):
    with open(log_path, "w", encoding="utf-8") as log:
        subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=True)


def parse_limits(log_path):
    patterns = {
        "observed_asimov": r"Observed Limit: r < ([0-9.eE+-]+)",
        "expected_2p5": r"Expected\s+2\.5%: r < ([0-9.eE+-]+)",
        "expected_16": r"Expected\s+16\.0%: r < ([0-9.eE+-]+)",
        "expected_50": r"Expected\s+50\.0%: r < ([0-9.eE+-]+)",
        "expected_84": r"Expected\s+84\.0%: r < ([0-9.eE+-]+)",
        "expected_97p5": r"Expected\s+97\.5%: r < ([0-9.eE+-]+)",
    }
    with open(log_path, encoding="utf-8") as handle:
        text = handle.read()
    values = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        if match is None:
            raise RuntimeError(f"Could not find {name} in {log_path}")
        values[name] = float(match.group(1))
    return values


def bin_quality(shape_path):
    root_file = ROOT.TFile.Open(shape_path)
    if not root_file or root_file.IsZombie():
        raise OSError(shape_path)
    signal = root_file.Get("SR/signal")
    nbins = signal.GetNbinsX()
    rows = []
    for index in range(1, nbins + 1):
        background = 0.0
        variance = 0.0
        for process in BACKGROUNDS:
            hist = root_file.Get(f"SR/{process}")
            scale = SCALED_PROCESSES.get(process, 1.0)
            background += scale * hist.GetBinContent(index)
            variance += (scale * hist.GetBinError(index)) ** 2
        sig = signal.GetBinContent(index)
        rows.append(
            {
                "low": signal.GetXaxis().GetBinLowEdge(index),
                "high": signal.GetXaxis().GetBinUpEdge(index),
                "signal": sig,
                "background": background,
                "s_over_b": sig / background if background > 0 else 0.0,
                "mc_neff": background * background / variance if variance > 0 else 0.0,
            }
        )
    root_file.Close()
    return rows


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-base", default="/eos/user/l/liwe/monoz_full_dnn")
    parser.add_argument("--out", required=True)
    parser.add_argument("--only", nargs="*", choices=sorted(CANDIDATES))
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    selected = args.only or list(CANDIDATES)
    summaries = []

    for tag in selected:
        edges = CANDIDATES[tag]
        candidate_dir = os.path.join(args.out, tag)
        logs_dir = os.path.join(candidate_dir, "logs")
        cards_dir = os.path.join(candidate_dir, "cards")
        workspaces_dir = os.path.join(candidate_dir, "workspaces")
        os.makedirs(logs_dir, exist_ok=True)
        os.makedirs(cards_dir, exist_ok=True)
        os.makedirs(workspaces_dir, exist_ok=True)
        shape_path = os.path.join(candidate_dir, "shapes.root")
        card_path = os.path.join(cards_dir, "final.txt")
        workspace_path = os.path.join(workspaces_dir, "final.root")

        print(f"[{tag}] edges = {edges}", flush=True)
        run(
            [
                sys.executable,
                SHAPE_SCRIPT,
                "--input-base",
                args.input_base,
                "--years",
                "2016",
                "2017",
                "2018",
                "--sr-bins",
                *[str(edge) for edge in edges],
                "--out",
                shape_path,
            ],
            args.out,
            os.path.join(logs_dir, "01_shapes.log"),
        )
        run(
            [
                sys.executable,
                CARD_SCRIPT,
                "--shapes",
                shape_path,
                "--tag",
                "final",
                "--channels",
                "SR",
                "DYCR",
                "CR3L",
                "EMUCR",
                "--k-factors",
                "k_Zjet",
                "k_WZ",
                "k_emu",
                "--out-dir",
                cards_dir,
            ],
            args.out,
            os.path.join(logs_dir, "02_card.log"),
        )
        run(
            ["text2workspace.py", card_path, "-o", workspace_path],
            candidate_dir,
            os.path.join(logs_dir, "03_workspace.log"),
        )
        limit_log = os.path.join(logs_dir, "04_limit.log")
        run(
            [
                "combine",
                "-M",
                "AsymptoticLimits",
                workspace_path,
                "-m",
                "125",
                "-t",
                "-1",
                "--expectSignal",
                "0",
                "--setParameters",
                "k_Zjet=1.27457654256,k_WZ=0.992271991508,k_emu=1.24039762208",
                "--redefineSignalPOIs",
                "r",
                "--rMin",
                "0",
                "--rMax",
                "5",
                "--cminDefaultMinimizerStrategy",
                "0",
                "-n",
                f".{tag}",
            ],
            candidate_dir,
            limit_log,
        )
        limits = parse_limits(limit_log)
        quality = bin_quality(shape_path)
        with open(os.path.join(candidate_dir, "bin_quality.json"), "w", encoding="utf-8") as handle:
            json.dump(quality, handle, indent=2)
            handle.write("\n")
        summary = {
            "tag": tag,
            "n_bins": len(edges) - 1,
            "edges": " ".join(str(edge) for edge in edges),
            **limits,
            "min_background": min(row["background"] for row in quality),
            "min_mc_neff": min(row["mc_neff"] for row in quality),
            "max_s_over_b": max(row["s_over_b"] for row in quality),
        }
        summaries.append(summary)
        print(
            f"[{tag}] median={summary['expected_50']:.6g}, "
            f"minB={summary['min_background']:.3g}, "
            f"minNeff={summary['min_mc_neff']:.3g}",
            flush=True,
        )

    summaries.sort(key=lambda row: row["expected_50"])
    fields = list(summaries[0])
    csv_path = os.path.join(args.out, "scan_summary.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries)
    with open(os.path.join(args.out, "scan_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summaries, handle, indent=2)
        handle.write("\n")
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
