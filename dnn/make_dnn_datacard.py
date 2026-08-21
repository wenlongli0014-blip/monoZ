#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess


BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]


CHANNELS = ["SR", "DYCR", "EMUCR", "CR3L"]


def write_channel_card(channel, has_signal, shapes_abs, cards_dir):
    processes = ["signal"] + BKG_PROCESSES if has_signal else BKG_PROCESSES[:]
    proc_ids = [0] + list(range(1, len(BKG_PROCESSES) + 1)) if has_signal else list(range(1, len(BKG_PROCESSES) + 1))
    lines = [
        "imax 1",
        "jmax *",
        "kmax *",
        "------------",
        f"shapes * {channel} {shapes_abs} {channel}/$PROCESS",
        "------------",
        f"bin {channel}",
        "observation -1",
        "------------",
        "bin " + " ".join([channel] * len(processes)),
        "process " + " ".join(processes),
        "process " + " ".join(str(pid) for pid in proc_ids),
        "rate " + " ".join(["-1"] * len(processes)),
        "------------",
        "",
    ]
    with open(os.path.join(cards_dir, f"{channel}.txt"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def normalize_combined_card_header(card_path):
    with open(card_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    for idx, line in enumerate(lines):
        if line.startswith("kmax "):
            lines[idx] = "kmax * number of nuisance parameters\n"
            break
    with open(card_path, "w", encoding="utf-8") as handle:
        handle.writelines(lines)


def parse_args():
    parser = argparse.ArgumentParser(description="Build a DNN score datacard without shape/systematic uncertainties.")
    parser.add_argument("--shapes", required=True)
    parser.add_argument("--tag", default="dnn_score")
    parser.add_argument("--outdir", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/cards")
    return parser.parse_args()


def main():
    args = parse_args()
    tag = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in args.tag)
    shapes_abs = os.path.abspath(args.shapes)
    cards_dir = os.path.join(args.outdir, tag)
    os.makedirs(cards_dir, exist_ok=True)

    write_channel_card("SR", True, shapes_abs, cards_dir)
    write_channel_card("DYCR", False, shapes_abs, cards_dir)
    write_channel_card("EMUCR", False, shapes_abs, cards_dir)
    write_channel_card("CR3L", False, shapes_abs, cards_dir)

    combined = os.path.join(args.outdir, f"combined_{tag}.txt")
    cmd = (
        f"combineCards.py SR={cards_dir}/SR.txt DYCR={cards_dir}/DYCR.txt "
        f"EMUCR={cards_dir}/EMUCR.txt CR3L={cards_dir}/CR3L.txt > {combined}"
    )
    subprocess.run(cmd, shell=True, check=True)
    normalize_combined_card_header(combined)

    with open(combined, "a", encoding="utf-8") as handle:
        handle.write("\n")
        handle.write("# Shared control-region normalization parameters.\n")
        handle.write("k_Zjet  rateParam *     DY    1 [0,5]\n")
        handle.write("\n")
        handle.write("k_WZ    rateParam *     WZ    1 [0,5]\n")
        handle.write("\n")
        handle.write("k_emu   rateParam *     ttbar 1 [0,5]\n")
        handle.write("k_emu   rateParam *     WW    1 [0,5]\n")
        handle.write("k_emu   rateParam *     ST    1 [0,5]\n")

    print(f"Wrote {combined}")


if __name__ == "__main__":
    main()
