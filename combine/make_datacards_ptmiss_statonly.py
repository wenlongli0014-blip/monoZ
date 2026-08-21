#!/usr/bin/env python3

import os
import subprocess

SHAPES_ABS = os.path.abspath("shapes_ptmiss_statonly.root")
CARDS_DIR = "cards"

# Keep the same process names as in shapes file.
BKG = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]


def write_channel_card(channel, has_signal):
    if has_signal:
        processes = ["signal"] + BKG
        proc_ids = [0] + list(range(1, len(BKG) + 1))
    else:
        processes = BKG[:]
        proc_ids = list(range(1, len(BKG) + 1))

    lines = []
    lines.append("imax 1")
    lines.append("jmax *")
    lines.append("kmax *")
    lines.append("------------")
    lines.append(f"shapes * {channel} {SHAPES_ABS} {channel}/$PROCESS")
    lines.append("------------")
    lines.append(f"bin {channel}")
    lines.append("observation -1")
    lines.append("------------")
    lines.append("bin " + " ".join([channel] * len(processes)))
    lines.append("process " + " ".join(processes))
    lines.append("process " + " ".join(str(x) for x in proc_ids))
    lines.append("rate " + " ".join(["-1"] * len(processes)))
    lines.append("------------")
    lines.append("* autoMCStats 0")
    lines.append("")

    out_path = os.path.join(CARDS_DIR, f"{channel}.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def normalize_combined_card_header(card_path):
    with open(card_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("kmax "):
            lines[i] = "kmax * number of nuisance parameters\n"
            break

    with open(card_path, "w", encoding="utf-8") as f:
        f.writelines(lines)


def main():
    os.makedirs(CARDS_DIR, exist_ok=True)

    write_channel_card("SR", has_signal=True)
    write_channel_card("DYCR", has_signal=False)
    write_channel_card("EMUCR", has_signal=False)
    write_channel_card("CR3L", has_signal=False)

    combined = os.path.join(CARDS_DIR, "combined_statonly.txt")
    cmd = (
        f"combineCards.py SR={CARDS_DIR}/SR.txt DYCR={CARDS_DIR}/DYCR.txt "
        f"EMUCR={CARDS_DIR}/EMUCR.txt CR3L={CARDS_DIR}/CR3L.txt > {combined}"
    )
    subprocess.run(cmd, shell=True, check=True)
    normalize_combined_card_header(combined)

    with open(combined, "a", encoding="utf-8") as f:
        f.write("\n")
        f.write("# Shared normalization parameters (user-requested minimal setup)\n")
        f.write("# DY_CR controls DY in SR (3l_CR does NOT share k_DY by request).\n")
        f.write("k_DY    rateParam SR    DY    1 [0,5]\n")
        f.write("k_DY    rateParam DYCR  DY    1 [0,5]\n")
        f.write("\n")
        f.write("# 3l_CR controls WZ in SR.\n")
        f.write("k_WZ    rateParam SR    WZ    1 [0,5]\n")
        f.write("k_WZ    rateParam CR3L  WZ    1 [0,5]\n")
        f.write("\n")
        f.write("# emu_CR controls only ttbar in SR.\n")
        f.write("k_TT    rateParam SR    ttbar 1 [0,5]\n")
        f.write("k_TT    rateParam EMUCR ttbar 1 [0,5]\n")

    print(f"Wrote {combined}")


if __name__ == "__main__":
    main()
