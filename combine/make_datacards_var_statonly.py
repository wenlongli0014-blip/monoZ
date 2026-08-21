#!/usr/bin/env python3

import argparse
import os
import subprocess

# Keep the same process names as in shapes file.
BKG = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
CHANNELS = ["SR", "DYCR", "EMUCR", "CR3L"]

# Flat normalization uncertainties requested as a temporary dummy model to
# absorb broad data/MC tension. Values are written as lnN nuisances.
FLAT_DUMMY_UNCERTAINTIES = [
    ("dummy_DY_norm", {"DY": "1.40"}),
    ("dummy_WZ_norm", {"WZ": "1.10"}),
    ("dummy_ZZ_norm", {"ggZZ": "1.15", "qqZZ": "1.15"}),
    ("dummy_ttbar_norm", {"ttbar": "1.15"}),
    ("dummy_WW_norm", {"WW": "1.15"}),
    ("dummy_ST_norm", {"ST": "1.15"}),
]


def sanitize_tag(text):
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def write_channel_card(channel, has_signal, shapes_abs, cards_dir):
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
    lines.append(f"shapes * {channel} {shapes_abs} {channel}/$PROCESS")
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

    out_path = os.path.join(cards_dir, f"{channel}.txt")
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


def combined_process_order():
    order = []
    for channel in CHANNELS:
        processes = ["signal"] + BKG if channel == "SR" else BKG[:]
        for process in processes:
            order.append((channel, process))
    return order


def write_flat_dummy_uncertainties(file_handle):
    process_order = combined_process_order()
    file_handle.write("\n")
    file_handle.write("# Flat dummy normalization uncertainties to absorb data/MC tension.\n")
    for nuisance_name, process_values in FLAT_DUMMY_UNCERTAINTIES:
        values = [process_values.get(process, "-") for _, process in process_order]
        file_handle.write(f"{nuisance_name:<18} lnN " + " ".join(values) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Build stat-only datacards from an input shapes ROOT file.")
    parser.add_argument("--shapes", required=True, help="Input shapes ROOT file")
    parser.add_argument("--tag", required=True, help="Tag used in output card names")
    parser.add_argument(
        "--flat-dummy-uncertainties",
        action="store_true",
        help="Add requested flat dummy normalization uncertainties as lnN nuisances.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    tag = sanitize_tag(args.tag)
    shapes_abs = os.path.abspath(args.shapes)

    cards_dir = os.path.join("cards", tag)
    os.makedirs(cards_dir, exist_ok=True)

    write_channel_card("SR", has_signal=True, shapes_abs=shapes_abs, cards_dir=cards_dir)
    write_channel_card("DYCR", has_signal=False, shapes_abs=shapes_abs, cards_dir=cards_dir)
    write_channel_card("EMUCR", has_signal=False, shapes_abs=shapes_abs, cards_dir=cards_dir)
    write_channel_card("CR3L", has_signal=False, shapes_abs=shapes_abs, cards_dir=cards_dir)

    combined = os.path.join("cards", f"combined_{tag}_statonly.txt")
    cmd = (
        f"combineCards.py SR={cards_dir}/SR.txt DYCR={cards_dir}/DYCR.txt "
        f"EMUCR={cards_dir}/EMUCR.txt CR3L={cards_dir}/CR3L.txt > {combined}"
    )
    subprocess.run(cmd, shell=True, check=True)
    normalize_combined_card_header(combined)

    with open(combined, "a", encoding="utf-8") as f:
        if args.flat_dummy_uncertainties:
            write_flat_dummy_uncertainties(f)

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
