#!/usr/bin/env python3

import argparse
import os
import subprocess


BKG = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
CHANNELS = ["SR", "DYCR", "EMUCR", "CR3L"]
MANUAL_SHAPES = [
    ("manual_DY_40", {"DY"}),
    ("manual_ZZ_15", {"ggZZ", "qqZZ"}),
    ("manual_WZ_10", {"WZ"}),
    ("manual_ttbar_15", {"ttbar"}),
    ("manual_WW_15", {"WW"}),
    ("manual_ST_15", {"ST"}),
]


def sanitize_tag(text):
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def write_channel_card(channel, has_signal, shapes_abs, cards_dir, manual_shapes):
    if has_signal:
        processes = ["signal"] + BKG
        proc_ids = [0] + list(range(1, len(BKG) + 1))
    else:
        processes = BKG[:]
        proc_ids = list(range(1, len(BKG) + 1))

    lines = [
        "imax 1",
        "jmax *",
        "kmax *",
        "------------",
        f"shapes * {channel} {shapes_abs} {channel}/$PROCESS {channel}/$PROCESS_$SYSTEMATIC",
        "------------",
        f"bin {channel}",
        "observation -1",
        "------------",
        "bin " + " ".join([channel] * len(processes)),
        "process " + " ".join(processes),
        "process " + " ".join(str(x) for x in proc_ids),
        "rate " + " ".join(["-1"] * len(processes)),
        "------------",
    ]

    for nuisance, affected in manual_shapes:
        values = ["1" if process in affected else "-" for process in processes]
        lines.append(f"{nuisance:<16s} shape " + " ".join(values))

    lines.extend(["* autoMCStats 0", ""])

    out_path = os.path.join(cards_dir, f"{channel}.txt")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def normalize_combined_card_header(card_path):
    with open(card_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    for i_line, line in enumerate(lines):
        if line.startswith("kmax "):
            lines[i_line] = "kmax * number of nuisance parameters\n"
            break

    with open(card_path, "w", encoding="utf-8") as handle:
        handle.writelines(lines)


def parse_args():
    parser = argparse.ArgumentParser(description="Build datacards with manual shape uncertainties.")
    parser.add_argument("--shapes", required=True, help="Input shapes ROOT file")
    parser.add_argument("--tag", required=True, help="Tag used in output card names")
    parser.add_argument(
        "--exclude-manual-shapes",
        nargs="*",
        default=[],
        help="Manual shape nuisance names to omit from the datacards, e.g. manual_DY_40.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    tag = sanitize_tag(args.tag)
    shapes_abs = os.path.abspath(args.shapes)
    cards_dir = os.path.join("cards", tag)
    os.makedirs(cards_dir, exist_ok=True)

    excluded = set(args.exclude_manual_shapes)
    manual_shapes = [(name, affected) for name, affected in MANUAL_SHAPES if name not in excluded]
    unknown = excluded - {name for name, _ in MANUAL_SHAPES}
    if unknown:
        raise ValueError(f"Unknown manual shape nuisance(s): {', '.join(sorted(unknown))}")

    write_channel_card("SR", has_signal=True, shapes_abs=shapes_abs, cards_dir=cards_dir, manual_shapes=manual_shapes)
    write_channel_card("DYCR", has_signal=False, shapes_abs=shapes_abs, cards_dir=cards_dir, manual_shapes=manual_shapes)
    write_channel_card("EMUCR", has_signal=False, shapes_abs=shapes_abs, cards_dir=cards_dir, manual_shapes=manual_shapes)
    write_channel_card("CR3L", has_signal=False, shapes_abs=shapes_abs, cards_dir=cards_dir, manual_shapes=manual_shapes)

    combined = os.path.join("cards", f"combined_{tag}_manualshape.txt")
    cmd = (
        f"combineCards.py SR={cards_dir}/SR.txt DYCR={cards_dir}/DYCR.txt "
        f"EMUCR={cards_dir}/EMUCR.txt CR3L={cards_dir}/CR3L.txt > {combined}"
    )
    subprocess.run(cmd, shell=True, check=True)
    normalize_combined_card_header(combined)

    with open(combined, "a", encoding="utf-8") as handle:
        handle.write("\n")
        handle.write("# Shared normalization parameters.\n")
        handle.write("k_DY    rateParam SR    DY    1 [0,5]\n")
        handle.write("k_DY    rateParam DYCR  DY    1 [0,5]\n")
        handle.write("\n")
        handle.write("k_WZ    rateParam SR    WZ    1 [0,5]\n")
        handle.write("k_WZ    rateParam CR3L  WZ    1 [0,5]\n")
        handle.write("\n")
        handle.write("k_TT    rateParam SR    ttbar 1 [0,5]\n")
        handle.write("k_TT    rateParam EMUCR ttbar 1 [0,5]\n")

    print(f"Wrote {combined}")


if __name__ == "__main__":
    main()
