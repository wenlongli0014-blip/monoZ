#!/usr/bin/env python3

import argparse
import os

BACKGROUNDS = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
CHANNELS = ["SR", "DYCR", "EMUCR", "CR3L"]

K_FACTOR_PROCESSES = {
    "k_Zjet": ["DY"],
    "k_WZ": ["WZ"],
    "k_emu": ["ttbar", "WW", "ST"],
}


def sanitize_tag(text):
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def process_list(channel, dummy_signal):
    if channel == "SR" or dummy_signal:
        return ["signal"] + BACKGROUNDS
    return BACKGROUNDS[:]


def process_ids(channel, dummy_signal):
    if channel == "SR" or dummy_signal:
        return [0] + list(range(1, len(BACKGROUNDS) + 1))
    return list(range(1, len(BACKGROUNDS) + 1))


def parse_args():
    parser = argparse.ArgumentParser(description="Build a combined datacard for one fit-check model.")
    parser.add_argument("--shapes", required=True, help="Input shapes ROOT file")
    parser.add_argument("--tag", required=True, help="Output tag")
    parser.add_argument("--channels", nargs="+", required=True, choices=CHANNELS)
    parser.add_argument(
        "--k-factors",
        nargs="*",
        default=[],
        choices=sorted(K_FACTOR_PROCESSES),
        help="Floating normalization factors to add as rateParams.",
    )
    parser.add_argument(
        "--out-dir",
        default="cards",
        help="Directory for the output datacard, relative to the current directory.",
    )
    parser.add_argument(
        "--dummy-signal",
        action="store_true",
        help="Include a negligible signal template in non-SR cards as a Combine technical placeholder.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    tag = sanitize_tag(args.tag)
    channels = args.channels
    k_factors = args.k_factors
    shapes_abs = os.path.abspath(args.shapes)
    os.makedirs(args.out_dir, exist_ok=True)
    card_path = os.path.join(args.out_dir, f"{tag}.txt")

    entries = []
    for channel in channels:
        use_dummy = args.dummy_signal and channel != "SR"
        for process, proc_id in zip(process_list(channel, use_dummy), process_ids(channel, use_dummy)):
            entries.append((channel, process, proc_id))

    lines = [
        f"# Fit-check datacard: {tag}",
        "# No flat lnN uncertainties, no manual shape uncertainties, no autoMCStats.",
        f"# Channels: {' '.join(channels)}",
        f"# Floating k-factors: {' '.join(k_factors) if k_factors else 'none'}",
        f"imax {len(channels)}",
        "jmax *",
        "kmax *",
        "------------",
    ]

    for channel in channels:
        lines.append(f"shapes * {channel} {shapes_abs} {channel}/$PROCESS")

    lines.extend(
        [
            "------------",
            "bin " + " ".join(channels),
            "observation " + " ".join(["-1"] * len(channels)),
            "------------",
            "bin " + " ".join(channel for channel, _, _ in entries),
            "process " + " ".join(process for _, process, _ in entries),
            "process " + " ".join(str(proc_id) for _, _, proc_id in entries),
            "rate " + " ".join(["-1"] * len(entries)),
            "------------",
        ]
    )

    if k_factors:
        lines.append("# Floating normalization parameters. They act on the named processes in every selected channel.")
        for k_factor in k_factors:
            for process in K_FACTOR_PROCESSES[k_factor]:
                lines.append(f"{k_factor:<8s} rateParam * {process:<5s} 1 [0,5]")

    lines.append("")

    with open(card_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    print(card_path)


if __name__ == "__main__":
    main()
