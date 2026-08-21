#!/usr/bin/env python3

import argparse
import os

import uproot


LABELS = {
    -1.0: "observed",
    0.025: "expected -2sigma",
    0.16: "expected -1sigma",
    0.5: "expected median",
    0.84: "expected +1sigma",
    0.975: "expected +2sigma",
}


def label_for_quantile(value):
    for quantile, label in LABELS.items():
        if abs(value - quantile) < 1e-4:
            return label
    return f"quantile {value:g}"


def read_limits(path):
    with uproot.open(path) as root_file:
        tree = root_file["limit"]
        quantiles = tree["quantileExpected"].array(library="np")
        limits = tree["limit"].array(library="np")
    return [(float(q), float(l)) for q, l in zip(quantiles, limits)]


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize Combine AsymptoticLimits ROOT files.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--observed", required=True)
    parser.add_argument("--expected", required=True)
    parser.add_argument("--outdir", default="limits")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for source, path in [("observed file", args.observed), ("expected file", args.expected)]:
        for quantile, limit in read_limits(path):
            rows.append((source, quantile, label_for_quantile(quantile), limit))

    out_path = os.path.join(args.outdir, f"limits_{args.tag}.txt")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(f"# {args.tag}\n")
        for source, quantile, label, limit in rows:
            handle.write(f"{source:13s} {quantile:8.3f} {label:18s} {limit:.6f}\n")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
