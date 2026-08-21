#!/usr/bin/env python3

"""Plot a compact expected/blinded CLs upper-limit summary."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import uproot


def read_limits(path):
    with uproot.open(path) as root_file:
        tree = root_file["limit"]
        limits = tree["limit"].array(library="np")
        quantiles = tree["quantileExpected"].array(library="np")
    return {
        float(quantile): float(limit)
        for quantile, limit in zip(quantiles, limits)
    }


def value_at(values, target):
    key = min(values, key=lambda candidate: abs(candidate - target))
    return values[key]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", required=True)
    parser.add_argument("--blinded", required=True)
    parser.add_argument("--output", required=True, help="Output path stem")
    args = parser.parse_args()

    expected = read_limits(args.expected)
    blinded = read_limits(args.blinded)
    q025 = value_at(expected, 0.025)
    q16 = value_at(expected, 0.16)
    q50 = value_at(expected, 0.50)
    q84 = value_at(expected, 0.84)
    q975 = value_at(expected, 0.975)
    blinded_value = value_at(blinded, -1.0)

    fig, axis = plt.subplots(figsize=(8.2, 4.3))
    axis.errorbar(
        q50, 1.0,
        xerr=np.array([[q50 - q025], [q975 - q50]]),
        fmt="none", ecolor="#F5C242", elinewidth=18, capsize=0,
        label=r"Expected $\pm2\sigma$",
    )
    axis.errorbar(
        q50, 1.0,
        xerr=np.array([[q50 - q16], [q84 - q50]]),
        fmt="none", ecolor="#59B85C", elinewidth=18, capsize=0,
        label=r"Expected $\pm1\sigma$",
    )
    axis.plot(q50, 1.0, marker="|", color="black", markersize=25,
              markeredgewidth=2.2, label="Expected median")
    axis.plot(blinded_value, 0.0, marker="o", color="#2B6CB0",
              markersize=8, label="Blinded result")

    axis.set_yticks([0.0, 1.0])
    axis.set_yticklabels([
        "SR Asimov + CR data",
        "Background-only expected",
    ])
    axis.set_ylim(-0.6, 1.6)
    axis.set_xlim(left=0.0)
    axis.set_xlabel(r"95% CL upper limit on $r$ ($\mathrm{BR}_{inv}$)")
    axis.set_title(r"Combined 2017+2018 $p_\mathrm{T}^{miss}$ shape analysis")
    axis.grid(axis="x", color="0.85", linewidth=0.8)
    axis.legend(frameon=False, loc="upper right", fontsize=9)
    fig.tight_layout()

    fig.savefig(f"{args.output}.png", dpi=180)
    fig.savefig(f"{args.output}.pdf")
    print(f"Wrote {args.output}.png")
    print(f"Wrote {args.output}.pdf")


if __name__ == "__main__":
    main()
