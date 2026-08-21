#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot


BASE = "/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/likelihood_scan"
ROOT_DIR = os.path.join(BASE, "root")
PLOT_DIR = os.path.join(BASE, "plots")
TABLE_DIR = os.path.join(BASE, "tables")


def read_scan(path):
    arr = uproot.open(path)["limit"].arrays(["r", "deltaNLL"], library="np")
    x = np.asarray(arr["r"], dtype=float)
    y = 2.0 * np.asarray(arr["deltaNLL"], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (y >= -1e-5) & (y < 20.0)
    x = x[mask]
    y = np.maximum(y[mask], 0.0)
    order = np.argsort(x)
    return x[order], y[order]


def crossing(x, y, level, side, reference=0.0):
    pairs = list(zip(x, y))
    if side == "right":
        pairs = [(a, b) for a, b in pairs if a >= reference]
        pairs.sort()
    else:
        pairs = [(a, b) for a, b in pairs if a <= reference]
        pairs.sort(reverse=True)
    for (x1, y1), (x2, y2) in zip(pairs[:-1], pairs[1:]):
        if (y1 - level) * (y2 - level) <= 0 and y1 != y2:
            return x1 + (level - y1) * (x2 - x1) / (y2 - y1)
    return float("nan")


def parse_singles(log_path):
    text = open(log_path, encoding="utf-8").read()
    match = re.search(
        r"\br\s*:\s*([+-]?\d+(?:\.\d+)?)\s*([+-]\d+(?:\.\d+)?)/([+-]\d+(?:\.\d+)?)",
        text,
    )
    if not match:
        return None
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def fmt_result(best, lo, hi, ul):
    if abs(best) < 5e-4:
        best = 0.0
    if np.isfinite(ul):
        suffix = rf", 95% crossing={ul:.3f}"
    else:
        suffix = r", no positive 95% crossing"
    return rf"{best:.3f}$^{{+{abs(hi):.3f}}}_{{-{abs(lo):.3f}}}$" + suffix


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(TABLE_DIR, exist_ok=True)

    scan_path = os.path.join(ROOT_DIR, "higgsCombine.Expected_dnn_score_grid.MultiDimFit.mH125.root")
    singles_log = os.path.join(BASE, "logs", "dnn_score_singles.log")
    x, y = read_scan(scan_path)
    best = float(x[np.argmin(y)])
    one_left = crossing(x, y, 1.0, "left", reference=best)
    one_right = crossing(x, y, 1.0, "right", reference=best)
    ul95_unconstrained = crossing(x, y, 3.84, "right", reference=best)
    ul95_physical = crossing(x, y, 3.84, "right", reference=max(best, 0.0))
    singles = parse_singles(singles_log)
    if singles is not None:
        best, err_lo, err_hi = singles
    else:
        err_lo = one_left - best
        err_hi = one_right - best
    if abs(best) < 5e-4:
        best = 0.0

    fig, ax = plt.subplots(figsize=(9.5, 7.0))
    ax.plot(x, y, color="#ff7f0e", lw=2.2, marker="o", ms=3.0, alpha=0.95, label=f"DNN score: {fmt_result(best, err_lo, err_hi, ul95_physical)}")
    if np.isfinite(ul95_unconstrained):
        ax.plot([ul95_unconstrained, ul95_unconstrained], [0.0, 3.84], color="#ff7f0e", ls="--", lw=2.0, alpha=0.45)
        ax.plot([ul95_unconstrained], [3.84], color="#ff7f0e", marker="o", ms=9)
    ax.axhline(1.0, color="red", ls=":", lw=1.2, alpha=0.65)
    ax.axhline(3.84, color="orange", ls=":", lw=1.2, alpha=0.8)
    ax.set_xlim(-1.0, 2.0)
    ax.set_ylim(0.0, 6.0)
    ax.set_xlabel(r"$BR_{inv}$", fontsize=18)
    ax.set_ylabel(r"$-2\Delta\ln L$", fontsize=18)
    ax.grid(True, ls="--", alpha=0.3)
    ax.text(0.0, 1.01, r"$\bf{CMS}$ $\it{Internal}$", transform=ax.transAxes, fontsize=18, va="bottom")
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, fontsize=12)
    ax.tick_params(axis="both", labelsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "likelihood_scan_expected_dnn_score.png"), dpi=150)
    fig.savefig(os.path.join(PLOT_DIR, "likelihood_scan_expected_dnn_score.pdf"))
    plt.close(fig)

    with open(os.path.join(TABLE_DIR, "likelihood_scan_dnn_score_summary.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "name",
            "best",
            "err_lo",
            "err_hi",
            "ul95_unconstrained",
            "ul95_physical_positive",
            "one_left",
            "one_right",
        ])
        writer.writeheader()
        writer.writerow({
            "name": "Expected_dnn_score",
            "best": best,
            "err_lo": err_lo,
            "err_hi": err_hi,
            "ul95_unconstrained": ul95_unconstrained,
            "ul95_physical_positive": ul95_physical,
            "one_left": one_left,
            "one_right": one_right,
        })


if __name__ == "__main__":
    main()
