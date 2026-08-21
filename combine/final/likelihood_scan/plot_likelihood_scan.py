#!/usr/bin/env python3
import csv
import math
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import uproot


BASE = os.environ.get(
    "MONOZ_SCAN_BASE",
    "/afs/cern.ch/user/l/liwe/hzz2l2nu/combine/final/likelihood_scan",
)
ROOT_DIR = os.path.join(BASE, "root")
PLOT_DIR = os.path.join(BASE, "plots")
TABLE_DIR = os.path.join(BASE, "tables")

SCAN_FILES = {
    "Expected_counting": os.path.join(
        ROOT_DIR, "higgsCombine.Expected_counting_grid.MultiDimFit.mH125.root"
    ),
    "Expected_shape": os.path.join(
        ROOT_DIR, "higgsCombine.Expected_shape_grid.MultiDimFit.mH125.root"
    ),
}

SINGLES_LOGS = {
    "Expected_counting": os.path.join(BASE, "logs", "counting_singles.log"),
    "Expected_shape": os.path.join(BASE, "logs", "shape_singles.log"),
}
ACTIVE_NAMES = (
    ["Expected_shape"]
    if os.environ.get("MONOZ_ONLY_SHAPE") == "1"
    else list(SCAN_FILES)
)


def read_scan(path):
    tree = uproot.open(path)["limit"]
    arr = tree.arrays(["r", "deltaNLL"], library="np")
    x = np.asarray(arr["r"], dtype=float)
    y = 2.0 * np.asarray(arr["deltaNLL"], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (y >= -1e-5) & (y < 20.0)
    x = x[mask]
    y = np.maximum(y[mask], 0.0)
    order = np.argsort(x)
    return x[order], y[order]


def crossing(x, y, level, side):
    pairs = list(zip(x, y))
    if side == "right":
        pairs = [(a, b) for a, b in pairs if a >= 0]
        pairs.sort()
    else:
        pairs = [(a, b) for a, b in pairs if a <= 0]
        pairs.sort(reverse=True)
    for (x1, y1), (x2, y2) in zip(pairs[:-1], pairs[1:]):
        if (y1 - level) * (y2 - level) <= 0 and y1 != y2:
            return x1 + (level - y1) * (x2 - x1) / (y2 - y1)
    return float("nan")


def parse_singles(log_path):
    text = open(log_path).read()
    m = re.search(
        r"r\s*:\s*([+-]?\d+(?:\.\d+)?)\s*"
        r"([+-]\d+(?:\.\d+)?)/([+-]\d+(?:\.\d+)?)",
        text,
    )
    if not m:
        return None
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


def fmt_result(best, lo, hi, ul):
    if abs(best) < 5e-4:
        best = 0.0
    down = abs(lo)
    up = abs(hi)
    return rf"{best:.3f}$^{{+{up:.3f}}}_{{-{down:.3f}}}$, 95% UL={ul:.3f}"


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(TABLE_DIR, exist_ok=True)

    colors = {"Expected_counting": "#1f77b4", "Expected_shape": "#ff7f0e"}
    labels = {"Expected_counting": "Expected_counting", "Expected_shape": "Expected_shape"}

    results = {}
    plt.figure(figsize=(10, 8))
    ax = plt.gca()

    for name in ACTIVE_NAMES:
        path = SCAN_FILES[name]
        x, y = read_scan(path)
        one_left = crossing(x, y, 1.0, "left")
        one_right = crossing(x, y, 1.0, "right")
        ul95 = crossing(x, y, 3.84, "right")
        best = x[np.argmin(y)]
        singles = parse_singles(SINGLES_LOGS[name])
        if singles is not None:
            best, err_lo, err_hi = singles
        else:
            err_lo = one_left - best
            err_hi = one_right - best
        if abs(best) < 5e-4:
            best = 0.0

        results[name] = {
            "best": best,
            "err_lo": err_lo,
            "err_hi": err_hi,
            "ul95": ul95,
            "one_left": one_left,
            "one_right": one_right,
        }

        ax.plot(x, y, color=colors[name], lw=2.2, marker="o", ms=3.0, alpha=0.95,
                label=f"{labels[name]}: {fmt_result(best, err_lo, err_hi, ul95)}")
        ax.plot([ul95, ul95], [0.0, 3.84], color=colors[name], ls="--", lw=2.0, alpha=0.45)
        ax.plot([ul95], [3.84], color=colors[name], marker="o", ms=9)

    ax.axhline(1.0, color="red", ls=":", lw=1.2, alpha=0.65)
    ax.axhline(3.84, color="orange", ls=":", lw=1.2, alpha=0.8)
    ax.set_xlim(-1.0, 2.0)
    ax.set_ylim(0.0, 6.0)
    ax.set_xlabel(r"$BR_{inv}$", fontsize=18)
    ax.set_ylabel(r"$-2\Delta\ln L$", fontsize=18)
    ax.grid(True, ls="--", alpha=0.3)
    ax.text(0.0, 1.01, r"$\bf{CMS}$ $\it{Internal}$",
            transform=ax.transAxes, fontsize=18, va="bottom")
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, fontsize=12)
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()

    stem = (
        "likelihood_scan_expected_shape"
        if ACTIVE_NAMES == ["Expected_shape"]
        else "likelihood_scan_expected_counting_shape"
    )
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(PLOT_DIR, f"{stem}.{ext}"))

    with open(os.path.join(TABLE_DIR, "likelihood_scan_summary.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["name", "best", "err_lo", "err_hi", "ul95", "one_left", "one_right"],
        )
        writer.writeheader()
        for name, values in results.items():
            writer.writerow({"name": name, **values})


if __name__ == "__main__":
    main()
