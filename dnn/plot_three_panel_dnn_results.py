#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot


BASE = Path("/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn")
OUT = BASE / "three_panel_results"
ROOT_DIR = OUT / "root"
LOG_DIR = OUT / "logs"
PLOT_DIR = OUT / "plots"
TABLE_DIR = OUT / "tables"

DATASETS = [
    {
        "key": "dnn2017",
        "title": "2017",
        "limit": BASE / "2017/limits/higgsCombine.Expected_dnn_score_2017.AsymptoticLimits.mH125.root",
        "params": ["k_Zjet", "k_WZ", "k_emu"],
    },
    {
        "key": "dnn2018",
        "title": "2018",
        "limit": BASE / "2018/limits/higgsCombine.Expected_dnn_score_2018.AsymptoticLimits.mH125.root",
        "params": ["k_Zjet", "k_WZ", "k_emu"],
    },
    {
        "key": "dnn2017_2018",
        "title": "2017+2018",
        "limit": BASE / "combined_2017_2018/limits/higgsCombine.Expected_dnn_score_2017_2018.AsymptoticLimits.mH125.root",
        "params": ["k_Zjet", "k_WZ", "k_emu"],
    },
]


def ensure_dirs() -> None:
    for path in [PLOT_DIR, TABLE_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def read_limit(path: Path) -> dict[str, float]:
    arr = uproot.open(path)["limit"].arrays(["limit", "quantileExpected"], library="np")
    values = {}
    for limit, quantile in zip(arr["limit"], arr["quantileExpected"]):
        q = float(quantile)
        if q < 0:
            values["observed_asimov"] = float(limit)
        elif abs(q - 0.025) < 1e-4:
            values["exp_m2"] = float(limit)
        elif abs(q - 0.16) < 1e-4:
            values["exp_m1"] = float(limit)
        elif abs(q - 0.5) < 1e-4:
            values["exp"] = float(limit)
        elif abs(q - 0.84) < 1e-4:
            values["exp_p1"] = float(limit)
        elif abs(q - 0.975) < 1e-4:
            values["exp_p2"] = float(limit)
    return values


def read_scan(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arr = uproot.open(path)["limit"].arrays(["r", "deltaNLL"], library="np")
    x = np.asarray(arr["r"], dtype=float)
    y = 2.0 * np.asarray(arr["deltaNLL"], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (y >= -1e-5) & (y < 25.0)
    x = x[mask]
    y = np.maximum(y[mask], 0.0)
    order = np.argsort(x)
    return x[order], y[order]


def parse_singles(log_path: Path) -> tuple[float, float, float] | None:
    if not log_path.exists():
        return None
    text = log_path.read_text(errors="ignore")
    match = re.search(
        r"\br\s*:\s*([+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)\s*([+-]\d+(?:\.\d+)?(?:e[+-]?\d+)?)/([+-]\d+(?:\.\d+)?(?:e[+-]?\d+)?)",
        text,
        re.IGNORECASE,
    )
    if not match:
        return None
    best = float(match.group(1))
    lo = float(match.group(2))
    hi = float(match.group(3))
    if abs(best) < 5e-4:
        best = 0.0
    return best, lo, hi


def crossing(x: np.ndarray, y: np.ndarray, level: float, side: str, reference: float = 0.0) -> float:
    pairs = list(zip(x, y))
    if side == "right":
        pairs = [(a, b) for a, b in pairs if a >= reference]
        pairs.sort()
    else:
        pairs = [(a, b) for a, b in pairs if a <= reference]
        pairs.sort(reverse=True)
    for (x1, y1), (x2, y2) in zip(pairs[:-1], pairs[1:]):
        if (y1 - level) * (y2 - level) <= 0 and y1 != y2:
            return float(x1 + (level - y1) * (x2 - x1) / (y2 - y1))
    return float("nan")


def read_impact(key: str, param: str) -> dict:
    path = ROOT_DIR / f"higgsCombine_paramFit_.{key}_impact_{param}.MultiDimFit.mH125.root"
    arr = uproot.open(path)["limit"].arrays([param, "r"], library="np")
    p = [float(v) for v in arr[param]]
    r = [float(v) for v in arr["r"]]
    if len(p) < 3 or len(r) < 3:
        raise RuntimeError(f"Unexpected impact entries in {path}")
    fit = [p[1], p[0], p[2]]
    rfit = [r[1], r[0], r[2]]
    impact = max(abs(rfit[0] - rfit[1]), abs(rfit[2] - rfit[1]))
    return {
        "name": param,
        "display_name": param.replace("_", " "),
        "components": [],
        "fit": fit,
        "r": rfit,
        "impact": impact,
    }


def with_single_component(impact: dict) -> dict:
    impact = dict(impact)
    impact["components"] = [{
        "name": impact["name"],
        "display_name": impact["display_name"],
        "fit": impact["fit"],
        "r": impact["r"],
    }]
    return impact


def draw_cls() -> dict[str, dict[str, float]]:
    values_by_key = {item["key"]: read_limit(item["limit"]) for item in DATASETS}
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.0), sharey=True)
    xmax = max(values["exp_p2"] for values in values_by_key.values()) * 1.18

    for ax, item in zip(axes, DATASETS):
        values = values_by_key[item["key"]]
        y = 0.0
        ax.barh(y, values["exp_p2"] - values["exp_m2"], left=values["exp_m2"], height=0.48, color="#f4d03f")
        ax.barh(y, values["exp_p1"] - values["exp_m1"], left=values["exp_m1"], height=0.30, color="#58d68d")
        ax.plot([values["exp"], values["exp"]], [y - 0.26, y + 0.26], color="black", lw=2)
        ax.plot(values["observed_asimov"], y, marker="o", ms=6, color="#ff7f0e")
        ax.text(values["exp"], y + 0.36, f"{values['exp']:.3f}", ha="center", va="bottom", fontsize=11)
        ax.set_title(item["title"], fontsize=14)
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.75, 0.75)
        ax.grid(axis="x", ls="--", alpha=0.32)
        ax.set_xlabel(r"95% CL upper limit on $BR_{inv}$", fontsize=12)
        ax.set_yticks([0])
        ax.set_yticklabels(["DNN score"])
    axes[0].text(0.0, 1.10, r"$\bf{CMS}$ $\it{Internal}$", transform=axes[0].transAxes, fontsize=17)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color="#f4d03f", label=r"Expected $\pm2\sigma$"),
        plt.Rectangle((0, 0), 1, 1, color="#58d68d", label=r"Expected $\pm1\sigma$"),
        plt.Line2D([0], [0], color="black", lw=2, label="Expected median"),
        plt.Line2D([0], [0], color="#ff7f0e", marker="o", lw=0, label="Asimov observed"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=11)
    fig.tight_layout(rect=[0, 0.12, 1, 1])
    fig.savefig(PLOT_DIR / "cls_limit_2017_2018_combined.png", dpi=170)
    fig.savefig(PLOT_DIR / "cls_limit_2017_2018_combined.pdf")
    plt.close(fig)
    return values_by_key


def draw_scan() -> list[dict]:
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    colors = ["#1f77b4", "#2ca02c", "#d62728"]
    for ax, item, color in zip(axes, DATASETS, colors):
        key = item["key"]
        x, y = read_scan(ROOT_DIR / f"higgsCombine.{key}_grid.MultiDimFit.mH125.root")
        singles = parse_singles(LOG_DIR / f"{key}_singles.log")
        if singles:
            best, err_lo, err_hi = singles
        else:
            best = float(x[np.argmin(y)])
            err_lo = crossing(x, y, 1.0, "left", best) - best
            err_hi = crossing(x, y, 1.0, "right", best) - best
        ul95 = crossing(x, y, 3.84, "right", max(best, 0.0))
        ax.plot(x, y, color=color, lw=2.2, marker="o", ms=2.5)
        ax.axhline(1.0, color="red", ls=":", lw=1.1, alpha=0.65)
        ax.axhline(3.84, color="orange", ls=":", lw=1.1, alpha=0.8)
        if math.isfinite(ul95):
            ax.plot([ul95, ul95], [0, 3.84], color=color, ls="--", lw=1.6, alpha=0.55)
            ax.plot([ul95], [3.84], color=color, marker="o", ms=6)
        ax.set_title(item["title"], fontsize=14)
        ax.set_xlim(-0.3, 0.7)
        ax.set_ylim(0, 6)
        ax.set_xlabel(r"$BR_{inv}$", fontsize=12)
        ax.grid(True, ls="--", alpha=0.3)
        label = rf"$\hat{{r}}={best:.3f}^{{+{abs(err_hi):.3f}}}_{{-{abs(err_lo):.3f}}}$"
        if math.isfinite(ul95):
            label += "\n" + rf"95% crossing={ul95:.3f}"
        ax.text(0.04, 0.95, label, transform=ax.transAxes, va="top", fontsize=11)
        rows.append({
            "key": key,
            "title": item["title"],
            "best": best,
            "err_lo": err_lo,
            "err_hi": err_hi,
            "ul95": ul95,
        })
    axes[0].set_ylabel(r"$-2\Delta\ln L$", fontsize=13)
    axes[0].text(0.0, 1.10, r"$\bf{CMS}$ $\it{Internal}$", transform=axes[0].transAxes, fontsize=17)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "likelihood_scan_2017_2018_combined.png", dpi=170)
    fig.savefig(PLOT_DIR / "likelihood_scan_2017_2018_combined.pdf")
    plt.close(fig)
    return rows


def draw_ranking() -> dict[str, list[dict]]:
    all_impacts = {}
    for item in DATASETS:
        raw_impacts = [read_impact(item["key"], param) for param in item["params"]]
        impacts = [with_single_component(impact) for impact in raw_impacts]
        impacts.sort(key=lambda p: p["impact"], reverse=True)
        all_impacts[item["key"]] = impacts

    max_impact = max(
        max(abs(p["r"][0] - p["r"][1]), abs(p["r"][2] - p["r"][1]))
        for impacts in all_impacts.values()
        for p in impacts
    )
    all_k_values = [
        value
        for impacts in all_impacts.values()
        for p in impacts
        for component in p["components"]
        for value in component["fit"]
    ]
    k_min = min(0.95, min(all_k_values) - 0.06)
    k_max = max(1.05, max(all_k_values) + 0.06)
    impact_lim = max_impact * 1.55 if max_impact > 0 else 0.1

    fig = plt.figure(figsize=(13.5, 10.2))
    grid = fig.add_gridspec(
        3,
        3,
        width_ratios=[1.15, 4.6, 3.0],
        height_ratios=[3, 3, 3],
        hspace=0.42,
        wspace=0.06,
    )

    fit_handle = None
    prefit_handle = None
    plus_handle = None
    minus_handle = None

    for row, item in enumerate(DATASETS):
        impacts = all_impacts[item["key"]]
        nrows = len(impacts)
        y = np.arange(nrows)
        rank_ax = fig.add_subplot(grid[row, 0])
        fit_ax = fig.add_subplot(grid[row, 1])
        impact_ax = fig.add_subplot(grid[row, 2], sharey=fit_ax)

        for ax in (rank_ax, fit_ax, impact_ax):
            ax.set_ylim(nrows - 0.5, -0.5)
            for idx in range(nrows):
                if idx % 2 == 0:
                    ax.axhspan(idx - 0.5, idx + 0.5, color="#e8e8e8", zorder=0)

        rank_ax.set_xlim(0, 1)
        rank_ax.axis("off")
        labels = [p["display_name"] for p in impacts]
        for idx, label in enumerate(labels):
            rank_ax.text(0.20, idx, str(idx + 1), ha="center", va="center", fontsize=12, fontweight="bold")
            rank_ax.text(0.96, idx, label, ha="right", va="center", fontsize=10, color="#6f83a5")

        fit_ax.set_yticks(y)
        fit_ax.set_yticklabels([])
        fit_ax.tick_params(axis="y", length=0)
        fit_ax.set_xlim(k_min, k_max)
        fit_ax.axvline(1.0, color="black", lw=1.2)
        fit_ax.grid(axis="x", ls=":", color="0.55", alpha=0.7)
        fit_ax.set_xlabel("k-factor value", fontsize=11)
        fit_ax.set_title(item["title"], fontsize=13, pad=8)

        impact_ax.set_xlim(-impact_lim, impact_lim)
        impact_ax.axvline(0.0, color="black", lw=1.2)
        impact_ax.grid(axis="x", ls=":", color="0.55", alpha=0.7)
        impact_ax.tick_params(axis="y", labelleft=False)
        impact_ax.set_xlabel(r"$\Delta r$", fontsize=14, loc="right")

        singles = parse_singles(LOG_DIR / f"{item['key']}_singles.log")
        if singles is not None:
            best, err_lo, err_hi = singles
            title = rf"$\hat{{r}}={best:.2f}^{{+{abs(err_hi):.2f}}}_{{-{abs(err_lo):.2f}}}$"
            impact_ax.set_title(title, fontsize=13, pad=8)

        for yy, impact in zip(y, impacts):
            offsets = [0.0] if len(impact["components"]) == 1 else [-0.16, 0.16]
            for component, offset in zip(impact["components"], offsets):
                k_low, k_best, k_high = component["fit"]
                xerr = np.array([[k_best - k_low], [k_high - k_best]])
                fit_handle = fit_ax.errorbar(
                    k_best,
                    yy + offset,
                    xerr=xerr,
                    fmt="o",
                    color="black",
                    ecolor="black",
                    elinewidth=1.2,
                    capsize=3.5,
                    ms=4.8,
                    zorder=4,
                    label="Fit",
                )
                prefit_handle = fit_ax.plot(
                    1.0,
                    yy + offset,
                    marker="x",
                    color="blue",
                    ms=6,
                    mew=1.5,
                    linestyle="None",
                    zorder=5,
                    label="Prefit",
                )[0]
                prefix = "" if len(impact["components"]) == 1 else f"{component['display_name']}: "
                fit_ax.text(
                    min(k_best + 0.04, k_max - 0.02),
                    yy + offset - 0.10,
                    prefix + rf"${k_best:.2f}^{{+{k_high - k_best:.2f}}}_{{-{k_best - k_low:.2f}}}$",
                    fontsize=8.5,
                    ha="left",
                    va="center",
                )

            r_low, r_best, r_high = impact["r"]
            minus_delta = r_low - r_best
            plus_delta = r_high - r_best
            plus_handle = impact_ax.barh(
                yy,
                plus_delta,
                left=0.0,
                height=0.72,
                color="#d98291",
                edgecolor="#d98291",
                alpha=0.95,
                label=r"$+1\sigma$ Impact",
            )
            minus_handle = impact_ax.barh(
                yy,
                minus_delta,
                left=0.0,
                height=0.72,
                color="#9bb7d9",
                edgecolor="#9bb7d9",
                alpha=0.95,
                label=r"$-1\sigma$ Impact",
            )
            impact_ax.text(
                impact_lim * 0.98,
                yy,
                f"{impact['impact']:.3f}",
                ha="right",
                va="center",
                fontsize=9,
            )

        if row == 0:
            fig.text(0.37, 0.965, r"$\bf{CMS}$ $\it{Internal}$", ha="center", fontsize=18)
            fit_ax.text(0.95, 1.16, "Expected DNN", transform=fit_ax.transAxes, ha="right", fontsize=13)

    handles = [
        fit_handle,
        prefit_handle,
        plus_handle,
        minus_handle,
    ]
    labels = ["Fit", "Prefit", r"$+1\sigma$ Impact", r"$-1\sigma$ Impact"]
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.08, 0.01), ncol=2, frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.09, top=0.90)
    fig.savefig(PLOT_DIR / "ranking_2017_2018_combined.png", dpi=170)
    fig.savefig(PLOT_DIR / "ranking_2017_2018_combined.pdf")
    plt.close(fig)
    return all_impacts


def write_tables(cls_values: dict, scan_rows: list[dict], impacts: dict) -> None:
    with open(TABLE_DIR / "cls_limit_2017_2018_combined.json", "w") as handle:
        json.dump(cls_values, handle, indent=2, sort_keys=True)
    with open(TABLE_DIR / "likelihood_scan_2017_2018_combined.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["key", "title", "best", "err_lo", "err_hi", "ul95"])
        writer.writeheader()
        writer.writerows(scan_rows)
    with open(TABLE_DIR / "ranking_2017_2018_combined.json", "w") as handle:
        json.dump(impacts, handle, indent=2, sort_keys=True)


def main() -> None:
    ensure_dirs()
    cls_values = draw_cls()
    scan_rows = draw_scan()
    impacts = draw_ranking()
    write_tables(cls_values, scan_rows, impacts)
    print(f"Wrote plots to {PLOT_DIR}")


if __name__ == "__main__":
    main()
