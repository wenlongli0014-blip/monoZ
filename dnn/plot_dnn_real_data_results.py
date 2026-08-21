#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


BASE = "/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/real_data_fit"
ROOT_DIR = os.path.join(BASE, "root")
PLOT_DIR = os.path.join(BASE, "plots")
TABLE_DIR = os.path.join(BASE, "tables")
JSON_DIR = os.path.join(BASE, "json")
PARAMS = ["k_Zjet", "k_WZ", "k_emu"]


def parse_singles(log_path):
    text = open(log_path, encoding="utf-8").read()
    match = re.search(
        r"\br\s*:\s*([+-]?\d+(?:\.\d+)?)\s*([+-]\d+(?:\.\d+)?)/([+-]\d+(?:\.\d+)?)",
        text,
    )
    if not match:
        return None
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def load_limit():
    path = os.path.join(ROOT_DIR, "higgsCombine.RealData_dnn_score.AsymptoticLimits.mH125.root")
    arr = uproot.open(path)["limit"].arrays(["limit", "quantileExpected"], library="np")
    out = {}
    for value, quantile in zip(arr["limit"], arr["quantileExpected"]):
        q = float(quantile)
        if q < 0:
            out["observed"] = float(value)
        elif abs(q - 0.025) < 1e-4:
            out["exp_m2"] = float(value)
        elif abs(q - 0.16) < 1e-4:
            out["exp_m1"] = float(value)
        elif abs(q - 0.5) < 1e-4:
            out["exp"] = float(value)
        elif abs(q - 0.84) < 1e-4:
            out["exp_p1"] = float(value)
        elif abs(q - 0.975) < 1e-4:
            out["exp_p2"] = float(value)
    return out


def draw_limit(values):
    fig, ax = plt.subplots(figsize=(8.8, 3.6))
    y = 0.0
    ax.barh(y, values["exp_p2"] - values["exp_m2"], left=values["exp_m2"], height=0.46, color="#f4d03f", label=r"Expected $\pm2\sigma$")
    ax.barh(y, values["exp_p1"] - values["exp_m1"], left=values["exp_m1"], height=0.30, color="#58d68d", label=r"Expected $\pm1\sigma$")
    ax.plot([values["exp"], values["exp"]], [y - 0.24, y + 0.24], color="black", lw=2.0, label="Expected median")
    ax.plot(values["observed"], y, marker="o", ms=7, color="#ff7f0e", label="Observed data")
    ax.text(values["exp_p2"] + 0.035, y, f"obs={values['observed']:.3f}, exp={values['exp']:.3f}", va="center", ha="left", fontsize=12)
    ax.set_yticks([0])
    ax.set_yticklabels(["DNN score"], fontsize=13)
    ax.set_xlabel(r"95% CL upper limit on $BR_{inv}$", fontsize=15)
    ax.set_xlim(0, max(values["exp_p2"], values["observed"]) * 1.38)
    ax.set_ylim(-0.65, 0.65)
    ax.grid(axis="x", ls="--", alpha=0.3)
    ax.text(0.0, 1.02, r"$\bf{CMS}$ $\it{Internal}$", transform=ax.transAxes, fontsize=17)
    ax.legend(loc="lower right", frameon=True, framealpha=0.95, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "cls_limit_real_data_dnn_score.png"), dpi=150)
    fig.savefig(os.path.join(PLOT_DIR, "cls_limit_real_data_dnn_score.pdf"))
    plt.close(fig)

    with open(os.path.join(TABLE_DIR, "cls_limit_real_data_dnn_score.json"), "w", encoding="utf-8") as handle:
        json.dump(values, handle, indent=2, sort_keys=True)


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


def draw_scan():
    scan_path = os.path.join(ROOT_DIR, "higgsCombine.RealData_dnn_score_grid.MultiDimFit.mH125.root")
    x, y = read_scan(scan_path)
    best = float(x[np.argmin(y)])
    one_left = crossing(x, y, 1.0, "left", reference=best)
    one_right = crossing(x, y, 1.0, "right", reference=best)
    ul95 = crossing(x, y, 3.84, "right", reference=best)
    ul95_physical = crossing(x, y, 3.84, "right", reference=max(best, 0.0))
    singles = parse_singles(os.path.join(BASE, "logs", "dnn_score_singles.log"))
    if singles is not None:
        best, err_lo, err_hi = singles
    else:
        err_lo = one_left - best
        err_hi = one_right - best
    if abs(best) < 5e-4:
        best = 0.0

    fig, ax = plt.subplots(figsize=(9.5, 7.0))
    label = rf"Real data: {best:.3f}$^{{+{abs(err_hi):.3f}}}_{{-{abs(err_lo):.3f}}}$, 95% crossing={ul95_physical:.3f}"
    ax.plot(x, y, color="#ff7f0e", lw=2.2, marker="o", ms=3.0, alpha=0.95, label=label)
    if np.isfinite(ul95):
        ax.plot([ul95, ul95], [0.0, 3.84], color="#ff7f0e", ls="--", lw=2.0, alpha=0.45)
        ax.plot([ul95], [3.84], color="#ff7f0e", marker="o", ms=9)
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
    fig.savefig(os.path.join(PLOT_DIR, "likelihood_scan_real_data_dnn_score.png"), dpi=150)
    fig.savefig(os.path.join(PLOT_DIR, "likelihood_scan_real_data_dnn_score.pdf"))
    plt.close(fig)

    with open(os.path.join(TABLE_DIR, "likelihood_scan_real_data_dnn_score_summary.csv"), "w", newline="") as handle:
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
            "name": "RealData_dnn_score",
            "best": best,
            "err_lo": err_lo,
            "err_hi": err_hi,
            "ul95_unconstrained": ul95,
            "ul95_physical_positive": ul95_physical,
            "one_left": one_left,
            "one_right": one_right,
        })


def load_impact(param):
    path = os.path.join(ROOT_DIR, f"higgsCombine_paramFit_.real_data_dnn_score_impact_{param}.MultiDimFit.mH125.root")
    arr = uproot.open(path)["limit"].arrays([param, "r"], library="np")
    p = [float(x) for x in arr[param]]
    r = [float(x) for x in arr["r"]]
    p_fit = [p[1], p[0], p[2]]
    r_fit = [r[1], r[0], r[2]]
    impact = max(abs(r_fit[0] - r_fit[1]), abs(r_fit[2] - r_fit[1]))
    return {"name": param, "fit": p_fit, "prefit": [1.0, 1.0, 1.0], "r": r_fit, "impact_r": impact, "type": "Unconstrained"}


def format_fit_value(low, central, high):
    return rf"{central:.2f}$^{{+{high-central:.2f}}}_{{-{central-low:.2f}}}$"


def draw_ranking():
    params = [load_impact(param) for param in PARAMS]
    params.sort(key=lambda x: x["impact_r"], reverse=True)
    with open(os.path.join(JSON_DIR, "impacts_real_data_dnn_score.json"), "w", encoding="utf-8") as handle:
        json.dump({"params": params, "method": "manual MultiDimFit impact scan on real data_obs"}, handle, indent=2, sort_keys=True)

    names = [p["name"] for p in params]
    centers = np.array([p["fit"][1] for p in params])
    lows = np.array([p["fit"][0] for p in params])
    highs = np.array([p["fit"][2] for p in params])
    impact_minus = np.array([p["r"][0] - p["r"][1] for p in params])
    impact_plus = np.array([p["r"][2] - p["r"][1] for p in params])
    impacts = np.maximum(np.abs(impact_minus), np.abs(impact_plus))
    y = np.arange(len(params))

    fig = plt.figure(figsize=(9.8, 4.4))
    grid = fig.add_gridspec(1, 2, left=0.34, right=0.96, top=0.78, bottom=0.25, width_ratios=[0.58, 0.42], wspace=0.0)
    ax_pull = fig.add_subplot(grid[0, 0])
    ax_impact = fig.add_subplot(grid[0, 1], sharey=ax_pull)
    for axis in (ax_pull, ax_impact):
        for yy in y:
            if yy % 2 == 0:
                axis.axhspan(yy - 0.5, yy + 0.5, color="#e3e3e3", zorder=0)
        axis.set_ylim(len(params) - 0.5, -0.5)

    fit_min = min(float(np.min(lows)), 1.0)
    fit_max = max(float(np.max(highs)), 1.0)
    fit_span = max(fit_max - fit_min, 0.1)
    ax_pull.set_xlim(fit_min - 0.18 * fit_span, fit_max + 0.25 * fit_span)
    ax_pull.axvline(1.0, color="black", lw=1.0)
    ax_pull.grid(axis="x", color="black", ls=":", alpha=0.65)
    ax_pull.errorbar(centers, y, xerr=[centers - lows, highs - centers], fmt="o", color="black", ecolor="black", capsize=3, ms=4, lw=1.2, zorder=4)
    ax_pull.scatter(np.ones_like(y), y, marker="x", color="blue", s=24, zorder=5)
    for idx, (name, low, central, high) in enumerate(zip(names, lows, centers, highs)):
        ax_pull.text(-0.55, idx, str(idx + 1), transform=ax_pull.get_yaxis_transform(), ha="center", va="center", fontsize=10, fontweight="bold", clip_on=False)
        ax_pull.text(-0.02, idx, name, transform=ax_pull.get_yaxis_transform(), ha="right", va="center", fontsize=9, color="#53658c", clip_on=False)
        ax_pull.text(central, idx - 0.25, format_fit_value(low, central, high), ha="center", va="center", fontsize=8.5)
    ax_pull.set_yticks([])
    ax_pull.set_xlabel("k-factor value", fontsize=12)
    ax_pull.tick_params(axis="x", direction="in", top=True, labelsize=11)

    impact_max = max(float(np.max(np.abs(impact_minus))), float(np.max(np.abs(impact_plus))), 1e-3)
    ax_impact.set_xlim(-1.28 * impact_max, 1.28 * impact_max)
    ax_impact.axvline(0.0, color="black", lw=1.0)
    ax_impact.grid(axis="x", color="black", ls=":", alpha=0.65)
    ax_impact.barh(y, impact_plus, color="#d98b91", height=0.82, label=r"+1$\sigma$ Impact")
    ax_impact.barh(y, impact_minus, color="#9db5d9", height=0.82, label=r"-1$\sigma$ Impact")
    for yy, val in zip(y, impacts):
        ax_impact.text(impact_max * 1.05, yy, f"{val:.3f}", va="center", fontsize=9)
    ax_impact.set_yticks(y)
    ax_impact.set_yticklabels([])
    ax_impact.tick_params(axis="x", direction="in", top=True, labelsize=11)
    ax_impact.set_xlabel(r"$\Delta r$", fontsize=15, loc="right")

    rhat, rdown, rup = parse_singles(os.path.join(BASE, "logs", "dnn_score_impacts_initial.log"))
    fig.text(0.34, 0.84, r"$\bf{CMS}$ $\it{Internal}$", fontsize=18, ha="left")
    fig.text(0.62, 0.84, "Real data DNN score", fontsize=13, ha="left")
    fig.text(0.96, 0.84, rf"$\hat{{r}} = {rhat:.2f}^{{+{abs(rup):.2f}}}_{{-{abs(rdown):.2f}}}$", fontsize=14, ha="right")
    fit_handles = [
        Line2D([0], [0], color="black", marker="o", lw=1.2, label="Fit"),
        Line2D([0], [0], color="blue", marker="x", lw=0, label="Prefit"),
        Patch(facecolor="#d98b91", edgecolor="black", label=r"+1$\sigma$ Impact"),
        Patch(facecolor="#9db5d9", edgecolor="black", label=r"-1$\sigma$ Impact"),
    ]
    fig.legend(handles=[Patch(facecolor="#8f93ad", edgecolor="black", label="Unconstrained")], loc="upper left", bbox_to_anchor=(0.12, 0.91), frameon=False, fontsize=10)
    fig.legend(handles=fit_handles, loc="lower left", bbox_to_anchor=(0.15, 0.03), frameon=False, ncol=2, fontsize=11)
    fig.savefig(os.path.join(PLOT_DIR, "ranking_real_data_dnn_score.png"), dpi=160)
    fig.savefig(os.path.join(PLOT_DIR, "ranking_real_data_dnn_score.pdf"))
    plt.close(fig)


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(TABLE_DIR, exist_ok=True)
    os.makedirs(JSON_DIR, exist_ok=True)
    draw_limit(load_limit())
    draw_scan()
    draw_ranking()


if __name__ == "__main__":
    main()
