#!/usr/bin/env python3
import csv
import json
import math
import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


BASE = os.environ.get(
    "MONOZ_CLS_BASE",
    "/afs/cern.ch/user/l/liwe/hzz2l2nu/combine/final/cls_limit_ranking",
)
ROOT_DIR = os.path.join(BASE, "root")
PLOT_DIR = os.path.join(BASE, "plots")
TABLE_DIR = os.path.join(BASE, "tables")
JSON_DIR = os.path.join(BASE, "json")

LABELS = {
    "counting": "Expected_counting",
    "shape": "Expected_shape",
}

COLORS = {
    "counting": "#1f77b4",
    "shape": "#ff7f0e",
}

PARAMS = ["k_Zjet", "k_WZ", "k_emu"]
ACTIVE_LABELS = (
    ["shape"]
    if os.environ.get("MONOZ_ONLY_SHAPE") == "1"
    else ["counting", "shape"]
)


def load_limit(label):
    path = os.path.join(ROOT_DIR, f"higgsCombine.Expected_{label}.AsymptoticLimits.mH125.root")
    tree = uproot.open(path)["limit"]
    arr = tree.arrays(["limit", "quantileExpected"], library="np")
    result = {}
    for limit, quantile in zip(arr["limit"], arr["quantileExpected"]):
        q = float(quantile)
        value = float(limit)
        if q < 0:
            result["observed_asimov"] = value
        elif abs(q - 0.025) < 1e-4:
            result["exp_m2"] = value
        elif abs(q - 0.16) < 1e-4:
            result["exp_m1"] = value
        elif abs(q - 0.5) < 1e-4:
            result["exp"] = value
        elif abs(q - 0.84) < 1e-4:
            result["exp_p1"] = value
        elif abs(q - 0.975) < 1e-4:
            result["exp_p2"] = value
    return result


def draw_limit_plot(limit_results):
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ypos = {
        label: float(len(ACTIVE_LABELS) - index - 1)
        for index, label in enumerate(ACTIVE_LABELS)
    }

    for label in ACTIVE_LABELS:
        values = limit_results[label]
        y = ypos[label]
        ax.barh(
            y,
            values["exp_p2"] - values["exp_m2"],
            left=values["exp_m2"],
            height=0.46,
            color="#f4d03f",
            edgecolor="none",
            label=r"Expected $\pm2\sigma$" if label == ACTIVE_LABELS[0] else None,
        )
        ax.barh(
            y,
            values["exp_p1"] - values["exp_m1"],
            left=values["exp_m1"],
            height=0.30,
            color="#58d68d",
            edgecolor="none",
            label=r"Expected $\pm1\sigma$" if label == ACTIVE_LABELS[0] else None,
        )
        ax.plot(
            [values["exp"], values["exp"]],
            [y - 0.24, y + 0.24],
            color="black",
            lw=2.0,
            label="Expected median" if label == ACTIVE_LABELS[0] else None,
        )
        ax.plot(
            values["observed_asimov"],
            y,
            marker="o",
            ms=7,
            color=COLORS[label],
            label="Asimov observed" if label == ACTIVE_LABELS[0] else None,
        )
        ax.text(
            values["exp_p2"] + 0.035,
            y,
            f"{values['exp']:.3f}",
            va="center",
            ha="left",
            fontsize=12,
        )

    ax.set_yticks([ypos[label] for label in ACTIVE_LABELS])
    ax.set_yticklabels([LABELS[label] for label in ACTIVE_LABELS], fontsize=13)
    ax.set_xlabel(r"95% CL upper limit on $BR_{inv}$", fontsize=15)
    ax.set_xlim(0.0, max(v["exp_p2"] for v in limit_results.values()) * 1.25)
    ax.set_ylim(-0.6, max(ypos.values()) + 0.6)
    ax.grid(axis="x", ls="--", alpha=0.3)
    ax.text(0.0, 1.02, r"$\bf{CMS}$ $\it{Internal}$", transform=ax.transAxes, fontsize=17)
    ax.legend(loc="lower right", frameon=True, framealpha=0.95, fontsize=10)
    fig.tight_layout()
    stem = (
        "cls_limit_expected_shape"
        if ACTIVE_LABELS == ["shape"]
        else "cls_limit_expected_counting_shape"
    )
    fig.savefig(os.path.join(PLOT_DIR, f"{stem}.png"), dpi=140)
    fig.savefig(os.path.join(PLOT_DIR, f"{stem}.pdf"))
    plt.close(fig)


def load_impact(label, param):
    path = os.path.join(
        ROOT_DIR, f"higgsCombine_paramFit_.{label}_impact_{param}.MultiDimFit.mH125.root"
    )
    tree = uproot.open(path)["limit"]
    branches = tree.arrays([param, "r"], library="np")
    p = [float(x) for x in branches[param]]
    r = [float(x) for x in branches["r"]]
    if len(p) != 3 or len(r) != 3:
        raise RuntimeError(f"Unexpected impact entries in {path}")
    # Combine writes entries as central, lower-param fit, upper-param fit.
    p_fit = [p[1], p[0], p[2]]
    r_fit = [r[1], r[0], r[2]]
    impact = max(abs(r_fit[0] - r_fit[1]), abs(r_fit[2] - r_fit[1]))
    return {
        "name": param,
        "fit": p_fit,
        "prefit": [1.0, 1.0, 1.0],
        "r": r_fit,
        "impact_r": impact,
        "type": "Unconstrained",
        "groups": [],
    }


def make_impact_json(label):
    params = [load_impact(label, param) for param in PARAMS]
    params.sort(key=lambda x: x["impact_r"], reverse=True)
    poi_vals = []
    for param in params:
        poi_vals.extend(param["r"])
    poi_central = 0.0
    data = {
        "POIs": [{"name": "r", "fit": [min(poi_vals), poi_central, max(poi_vals)]}],
        "method": "manual MultiDimFit impact scan on asimovData_1",
        "params": params,
    }
    path = os.path.join(JSON_DIR, f"impacts_{label}.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
    return data


def read_poi_fit(label):
    path = os.path.join(BASE, "logs", f"{label}_impacts_initial.log")
    text = open(path, encoding="utf-8").read()
    match = re.search(
        r"\br\s*:\s*([+-]?\d+(?:\.\d+)?)\s*"
        r"([+-]\d+(?:\.\d+)?)/([+-]\d+(?:\.\d+)?)",
        text,
    )
    if not match:
        return 0.0, 0.0, 0.0
    central = float(match.group(1))
    down = abs(float(match.group(2)))
    up = abs(float(match.group(3)))
    if abs(central) < 5e-4:
        central = 0.0
    return central, down, up


def format_fit_value(low, central, high):
    up = high - central
    down = central - low
    return rf"{central:.2f}$^{{+{up:.2f}}}_{{-{down:.2f}}}$"


def draw_single_ranking(label, data):
    params = data["params"]
    names = [p["name"] for p in params]
    centers = np.array([p["fit"][1] for p in params])
    lows = np.array([p["fit"][0] for p in params])
    highs = np.array([p["fit"][2] for p in params])
    impact_minus = np.array([p["r"][0] - p["r"][1] for p in params])
    impact_plus = np.array([p["r"][2] - p["r"][1] for p in params])
    impacts = np.maximum(np.abs(impact_minus), np.abs(impact_plus))
    y = np.arange(len(params))

    fig = plt.figure(figsize=(9.8, 4.4))
    grid = fig.add_gridspec(
        1,
        2,
        left=0.34,
        right=0.96,
        top=0.78,
        bottom=0.25,
        width_ratios=[0.58, 0.42],
        wspace=0.0,
    )
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
    ax_pull.errorbar(
        centers,
        y,
        xerr=[centers - lows, highs - centers],
        fmt="o",
        color="black",
        ecolor="black",
        capsize=3,
        ms=4,
        lw=1.2,
        zorder=4,
    )
    ax_pull.scatter(np.ones_like(y), y, marker="x", color="blue", s=24, zorder=5)
    for idx, (name, low, central, high) in enumerate(zip(names, lows, centers, highs)):
        ax_pull.text(
            -0.55,
            idx,
            str(idx + 1),
            transform=ax_pull.get_yaxis_transform(),
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            clip_on=False,
        )
        ax_pull.text(
            -0.02,
            idx,
            name,
            transform=ax_pull.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=9,
            color="#53658c" if name.startswith("k_") else "black",
            clip_on=False,
        )
        ax_pull.text(
            central,
            idx - 0.25,
            format_fit_value(low, central, high),
            ha="center",
            va="center",
            fontsize=8.5,
        )
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

    rhat, rdown, rup = read_poi_fit(label)
    fig.text(0.34, 0.84, r"$\bf{CMS}$ $\it{Internal}$", fontsize=18, ha="left")
    fig.text(0.62, 0.84, LABELS[label], fontsize=13, ha="left")
    fig.text(
        0.96,
        0.84,
        rf"$\hat{{r}} = {rhat:.2f}^{{+{rup:.2f}}}_{{-{rdown:.2f}}}$",
        fontsize=14,
        ha="right",
    )

    type_handles = [Patch(facecolor="#8f93ad", edgecolor="black", label="Unconstrained")]
    fit_handles = [
        Line2D([0], [0], color="black", marker="o", lw=1.2, label="Fit"),
        Line2D([0], [0], color="blue", marker="x", lw=0, label="Prefit"),
        Patch(facecolor="#d98b91", edgecolor="black", label=r"+1$\sigma$ Impact"),
        Patch(facecolor="#9db5d9", edgecolor="black", label=r"-1$\sigma$ Impact"),
    ]
    fig.legend(handles=type_handles, loc="upper left", bbox_to_anchor=(0.12, 0.91), frameon=False, fontsize=10)
    fig.legend(handles=fit_handles, loc="lower left", bbox_to_anchor=(0.15, 0.03), frameon=False, ncol=2, fontsize=11)

    fig.savefig(os.path.join(PLOT_DIR, f"ranking_{label}.png"), dpi=160)
    fig.savefig(os.path.join(PLOT_DIR, f"ranking_{label}.pdf"))
    plt.close(fig)


def draw_combined_ranking(all_data):
    images = [
        plt.imread(os.path.join(PLOT_DIR, "ranking_counting.png")),
        plt.imread(os.path.join(PLOT_DIR, "ranking_shape.png")),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(9.8, 8.8))
    for ax, image in zip(axes, images):
        ax.imshow(image)
        ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0.02)
    fig.savefig(os.path.join(PLOT_DIR, "ranking_expected_counting_shape.png"), dpi=150)
    fig.savefig(os.path.join(PLOT_DIR, "ranking_expected_counting_shape.pdf"))
    plt.close(fig)


def write_limit_table(limit_results):
    path = os.path.join(TABLE_DIR, "cls_limit_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "label",
                "observed_asimov",
                "expected_median",
                "expected_m2",
                "expected_m1",
                "expected_p1",
                "expected_p2",
            ],
        )
        writer.writeheader()
        for label in ACTIVE_LABELS:
            v = limit_results[label]
            writer.writerow(
                {
                    "label": LABELS[label],
                    "observed_asimov": v["observed_asimov"],
                    "expected_median": v["exp"],
                    "expected_m2": v["exp_m2"],
                    "expected_m1": v["exp_m1"],
                    "expected_p1": v["exp_p1"],
                    "expected_p2": v["exp_p2"],
                }
            )


def write_impact_table(all_data):
    path = os.path.join(TABLE_DIR, "ranking_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "label",
                "parameter",
                "postfit_low",
                "postfit_central",
                "postfit_high",
                "r_at_low",
                "r_central",
                "r_at_high",
                "impact",
            ],
        )
        writer.writeheader()
        for label in ACTIVE_LABELS:
            for param in all_data[label]["params"]:
                writer.writerow(
                    {
                        "label": LABELS[label],
                        "parameter": param["name"],
                        "postfit_low": param["fit"][0],
                        "postfit_central": param["fit"][1],
                        "postfit_high": param["fit"][2],
                        "r_at_low": param["r"][0],
                        "r_central": param["r"][1],
                        "r_at_high": param["r"][2],
                        "impact": param["impact_r"],
                    }
                )


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(TABLE_DIR, exist_ok=True)
    os.makedirs(JSON_DIR, exist_ok=True)

    limit_results = {label: load_limit(label) for label in ACTIVE_LABELS}
    draw_limit_plot(limit_results)
    write_limit_table(limit_results)

    all_data = {label: make_impact_json(label) for label in ACTIVE_LABELS}
    for label in ACTIVE_LABELS:
        draw_single_ranking(label, all_data[label])
    if len(ACTIVE_LABELS) > 1:
        draw_combined_ranking(all_data)
    write_impact_table(all_data)


if __name__ == "__main__":
    main()
