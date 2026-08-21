#!/usr/bin/env python3

import argparse
import json
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def clean_name(name):
    return name.replace("manual_", "").replace("prop_bin", "mcstat_")


def impact_components(param, poi_central):
    r_values = param.get("r", [])
    if len(r_values) != 3:
        impact = abs(float(param.get("impact_r", 0.0)))
        return impact, impact

    low, _, high = [float(x) for x in r_values]
    down = abs(poi_central - low)
    up = abs(high - poi_central)
    return down, up


def postfit_width(param):
    fit = param.get("fit", [])
    if len(fit) != 3:
        return 0.0
    return 0.5 * abs(float(fit[2]) - float(fit[0]))


def finite_param(param):
    values = []
    values.extend(param.get("fit", []))
    values.extend(param.get("r", []))
    values.append(param.get("impact_r", 0.0))
    return all(math.isfinite(float(value)) for value in values)


def draw_page(pdf, params, title, page_index, n_pages, poi_central):
    names = [clean_name(param["name"]) for param in params]
    centers = np.array([float(param["fit"][1]) for param in params])
    lows = np.array([float(param["fit"][0]) for param in params])
    highs = np.array([float(param["fit"][2]) for param in params])
    err_low = centers - lows
    err_high = highs - centers

    impact_down = []
    impact_up = []
    for param in params:
        down, up = impact_components(param, poi_central)
        impact_down.append(down)
        impact_up.append(up)
    impact_down = np.array(impact_down)
    impact_up = np.array(impact_up)

    y = np.arange(len(params))
    fig, (ax_pull, ax_impact) = plt.subplots(
        1,
        2,
        figsize=(14.0, max(5.0, 0.48 * len(params) + 1.6)),
        gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.18},
    )

    ax_pull.axvspan(-1, 1, color="#e9ecef", alpha=0.9, label="prefit +/-1")
    ax_pull.axvline(0, color="black", linewidth=1.0)
    ax_pull.errorbar(
        centers,
        y,
        xerr=[err_low, err_high],
        fmt="o",
        color="#1f4e79",
        ecolor="#1f4e79",
        capsize=3,
        markersize=4,
    )
    ax_pull.set_yticks(y)
    ax_pull.set_yticklabels(names, fontsize=8)
    ax_pull.invert_yaxis()
    ax_pull.grid(axis="x", linestyle=":", alpha=0.5)
    ax_pull.set_xlabel("postfit nuisance value")

    finite_x = np.concatenate([lows, highs, np.array([-1.0, 1.0])])
    xmin = float(np.min(finite_x))
    xmax = float(np.max(finite_x))
    span = max(xmax - xmin, 2.0)
    ax_pull.set_xlim(xmin - 0.08 * span, xmax + 0.08 * span)

    ax_impact.barh(y, impact_down, color="#ef8a62", alpha=0.85, label="down")
    ax_impact.barh(y, impact_up, left=impact_down, color="#67a9cf", alpha=0.85, label="up")
    impact_total = impact_down + impact_up
    for y_value, value in zip(y, impact_total):
        label = f"{value:.2e}" if value < 0.01 else f"{value:.3f}"
        ax_impact.text(
            value,
            y_value,
            "  " + label,
            va="center",
            ha="left",
            fontsize=7,
            color="black",
        )
    ax_impact.set_yticks(y)
    ax_impact.set_yticklabels([])
    ax_impact.invert_yaxis()
    ax_impact.grid(axis="x", linestyle=":", alpha=0.5)
    ax_impact.set_xlabel("impact on r")
    ax_impact.legend(frameon=False, fontsize=8, loc="lower right")

    max_impact = float(np.max(impact_down + impact_up)) if len(params) else 0.0
    ax_impact.set_xlim(0, max(max_impact * 1.45, 1e-3))

    fig.suptitle(f"{title} ranking, page {page_index + 1}/{n_pages}", fontsize=12, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Draw robust ranking plots from Combine impacts JSON.")
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--per-page", type=int, default=14)
    parser.add_argument("--sort", choices=["impact", "pull", "constraint"], default="impact")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.input, encoding="utf-8") as handle:
        data = json.load(handle)

    poi_fit = data.get("POIs", [{}])[0].get("fit", [0.0, 0.0, 0.0])
    poi_central = float(poi_fit[1]) if len(poi_fit) == 3 else 0.0

    params = [param for param in data["params"] if finite_param(param)]
    if args.sort == "impact":
        params.sort(key=lambda param: abs(float(param.get("impact_r", 0.0))), reverse=True)
    elif args.sort == "pull":
        params.sort(key=lambda param: abs(float(param["fit"][1])), reverse=True)
    else:
        params.sort(key=postfit_width)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    n_pages = max(1, math.ceil(len(params) / args.per_page))
    with PdfPages(args.output) as pdf:
        for i_page in range(n_pages):
            page_params = params[i_page * args.per_page : (i_page + 1) * args.per_page]
            draw_page(pdf, page_params, args.title, i_page, n_pages, poi_central)

    print(f"Wrote {args.output} with {len(params)} parameters over {n_pages} pages")


if __name__ == "__main__":
    main()
