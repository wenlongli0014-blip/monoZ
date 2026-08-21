#!/usr/bin/env python3

import argparse
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.backends.backend_pdf import PdfPages


BACKGROUND_ORDER = ["DY", "WZ", "qqZZ", "ggZZ", "ttbar", "WW", "ST", "VVV", "Other"]
COLORS = {
    "DY": "#4C78A8",
    "WZ": "#F58518",
    "qqZZ": "#54A24B",
    "ggZZ": "#B279A2",
    "ttbar": "#E45756",
    "WW": "#72B7B2",
    "ST": "#FF9DA6",
    "VVV": "#9D755D",
    "Other": "#BAB0AC",
}
FIT_DIRS = {"prefit": "shapes_prefit", "postfit": "shapes_fit_s"}


def values_and_errors(directory, name):
    obj = directory.get(name)
    if obj is None:
        return None, None
    return np.asarray(obj.values(), dtype=float), np.asarray(obj.errors(), dtype=float)


def graph_values(graph):
    x_values, y_values = graph.values()
    _, yerr_low = graph.errors("low")
    _, yerr_high = graph.errors("high")
    return (
        np.asarray(x_values, dtype=float),
        np.asarray(y_values, dtype=float),
        np.asarray(yerr_low, dtype=float),
        np.asarray(yerr_high, dtype=float),
    )


def format_edge(value):
    if abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    return f"{value:g}"


def bin_labels_from_edges(edges):
    labels = []
    for i_edge, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
        if i_edge == len(edges) - 2:
            labels.append(f">={format_edge(low)}")
        else:
            labels.append(f"{format_edge(low)}-{format_edge(high)}")
    return labels


def draw_channel(root_file, shapes_file, fit_dir, tag, fit_label, channel,
                 observable, subtitle, is_asimov=False):
    directory = root_file[f"{fit_dir}/{channel}"]
    total_bkg, total_bkg_err = values_and_errors(directory, "total_background")
    edges = np.asarray(shapes_file[f"{channel}/data_obs"].axes[0].edges(), dtype=float)
    labels = bin_labels_from_edges(edges)
    n_bins = min(len(total_bkg), len(labels))

    total_bkg = total_bkg[:n_bins]
    total_bkg_err = total_bkg_err[:n_bins]
    labels = labels[:n_bins]
    bins = np.arange(n_bins)
    centers = bins + 0.5

    _, data_y, data_ey_low, data_ey_high = graph_values(directory["data"])
    data_y = data_y[:n_bins]
    data_ey_low = data_ey_low[:n_bins]
    data_ey_high = data_ey_high[:n_bins]

    fig, (axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(8.4, 8.0),
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.05},
        sharex=True,
    )

    bottom = np.zeros(n_bins)
    for process in BACKGROUND_ORDER:
        values, _ = values_and_errors(directory, process)
        if values is None:
            continue
        values = values[:n_bins]
        axis.bar(
            centers,
            values,
            bottom=bottom,
            width=1.0,
            color=COLORS.get(process, "#999999"),
            edgecolor="black",
            linewidth=0.35,
            label=process,
        )
        bottom += values

    axis.fill_between(
        np.r_[bins, n_bins],
        np.r_[total_bkg - total_bkg_err, total_bkg[-1] - total_bkg_err[-1]],
        np.r_[total_bkg + total_bkg_err, total_bkg[-1] + total_bkg_err[-1]],
        step="post",
        color="0.35",
        alpha=0.30,
        label="Bkg. unc.",
    )

    signal, _ = values_and_errors(directory, "signal")
    if signal is not None and np.any(signal[:n_bins] > 0):
        signal = signal[:n_bins]
        axis.step(
            np.r_[bins, n_bins],
            np.r_[signal, signal[-1]],
            where="post",
            color="crimson",
            linewidth=2.0,
            label="signal",
        )

    axis.errorbar(
        centers,
        data_y,
        yerr=[data_ey_low, data_ey_high],
        fmt="o",
        color="black",
        markersize=4,
        linewidth=1,
        label="Pseudodata" if is_asimov else "Data",
    )

    ymax = max(float(np.max(total_bkg + total_bkg_err)), float(np.max(data_y + data_ey_high)), 1.0)
    if signal is not None:
        ymax = max(ymax, float(np.max(signal)))
    axis.set_ylim(0, ymax * 1.55)
    axis.set_ylabel("Events")
    axis.text(0.02, 0.95, f"{tag}  {channel}  {fit_label}", transform=axis.transAxes, va="top", fontsize=12)
    axis.text(0.02, 0.89, subtitle, transform=axis.transAxes, va="top", fontsize=10)
    axis.legend(ncol=2, fontsize=8, frameon=False, loc="upper right")

    ratio = np.divide(data_y, total_bkg, out=np.zeros_like(data_y), where=total_bkg > 0)
    ratio_low = np.divide(data_ey_low, total_bkg, out=np.zeros_like(data_y), where=total_bkg > 0)
    ratio_high = np.divide(data_ey_high, total_bkg, out=np.zeros_like(data_y), where=total_bkg > 0)
    rel_unc = np.divide(total_bkg_err, total_bkg, out=np.zeros_like(total_bkg), where=total_bkg > 0)

    ratio_axis.fill_between(
        np.r_[bins, n_bins],
        np.r_[1 - rel_unc, 1 - rel_unc[-1]],
        np.r_[1 + rel_unc, 1 + rel_unc[-1]],
        step="post",
        color="0.35",
        alpha=0.30,
    )
    ratio_axis.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ratio_axis.errorbar(centers, ratio, yerr=[ratio_low, ratio_high], fmt="o", color="black", markersize=4)
    ratio_axis.set_ylim(0, 2)
    ratio_axis.set_ylabel("Pseudo./Bkg" if is_asimov else "Data/Bkg")
    ratio_axis.set_xlabel(observable)
    ratio_axis.set_xlim(0, n_bins)
    ratio_axis.set_xticks(centers)
    rotation = 90 if n_bins > 10 else 35
    ratio_axis.set_xticklabels(labels, rotation=rotation, ha="right")
    return fig


def parse_args():
    parser = argparse.ArgumentParser(description="Draw multipage prefit and postfit PDFs from fitDiagnostics.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--fitdiagnostics", required=True)
    parser.add_argument("--shapes", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--channels", nargs="+", default=["SR", "DYCR", "EMUCR", "CR3L"])
    parser.add_argument("--observable", required=True)
    parser.add_argument("--subtitle", default="Stat-only model")
    parser.add_argument(
        "--asimov-channels",
        nargs="*",
        default=[],
        help="Channels whose data points are background-only pseudodata.",
    )
    parser.add_argument(
        "--png",
        action="store_true",
        help="Also write a separate PNG for each channel.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with uproot.open(args.fitdiagnostics) as root_file, uproot.open(args.shapes) as shapes_file:
        for fit_label, fit_dir in FIT_DIRS.items():
            out_pdf = os.path.join(args.outdir, f"{args.tag}_{fit_label}.pdf")
            with PdfPages(out_pdf) as pdf:
                for channel in args.channels:
                    fig = draw_channel(
                        root_file,
                        shapes_file,
                        fit_dir,
                        args.tag,
                        fit_label,
                        channel,
                        args.observable,
                        args.subtitle,
                        is_asimov=channel in args.asimov_channels,
                    )
                    pdf.savefig(fig, bbox_inches="tight")
                    if args.png:
                        png_path = os.path.join(
                            args.outdir,
                            f"{args.tag}_{fit_label}_{channel}.png",
                        )
                        fig.savefig(png_path, bbox_inches="tight", dpi=160)
                        print(f"Wrote {png_path}")
                    plt.close(fig)
            print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()
