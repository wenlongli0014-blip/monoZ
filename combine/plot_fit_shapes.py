#!/usr/bin/env python3

import argparse
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot


DEFAULT_TAGS = ["ptmiss_repro", "mt", "ptsigcorr"]
FIT_DIRS = {
    "prefit": "shapes_prefit",
    "postfit_b": "shapes_fit_b",
    "postfit_s": "shapes_fit_s",
}
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


def values_and_errors(directory, name):
    obj = directory.get(name)
    if obj is None:
        return None, None
    return np.asarray(obj.values(), dtype=float), np.asarray(obj.errors(), dtype=float)


def graph_values(graph):
    x, y = graph.values()
    _, yerr_low = graph.errors("low")
    _, yerr_high = graph.errors("high")
    return (
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(yerr_low, dtype=float),
        np.asarray(yerr_high, dtype=float),
    )


def active_nbins(directory):
    total, total_err = values_and_errors(directory, "total_background")
    data_x, data_y, _, _ = graph_values(directory["data"])
    n_bins = len(total)
    active = (total > 0) | (total_err > 0) | (data_y[:n_bins] > 0)

    for proc in BACKGROUND_ORDER + ["signal"]:
        vals, errs = values_and_errors(directory, proc)
        if vals is not None:
            active |= (vals > 0) | (errs > 0)

    nonzero = np.where(active)[0]
    return int(nonzero[-1] + 1) if len(nonzero) else 1


def draw_one(root_file, tag, fit_label, fit_dir, channel, outdir):
    path = f"{fit_dir}/{channel}"
    if path not in root_file:
        return

    directory = root_file[path]
    n_bins = active_nbins(directory)
    bins = np.arange(n_bins)
    centers = bins + 0.5

    total_bkg, total_bkg_err = values_and_errors(directory, "total_background")
    total_bkg = total_bkg[:n_bins]
    total_bkg_err = total_bkg_err[:n_bins]

    data_x, data_y, data_ey_low, data_ey_high = graph_values(directory["data"])
    data_y = data_y[:n_bins]
    data_ey_low = data_ey_low[:n_bins]
    data_ey_high = data_ey_high[:n_bins]

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(8, 8),
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.05},
        sharex=True,
    )

    bottom = np.zeros(n_bins)
    for proc in BACKGROUND_ORDER:
        vals, _ = values_and_errors(directory, proc)
        if vals is None:
            continue
        vals = vals[:n_bins]
        ax.bar(
            centers,
            vals,
            bottom=bottom,
            width=1.0,
            color=COLORS.get(proc, "#999999"),
            edgecolor="black",
            linewidth=0.4,
            label=proc,
        )
        bottom += vals

    ax.fill_between(
        bins.tolist() + [n_bins],
        np.r_[total_bkg - total_bkg_err, total_bkg[-1] - total_bkg_err[-1]],
        np.r_[total_bkg + total_bkg_err, total_bkg[-1] + total_bkg_err[-1]],
        step="post",
        color="0.35",
        alpha=0.30,
        label="Bkg. unc.",
    )

    sig, _ = values_and_errors(directory, "signal")
    if sig is not None and np.any(sig[:n_bins] > 0):
        sig = sig[:n_bins]
        ax.step(np.r_[bins, n_bins], np.r_[sig, sig[-1]], where="post", color="crimson", linewidth=2.0, label="signal")

    ax.errorbar(
        centers,
        data_y,
        yerr=[data_ey_low, data_ey_high],
        fmt="o",
        color="black",
        markersize=4,
        linewidth=1,
        label="Data",
    )

    ymax = max(float(np.max(total_bkg + total_bkg_err)), float(np.max(data_y + data_ey_high)), 1.0)
    if sig is not None:
        ymax = max(ymax, float(np.max(sig)))
    ax.set_ylim(0, ymax * 1.55)
    ax.set_ylabel("Events")
    ax.text(0.02, 0.95, f"{tag}  {channel}  {fit_label}", transform=ax.transAxes, va="top", fontsize=12)
    ax.text(0.02, 0.89, "Stat-only model", transform=ax.transAxes, va="top", fontsize=10)
    ax.legend(ncol=2, fontsize=8, frameon=False, loc="upper right")

    ratio = np.divide(data_y, total_bkg, out=np.zeros_like(data_y), where=total_bkg > 0)
    ratio_low = np.divide(data_ey_low, total_bkg, out=np.zeros_like(data_y), where=total_bkg > 0)
    ratio_high = np.divide(data_ey_high, total_bkg, out=np.zeros_like(data_y), where=total_bkg > 0)
    rel_unc = np.divide(total_bkg_err, total_bkg, out=np.zeros_like(total_bkg), where=total_bkg > 0)

    rax.fill_between(
        bins.tolist() + [n_bins],
        np.r_[1 - rel_unc, 1 - rel_unc[-1]],
        np.r_[1 + rel_unc, 1 + rel_unc[-1]],
        step="post",
        color="0.35",
        alpha=0.30,
    )
    rax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    rax.errorbar(centers, ratio, yerr=[ratio_low, ratio_high], fmt="o", color="black", markersize=4, linewidth=1)
    rax.set_ylim(0, 2)
    rax.set_ylabel("Data/Bkg")
    rax.set_xlabel("Template bin")
    rax.set_xlim(0, n_bins)
    rax.set_xticks(centers)
    rax.set_xticklabels([str(i) for i in range(1, n_bins + 1)])

    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, f"{tag}_{channel}_{fit_label}")
    fig.savefig(base + ".png", dpi=160, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Draw prefit/postfit plots from fitDiagnostics ROOT files.")
    parser.add_argument("--outdir", default="plots/fit_shapes", help="Output directory")
    parser.add_argument("--tags", nargs="*", default=DEFAULT_TAGS, help="Fit tags to plot")
    parser.add_argument("--channels", nargs="*", default=["SR", "DYCR", "EMUCR", "CR3L"], help="Channels to plot")
    args = parser.parse_args()

    for tag in args.tags:
        fit_path = f"fitDiagnostics.{tag}_fit.root"
        with uproot.open(fit_path) as root_file:
            for fit_label, fit_dir in FIT_DIRS.items():
                for channel in args.channels:
                    draw_one(root_file, tag, fit_label, fit_dir, channel, args.outdir)


if __name__ == "__main__":
    main()
