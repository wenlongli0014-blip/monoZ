#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot

from dnn_common import load_metadata


BASE = Path("/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn")
OUTDIR = BASE / "plots"


def save(fig, name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        fig.savefig(OUTDIR / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


def plot_training_history():
    hist = np.genfromtxt(BASE / "output/training_history.csv", delimiter=",", names=True)
    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax1.plot(hist["epoch"], hist["train_loss"], label="train loss", color="tab:blue")
    ax1.plot(hist["epoch"], hist["val_loss"], label="validation loss", color="tab:cyan")
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Weighted cross entropy loss")
    ax2 = ax1.twinx()
    ax2.plot(hist["epoch"], hist["val_auc"], label="validation AUC", color="tab:red")
    ax2.set_ylabel("Validation AUC")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="center right")
    ax1.grid(True, alpha=0.25)
    save(fig, "training_history")


def plot_test_roc():
    roc = np.genfromtxt(BASE / "output/test_roc.csv", delimiter=",", names=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax in axes:
        ax.plot(roc["fpr"], roc["tpr"], label="DNN")
        ax.plot([0, 1], [0, 1], "--", color="0.5", label="random")
        ax.set_xlabel("Background efficiency")
        ax.set_ylabel("Signal efficiency")
        ax.grid(True, alpha=0.25)
        ax.legend()
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[1].set_xlim(1e-4, 1)
    axes[1].set_ylim(0, 1)
    axes[1].set_xscale("log")
    save(fig, "test_roc")


def plot_test_score():
    data = np.load(BASE / "output/test_scores.npz")
    score = data["score"]
    label = data["label"]
    weight = data["train_weight"]

    bins = np.linspace(0, 1, 41)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(score[label == 1], bins=bins, weights=weight[label == 1], density=True, histtype="step", linewidth=1.8, label="signal")
    ax.hist(score[label == 0], bins=bins, weights=weight[label == 0], density=True, histtype="step", linewidth=1.8, label="background")
    ax.set_xlabel("DNN score")
    ax.set_ylabel("Normalized events")
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "test_score_distribution")


def plot_input_distributions():
    data = np.load(BASE / "output/test_scores.npz")
    metadata, _, _ = load_metadata(BASE / "output/dnn_metadata.json")
    feature_names = metadata["feature_names"]
    x = data["raw_features"]
    label = data["label"]
    weight = data["train_weight"]

    preferred = [
        "ptmiss",
        "mT",
        "ptmiss_significance_corrected",
        "ll_pt",
        "ll_mass",
        "dphi_visibles_ptmiss",
        "jet_size",
        "dijet_mass",
        "jet_pt_0",
        "jet_pt_1",
        "num_pv_good",
        "ptmiss_phi",
    ]
    indices = [feature_names.index(name) for name in preferred if name in feature_names]
    indices = indices[:12]

    ncols = 4
    nrows = int(np.ceil(len(indices) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.5 * nrows))
    axes = np.asarray(axes).ravel()
    for ax in axes[len(indices):]:
        ax.axis("off")
    for ax, idx in zip(axes, indices):
        name = feature_names[idx]
        lo, hi = np.nanpercentile(x[:, idx], [1, 99])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = float(np.nanmin(x[:, idx])), float(np.nanmax(x[:, idx]) + 1.0)
        bins = np.linspace(lo, hi, 45)
        ax.hist(x[label == 1, idx], bins=bins, weights=weight[label == 1], density=True, histtype="step", linewidth=1.6, label="signal")
        ax.hist(x[label == 0, idx], bins=bins, weights=weight[label == 0], density=True, histtype="step", linewidth=1.6, label="background")
        ax.set_xlabel(name)
        ax.set_ylabel("Normalized events")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    save(fig, "test_input_distributions")


def hist_values(root_file, channel, process):
    with uproot.open(root_file) as f:
        hist = f[f"{channel}/{process}"]
        values, edges = hist.to_numpy(flow=False)
    return values, edges


def plot_sr_templates():
    root_file = BASE / "shapes_dnn_score.root"
    signal, edges = hist_values(root_file, "SR", "signal")
    bkg_total = np.zeros_like(signal)
    for process in ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]:
        vals, _ = hist_values(root_file, "SR", process)
        bkg_total += vals

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(centers, bkg_total, width=widths, align="center", alpha=0.45, label="total background")
    ax.step(edges[:-1], signal, where="post", linewidth=1.8, label="signal")
    ax.set_xlabel("DNN score")
    ax.set_ylabel("Expected yield")
    ax.set_xlim(0, 1)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save(fig, "sr_dnn_template_signal_vs_background")

    ratio = np.divide(signal, bkg_total, out=np.zeros_like(signal), where=bkg_total > 0)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(centers, ratio, width=widths, align="center", alpha=0.7)
    ax.set_xlabel("DNN score")
    ax.set_ylabel("S / B")
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.25)
    save(fig, "sr_dnn_template_s_over_b")


def plot_channel_data_mc():
    root_file = BASE / "shapes_dnn_score.root"
    processes = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]
    channels = ["SR", "DYCR", "EMUCR", "CR3L"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, channel in zip(axes.ravel(), channels):
        data, edges = hist_values(root_file, channel, "data_obs")
        bkg = np.zeros_like(data)
        for process in processes:
            vals, _ = hist_values(root_file, channel, process)
            bkg += vals
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = np.diff(edges)
        ax.bar(centers, bkg, width=widths, align="center", alpha=0.45, label="MC bkg")
        ax.errorbar(centers, data, yerr=np.sqrt(np.maximum(data, 1.0)), fmt="o", color="black", label="data")
        ax.set_title(channel)
        ax.set_xlabel("DNN score")
        ax.set_ylabel("Events")
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    save(fig, "channel_dnn_data_mc")


def main():
    plot_training_history()
    plot_test_roc()
    plot_test_score()
    plot_input_distributions()
    plot_sr_templates()
    plot_channel_data_mc()
    print(f"Wrote plots to {OUTDIR}")


if __name__ == "__main__":
    main()
