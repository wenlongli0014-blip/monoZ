#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import uproot


BASE = Path("/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn")
OUT = BASE / "training_validation_plots"
YEARS = ["2017", "2018"]
BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]


def ensure_out() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, name: str) -> None:
    ensure_out()
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / f"{name}.{ext}", dpi=170 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def load_metadata(year: str) -> dict:
    with open(BASE / year / "output" / "dnn_metadata.json") as handle:
        return json.load(handle)


def load_history(year: str):
    return np.genfromtxt(BASE / year / "output" / "training_history.csv", delimiter=",", names=True)


def load_test_scores(year: str):
    return np.load(BASE / year / "output" / "test_scores.npz")


def plot_history() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colors = {"2017": "#1f77b4", "2018": "#d62728"}
    for year in YEARS:
        hist = load_history(year)
        metadata = load_metadata(year)
        summary = metadata["training_summary"]
        axes[0].plot(hist["epoch"], hist["train_loss"], color=colors[year], lw=2, label=f"{year} train")
        axes[0].plot(hist["epoch"], hist["val_loss"], color=colors[year], lw=2, ls="--", label=f"{year} validation")
        axes[1].plot(
            hist["epoch"],
            hist["val_auc"],
            color=colors[year],
            lw=2,
            label=f"{year}: best={summary['best_val_auc']:.3f}, test={summary['test_auc']:.3f}",
        )
        best_idx = int(np.argmax(hist["val_auc"]))
        axes[1].plot(hist["epoch"][best_idx], hist["val_auc"][best_idx], marker="o", color=colors[year], ms=5)

    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Weighted cross entropy loss")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=9)
    axes[0].set_title("Training and validation loss")

    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Validation AUC")
    axes[1].set_ylim(0.80, 0.88)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9)
    axes[1].set_title("Validation AUC")

    fig.suptitle("DNN training convergence", fontsize=15)
    save(fig, "training_loss_auc_2017_2018")


def plot_roc() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colors = {"2017": "#1f77b4", "2018": "#d62728"}
    for year in YEARS:
        roc = np.genfromtxt(BASE / year / "output" / "test_roc.csv", delimiter=",", names=True)
        test_auc = load_metadata(year)["training_summary"]["test_auc"]
        label = f"{year} test AUC={test_auc:.3f}"
        for ax in axes:
            ax.plot(roc["fpr"], roc["tpr"], lw=2, color=colors[year], label=label)

    for ax in axes:
        ax.plot([0, 1], [0, 1], color="0.55", ls="--", lw=1.3, label="random")
        ax.set_xlabel("Background efficiency")
        ax.set_ylabel("Signal efficiency")
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
    axes[0].set_xlim(0, 1)
    axes[0].set_title("Linear background efficiency")
    axes[1].set_xlim(1e-4, 1)
    axes[1].set_xscale("log")
    axes[1].set_title("Low background efficiency")
    fig.suptitle("DNN ROC curves on held-out test sample", fontsize=15)
    save(fig, "test_roc_2017_2018")


def plot_test_score() -> None:
    bins = np.linspace(0, 1, 51)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, year in zip(axes, YEARS):
        data = load_test_scores(year)
        score = data["score"]
        label = data["label"]
        weight = data["train_weight"]
        ax.hist(
            score[label == 0],
            bins=bins,
            weights=weight[label == 0],
            density=True,
            histtype="stepfilled",
            alpha=0.25,
            color="#4c78a8",
            label="background",
        )
        ax.hist(
            score[label == 1],
            bins=bins,
            weights=weight[label == 1],
            density=True,
            histtype="step",
            lw=2.0,
            color="#e45756",
            label="signal",
        )
        ax.set_xlabel("DNN score")
        ax.set_ylabel("Normalized weighted events")
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.set_title(year)
        ax.legend(fontsize=9)
    fig.suptitle("DNN score separation on held-out test sample", fontsize=15)
    save(fig, "test_score_signal_background_2017_2018")


def arrays_from_scored_root(path: Path, branches: list[str]) -> dict[str, np.ndarray]:
    with uproot.open(path)["Vars"] as tree:
        arrays = tree.arrays(branches, library="np")
    return {key: np.asarray(value) for key, value in arrays.items()}


def selected_scores(path: Path, weighted: bool = True) -> tuple[np.ndarray, np.ndarray]:
    branches = ["dnn_score", "lepton_cat", "ptmiss"]
    if weighted:
        branches.append("weight")
    arrays = arrays_from_scored_root(path, branches)
    mask = (
        (arrays["lepton_cat"] != 2)
        & (arrays["ptmiss"] >= 100.0)
        & np.isfinite(arrays["dnn_score"])
        & (arrays["dnn_score"] >= 0.0)
        & (arrays["dnn_score"] <= 1.0)
    )
    score = arrays["dnn_score"][mask].astype(float)
    if weighted:
        weight = arrays["weight"][mask].astype(float)
        weight = np.where(np.isfinite(weight), weight, 0.0)
    else:
        weight = np.ones_like(score)
    return score, weight


def scored_base(year: str) -> Path:
    return Path(f"/eos/user/l/liwe/monoZ_combine_17_18_dnn/{year}_dnn")


def plot_sr_score_templates() -> None:
    bins = np.linspace(0, 1, 41)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, year in zip(axes, YEARS):
        base = scored_base(year)
        sig_score, sig_w = selected_scores(base / "SR" / "signal.root", True)
        bkg_scores = []
        bkg_weights = []
        for process in BKG_PROCESSES:
            score, weight = selected_scores(base / "SR" / f"{process}.root", True)
            bkg_scores.append(score)
            bkg_weights.append(weight)
        bkg_score = np.concatenate(bkg_scores)
        bkg_w = np.concatenate(bkg_weights)

        sig_hist, edges = np.histogram(sig_score, bins=bins, weights=sig_w)
        bkg_hist, _ = np.histogram(bkg_score, bins=bins, weights=bkg_w)
        widths = np.diff(edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.bar(centers, np.maximum(bkg_hist, 1e-8), width=widths, color="#4c78a8", alpha=0.35, label="total background")
        ax.stairs(np.maximum(sig_hist, 1e-8), edges, color="#e45756", lw=2.0, label="signal")
        ax.set_yscale("log")
        ax.set_xlabel("DNN score")
        ax.set_ylabel("Expected weighted yield")
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.set_title(year)
        ax.legend(fontsize=9)
    fig.suptitle("SR DNN score templates from scored ROOT files", fontsize=15)
    save(fig, "sr_dnn_score_signal_vs_background_2017_2018")


def plot_cut_significance() -> None:
    thresholds = np.linspace(0.0, 0.98, 100)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, year in zip(axes, YEARS):
        base = scored_base(year)
        sig_score, sig_w = selected_scores(base / "SR" / "signal.root", True)
        bkg_scores = []
        bkg_weights = []
        for process in BKG_PROCESSES:
            score, weight = selected_scores(base / "SR" / f"{process}.root", True)
            bkg_scores.append(score)
            bkg_weights.append(weight)
        bkg_score = np.concatenate(bkg_scores)
        bkg_w = np.concatenate(bkg_weights)

        s_vals = np.array([sig_w[sig_score >= cut].sum() for cut in thresholds])
        b_vals = np.array([bkg_w[bkg_score >= cut].sum() for cut in thresholds])
        z_vals = np.divide(s_vals, np.sqrt(np.maximum(b_vals, 1e-9)))
        best = int(np.nanargmax(z_vals))
        ax.plot(thresholds, z_vals, color="#2ca02c", lw=2)
        ax.axvline(thresholds[best], color="0.25", ls="--", lw=1.2)
        ax.text(
            0.04,
            0.93,
            f"best cut={thresholds[best]:.2f}\nS={s_vals[best]:.1f}, B={b_vals[best]:.1f}\nS/sqrt(B)={z_vals[best]:.2f}",
            transform=ax.transAxes,
            va="top",
            fontsize=10,
        )
        ax.set_xlabel("DNN score cut")
        ax.set_ylabel(r"Approx. $S/\sqrt{B}$ in SR")
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.set_title(year)
    fig.suptitle("Approximate SR sensitivity versus DNN score threshold", fontsize=15)
    save(fig, "sr_cut_significance_2017_2018")


def plot_input_examples() -> None:
    preferred = ["ptmiss", "mT", "ptmiss_significance_corrected", "dphi_visibles_ptmiss", "ll_pt", "jet_pt_0"]
    fig, axes = plt.subplots(len(preferred), 2, figsize=(11.5, 3.0 * len(preferred)), sharey=False)
    for col, year in enumerate(YEARS):
        data = load_test_scores(year)
        metadata = load_metadata(year)
        names = metadata["feature_names"]
        x = data["raw_features"]
        label = data["label"]
        weight = data["train_weight"]
        for row, name in enumerate(preferred):
            ax = axes[row, col]
            if name not in names:
                ax.axis("off")
                continue
            idx = names.index(name)
            lo, hi = np.nanpercentile(x[:, idx], [1, 99])
            if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
                lo, hi = np.nanmin(x[:, idx]), np.nanmax(x[:, idx]) + 1.0
            bins = np.linspace(lo, hi, 45)
            ax.hist(x[label == 0, idx], bins=bins, weights=weight[label == 0], density=True, histtype="stepfilled", alpha=0.25, color="#4c78a8", label="background")
            ax.hist(x[label == 1, idx], bins=bins, weights=weight[label == 1], density=True, histtype="step", lw=1.7, color="#e45756", label="signal")
            ax.set_xlabel(name)
            ax.set_ylabel("Normalized")
            ax.grid(True, alpha=0.25)
            if row == 0:
                ax.set_title(year)
                ax.legend(fontsize=8)
    fig.suptitle("Representative input-variable separation on held-out test sample", fontsize=15, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    save(fig, "input_variable_examples_2017_2018")


def write_summary() -> None:
    lines = [
        "# DNN Training Validation Plots",
        "",
        "These plots are intended to document the DNN training behavior and the",
        "separation power of the trained score.",
        "",
        "| year | train candidates | best validation AUC | test AUC | test loss |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for year in YEARS:
        summary = load_metadata(year)["training_summary"]
        lines.append(
            f"| {year} | {summary['num_training_candidates']} | "
            f"{summary['best_val_auc']:.6f} | {summary['test_auc']:.6f} | {summary['test_loss']:.6f} |"
        )
    lines.extend([
        "",
        "Generated figures:",
        "",
        "- `training_loss_auc_2017_2018.png`: training/validation loss and validation AUC.",
        "- `test_roc_2017_2018.png`: ROC curves on held-out test samples.",
        "- `test_score_signal_background_2017_2018.png`: test-sample DNN score separation.",
        "- `sr_dnn_score_signal_vs_background_2017_2018.png`: weighted SR score templates.",
        "- `sr_cut_significance_2017_2018.png`: approximate cumulative `S/sqrt(B)` versus score cut.",
        "- `input_variable_examples_2017_2018.png`: representative input-variable separation.",
    ])
    ensure_out()
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    plot_history()
    plot_roc()
    plot_test_score()
    plot_sr_score_templates()
    plot_cut_significance()
    plot_input_examples()
    write_summary()
    print(f"Wrote DNN validation plots to {OUT}")


if __name__ == "__main__":
    main()
