#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve


def main():
    parser = argparse.ArgumentParser(description="Plot one HZZ2l2nu DNN training output directory.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--title", default="HZZ2l2nu DNN")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    history = np.loadtxt(outdir / "training_history.csv", delimiter=",", skiprows=1)
    if history.ndim == 1:
        history = history[None, :]
    scores = np.load(outdir / "split_scores.npz")
    ranking = json.loads((outdir / "feature_significance.json").read_text())
    metadata = json.loads((outdir / "dnn_metadata.json").read_text())
    summary = metadata["training_summary"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.4))

    ax = axes[0]
    ax.plot(history[:, 0], history[:, 1], label="train BCE", color="#3f90da")
    ax.plot(history[:, 0], history[:, 2], label="validation BCE", color="#bd1f01")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Weighted BCEWithLogits loss")
    ax.grid(alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(history[:, 0], history[:, 3], label="validation AUC", color="#00a087", alpha=0.85)
    ax2.set_ylabel("Validation AUC")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="best", fontsize=9)
    ax.set_title("Training history")

    ax = axes[1]
    for split, color in [("train", "#3f90da"), ("val", "#ffa90e"), ("test", "#bd1f01")]:
        y = scores[f"{split}_label"]
        score = scores[f"{split}_score"]
        weight = scores[f"{split}_weight"]
        fpr, tpr, _ = roc_curve(y, score, sample_weight=weight)
        auc_key = "validation_auc" if split == "val" else f"{split}_auc"
        ax.plot(fpr, tpr, color=color, label=f"{split} AUC={summary[auc_key]:.4f}")
    ax.plot([0, 1], [0, 1], "--", color="0.5", linewidth=1)
    ax.set_xlabel("Background efficiency")
    ax.set_ylabel("Signal efficiency")
    ax.set_title("Weighted ROC")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=9)

    ax = axes[2]
    shown = ranking[:15][::-1]
    ax.barh(
        [row["feature"] for row in shown],
        [row["asimov_z_syst"] for row in shown],
        color=["#bd1f01" if row["rank"] <= 10 else "#9e9e9e" for row in shown],
    )
    ax.set_xlabel("Binned Asimov Z with 20% background uncertainty")
    ax.set_title("Full-dataset feature ranking\n(known test-leakage issue)")
    ax.grid(axis="x", alpha=0.25)

    fig.suptitle(args.title, fontsize=15)
    fig.tight_layout()
    fig.savefig(outdir / "training_diagnostics.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(outdir / "training_diagnostics.png")


if __name__ == "__main__":
    main()
