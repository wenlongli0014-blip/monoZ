#!/usr/bin/env python3

"""Report the pre-training Tiwari-style AUC/Asimov feature rankings.

The ranking itself is now computed inside train_dnn_score.py before the
train/validation/test split. This script only combines and displays those
saved results; it does not run the legacy post-training permutation ranking.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


BASE = Path("/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn")
DEFAULT_OUTPUTS = {
    "2017": BASE / "2017" / "output",
    "2018": BASE / "2018" / "output",
    "combined_2017_2018": BASE / "combined_2017_2018" / "output",
}


def main():
    parser = argparse.ArgumentParser(description="Collect saved pre-training Asimov feature rankings.")
    parser.add_argument(
        "--out",
        default=str(BASE / "feature_importance" / "asimov_feature_ranking_2017_2018_combined.csv"),
    )
    parser.add_argument("--top", type=int, default=15)
    args = parser.parse_args()

    all_rows = []
    for dataset, outdir in DEFAULT_OUTPUTS.items():
        rows = json.loads((outdir / "feature_significance.json").read_text())
        print(f"\n{dataset}")
        for row in rows:
            enriched = {"dataset": dataset, **row}
            all_rows.append(enriched)
        for row in rows[: args.top]:
            print(
                f"{row['rank']:2d}. {row['feature']:32s} "
                f"AUC={row['auc']:.6f} Z20={row['asimov_z_syst']:.6f}"
            )

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(all_rows[0].keys())
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
