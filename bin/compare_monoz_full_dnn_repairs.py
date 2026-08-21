#!/usr/bin/env python3

"""Compare repaired dataset outputs with the incomplete initial production."""

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import uproot


DATASETS = [
    ("2016HIPM", "DYJetsToLL_PtZ-50To100"),
    ("2016HIPM", "DYJetsToLL_PtZ-250To400"),
    ("2016noHIPM", "DYJetsToLL_PtZ-100To250"),
]
TASKS = ["dilepton", "trilepton"]


def read_events(paths):
    identifiers = []
    weight_sum = 0.0
    for path in paths:
        with uproot.open(f"{path}:Vars") as tree:
            arrays = tree.arrays(["run", "lumi", "event"], library="np")
            identifiers.extend(
                zip(
                    arrays["run"].astype(np.int64),
                    arrays["lumi"].astype(np.int64),
                    arrays["event"].astype(np.int64),
                )
            )
            if "weight" in tree:
                weight_sum += float(
                    np.sum(tree["weight"].array(library="np"), dtype=np.float64)
                )
    event_set = set(identifiers)
    return {
        "entries": len(identifiers),
        "unique_events": len(event_set),
        "duplicate_entries": len(identifiers) - len(event_set),
        "events": event_set,
        "weight_sum": weight_sum,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--production",
        default="/eos/user/l/liwe/monoz_full_dnn/production",
    )
    parser.add_argument(
        "--repairs", default="/eos/user/l/liwe/monoz_full_dnn/repairs"
    )
    parser.add_argument(
        "--output",
        default=(
            "/eos/user/l/liwe/monoz_full_dnn/metadata/"
            "repair_vs_initial_production.json"
        ),
    )
    args = parser.parse_args()

    report = {}
    for period, dataset in DATASETS:
        for task in TASKS:
            old_paths = sorted(
                glob.glob(
                    f"{args.production}/{period}/{task}/output/"
                    f"{dataset}_[0-9]*.root"
                )
            )
            repaired_paths = sorted(
                glob.glob(
                    f"{args.repairs}/{period}/{task}/output/"
                    f"{dataset}_[0-9]*.root"
                )
            )
            old = read_events(old_paths)
            repaired = read_events(repaired_paths)
            key = f"{period}/{task}/{dataset}"
            report[key] = {
                "initial_files": len(old_paths),
                "repaired_files": len(repaired_paths),
                "initial_entries": old["entries"],
                "repaired_entries": repaired["entries"],
                "initial_duplicate_entries": old["duplicate_entries"],
                "repaired_duplicate_entries": repaired["duplicate_entries"],
                "intersection": len(old["events"] & repaired["events"]),
                "initial_only": len(old["events"] - repaired["events"]),
                "repaired_only": len(repaired["events"] - old["events"]),
                "initial_weight_sum": old["weight_sum"],
                "repaired_weight_sum": repaired["weight_sum"],
            }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
