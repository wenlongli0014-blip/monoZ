#!/usr/bin/env python3

"""Audit Events and Runs normalization metadata for tau-study DDFs."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "python"))

from hzz.dataset import parse_datasets_file  # noqa: E402


SAMPLE_LISTS = [
    ("2017", "config/samples_dilepton_2017.txt",
     "dilepton/2017-ul.yaml"),
    ("2018", "config/samples_dilepton_2018.txt",
     "dilepton/2018-ul.yaml"),
    ("2017", "config/samples_zh_invisible_2017_shears_skim.txt",
     "dilepton/2017-ul-signal-monoZ.yaml"),
    ("2018", "config/samples_zh_invisible_2018_shears_skim.txt",
     "dilepton/2018-ul-signal-monoZ.yaml"),
]


def read_file_metadata(path, is_simulation):
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise OSError(f"Cannot open {path}")
    events = root_file.Get("Events")
    runs = root_file.Get("Runs")
    if not events or not runs:
        root_file.Close()
        raise RuntimeError(f"Events or Runs is missing in {path}")

    num_selected_events = int(events.GetEntries())
    gen_event_sumw = 0.0
    runs_signature = None
    if is_simulation:
        if not runs.GetBranch("genEventSumw"):
            root_file.Close()
            raise RuntimeError(f"Runs.genEventSumw is missing in {path}")
        signature_rows = []
        for entry in runs:
            sumw = float(entry.genEventSumw)
            gen_event_sumw += sumw
            signature_rows.append((
                int(entry.run),
                int(entry.genEventCount),
                sumw,
                float(entry.genEventSumw2),
            ))
        runs_signature = tuple(signature_rows)
    root_file.Close()
    return num_selected_events, gen_event_sumw, runs_signature


def relative_difference(actual, expected):
    scale = max(abs(actual), abs(expected), 1.0)
    return abs(actual - expected) / scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-prefix",
        default=(
            "/eos/user/l/liwe/hadronic_tau_test/metadata/"
            "ddf_normalization_audit"
        ),
    )
    parser.add_argument(
        "--sumw-tolerance",
        type=float,
        default=5e-5,
    )
    args = parser.parse_args()

    os.environ.setdefault("HZZ2L2NU_BASE", str(REPO))
    datasets = {}
    for year, sample_list, config in SAMPLE_LISTS:
        parsed = parse_datasets_file(
            str(REPO / sample_list), config_path=config
        )
        for dataset in parsed:
            datasets[(year, dataset.path)] = dataset

    rows = []
    for index, ((year, ddf_path), dataset) in enumerate(
        sorted(datasets.items()), start=1
    ):
        actual_selected = 0
        actual_sumw = 0.0
        raw_sumw = 0.0
        runs_signatures = set()
        duplicate_runs_files = 0
        errors = []
        for path in dataset.files:
            try:
                selected, sumw, runs_signature = read_file_metadata(
                    path, dataset.is_sim
                )
                actual_selected += selected
                raw_sumw += sumw
                if runs_signature is not None:
                    if runs_signature in runs_signatures:
                        duplicate_runs_files += 1
                    else:
                        runs_signatures.add(runs_signature)
                        actual_sumw += sumw
            except Exception as error:
                errors.append(str(error))

        expected_selected = dataset.parameters.get(
            "num_selected_events"
        )
        selected_match = (
            expected_selected is not None
            and not errors
            and actual_selected == int(expected_selected)
        )

        expected_sumw = None
        sumw_relative_difference = None
        raw_sumw_relative_difference = None
        deduplicated_sumw_relative_difference = None
        sumw_mode = None
        sumw_match = None
        if dataset.is_sim:
            num_events = dataset.parameters.get("num_events")
            mean_weight = dataset.parameters.get("mean_weight")
            if num_events is not None and mean_weight is not None:
                expected_sumw = float(num_events) * float(mean_weight)
                raw_sumw_relative_difference = relative_difference(
                    raw_sumw, expected_sumw
                )
                deduplicated_sumw_relative_difference = relative_difference(
                    actual_sumw, expected_sumw
                )
                if (
                    raw_sumw_relative_difference
                    <= deduplicated_sumw_relative_difference
                ):
                    actual_sumw = raw_sumw
                    sumw_mode = "raw"
                    sumw_relative_difference = (
                        raw_sumw_relative_difference
                    )
                else:
                    sumw_mode = "deduplicated"
                    sumw_relative_difference = (
                        deduplicated_sumw_relative_difference
                    )
                sumw_match = (
                    not errors
                    and math.isfinite(actual_sumw)
                    and sumw_relative_difference <= args.sumw_tolerance
                )

        row = {
            "year": year,
            "name": dataset.name,
            "ddf": ddf_path,
            "is_sim": dataset.is_sim,
            "num_files": len(dataset.files),
            "expected_selected_events": expected_selected,
            "actual_selected_events": actual_selected,
            "selected_events_match": selected_match,
            "expected_gen_event_sumw": expected_sumw,
            "actual_gen_event_sumw": (
                actual_sumw if dataset.is_sim else None
            ),
            "raw_gen_event_sumw_before_deduplication": (
                raw_sumw if dataset.is_sim else None
            ),
            "unique_runs_metadata": (
                len(runs_signatures) if dataset.is_sim else None
            ),
            "duplicate_runs_files": (
                duplicate_runs_files if dataset.is_sim else None
            ),
            "sumw_mode": sumw_mode,
            "raw_sumw_relative_difference": raw_sumw_relative_difference,
            "deduplicated_sumw_relative_difference": (
                deduplicated_sumw_relative_difference
            ),
            "sumw_relative_difference": sumw_relative_difference,
            "sumw_match": sumw_match,
            "num_errors": len(errors),
            "errors": errors,
        }
        rows.append(row)
        print(
            f"[{index}/{len(datasets)}] {year} {dataset.name}: "
            f"selected={selected_match} sumw={sumw_match} "
            f"errors={len(errors)}"
        )

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(output_prefix.with_suffix(".json"), "w") as stream:
        json.dump(rows, stream, indent=2)
        stream.write("\n")
    with open(output_prefix.with_suffix(".csv"), "w", newline="") as stream:
        csv_rows = [
            {**row, "errors": " | ".join(row["errors"])}
            for row in rows
        ]
        writer = csv.DictWriter(
            stream, fieldnames=list(csv_rows[0])
        )
        writer.writeheader()
        writer.writerows(csv_rows)

    failures = [
        row for row in rows
        if row["num_errors"]
        or not row["selected_events_match"]
        or (row["is_sim"] and row["sumw_match"] is not True)
    ]
    print(
        f"Wrote {len(rows)} DDF rows; failures={len(failures)}: "
        f"{output_prefix}.json and .csv"
    )
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
