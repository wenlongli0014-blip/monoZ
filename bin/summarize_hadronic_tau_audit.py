#!/usr/bin/env python3

"""Summarize TauVetoAudit and final Vars yields for every production dataset."""

import argparse
import csv
import json
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

DATA_DATASETS = {
    "DoubleEG", "DoubleMuon", "EGamma", "MuonEG", "SingleElectron",
    "SingleMuon",
}


def dataset_from_job_name(job_name):
    dataset, separator, job_id = job_name.rpartition("_")
    if not separator or not job_id.isdigit():
        raise ValueError(f"Cannot parse job name {job_name}")
    return dataset


def empty_summary():
    return {
        "num_jobs": 0,
        "audit_entries": 0,
        "audit_sum_weight": 0.0,
        "audit_sum_weight2": 0.0,
        "tau_rejected_events": 0,
        "tau_rejected_sum_weight": 0.0,
        "tau_rejected_sum_weight2": 0.0,
        "tau_preselected_candidates": 0,
        "tau_selected_candidates": 0,
        "tau_vvvloose_candidates": 0,
        "tau_vloose_candidates": 0,
        "tau_medium_candidates": 0,
        "tau_genuine_candidates": 0,
        "tau_electron_candidates": 0,
        "tau_muon_candidates": 0,
        "tau_unmatched_candidates": 0,
        "vars_entries": 0,
        "vars_sum_weight": 0.0,
        "vars_sum_weight2": 0.0,
    }


def summarize_dataset(paths, is_data):
    summary = empty_summary()
    summary["num_jobs"] = len(paths)

    audit_chain = ROOT.TChain("TauVetoAudit")
    vars_chain = ROOT.TChain("Vars")
    for path in paths:
        if audit_chain.Add(str(path)) != 1:
            raise OSError(f"Could not add TauVetoAudit from {path}")
        if vars_chain.Add(str(path)) != 1:
            raise OSError(f"Could not add Vars from {path}")

    audit = ROOT.RDataFrame(audit_chain).Define(
        "tau_audit_weight2", "weight * weight"
    )
    rejected = audit.Filter("!pass_tau_veto")
    actions = {
        "audit_entries": audit.Count(),
        "audit_sum_weight": audit.Sum("weight"),
        "audit_sum_weight2": audit.Sum("tau_audit_weight2"),
        "tau_rejected_events": rejected.Count(),
        "tau_rejected_sum_weight": rejected.Sum("weight"),
        "tau_rejected_sum_weight2": rejected.Sum("tau_audit_weight2"),
        "tau_preselected_candidates": audit.Sum("n_tau_preselected"),
        "tau_selected_candidates": audit.Sum("n_tau_selected"),
        "tau_vvvloose_candidates": audit.Sum("n_tau_vvvloose"),
        "tau_vloose_candidates": audit.Sum("n_tau_vloose"),
        "tau_medium_candidates": audit.Sum("n_tau_medium"),
    }
    if not is_data:
        actions.update({
            "tau_genuine_candidates": audit.Sum("n_tau_genuine"),
            "tau_electron_candidates": audit.Sum("n_tau_electron"),
            "tau_muon_candidates": audit.Sum("n_tau_muon"),
            "tau_unmatched_candidates": audit.Sum("n_tau_unmatched"),
        })
    else:
        for key in [
            "tau_genuine_candidates",
            "tau_electron_candidates",
            "tau_muon_candidates",
            "tau_unmatched_candidates",
        ]:
            summary[key] = None

    for key, action in actions.items():
        value = action.GetValue()
        if "weight" in key:
            summary[key] = float(value)
        else:
            summary[key] = int(value)

    vars_frame = ROOT.RDataFrame(vars_chain)
    summary["vars_entries"] = int(vars_frame.Count().GetValue())
    if is_data:
        summary["vars_sum_weight"] = float(summary["vars_entries"])
        summary["vars_sum_weight2"] = float(summary["vars_entries"])
    else:
        weighted_vars = vars_frame.Define(
            "tau_study_vars_weight2", "weight * weight"
        )
        summary["vars_sum_weight"] = float(weighted_vars.Sum("weight").GetValue())
        summary["vars_sum_weight2"] = float(
            weighted_vars.Sum("tau_study_vars_weight2").GetValue()
        )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default=(
            "/eos/user/l/liwe/hadronic_tau_test/metadata/"
            "production_manifest.json"
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default="/eos/user/l/liwe/hadronic_tau_test/metadata/tau_audit_summary",
    )
    parser.add_argument("--variants", nargs="+")
    parser.add_argument("--years", nargs="+")
    parser.add_argument("--tasks", nargs="+")
    parser.add_argument("--datasets", nargs="+")
    parser.add_argument(
        "--merge-existing",
        action="store_true",
        help="Replace selected rows in an existing output JSON.",
    )
    args = parser.parse_args()

    with open(args.manifest) as stream:
        production = json.load(stream)

    dataset_paths = {}
    for task in production["tasks"]:
        if args.variants and task["variant"] not in args.variants:
            continue
        if args.years and task["year"] not in args.years:
            continue
        if args.tasks and task["task"] not in args.tasks:
            continue
        task_dir = Path(task["task_dir"])
        output_dir = Path(task["output_dir"])
        with open(task_dir / "job_names.dat") as stream:
            job_names = [line.strip() for line in stream if line.strip()]

        for job_name in job_names:
            dataset = dataset_from_job_name(job_name)
            if args.datasets and dataset not in args.datasets:
                continue
            key = (
                task["variant"],
                task["year"],
                task["task"],
                dataset,
            )
            path = output_dir / f"{job_name}.root"
            if not path.is_file():
                raise FileNotFoundError(path)
            dataset_paths.setdefault(key, []).append(path)

    summaries = {}
    for index, (key, paths) in enumerate(sorted(dataset_paths.items()), start=1):
        print(
            f"[{index}/{len(dataset_paths)}] "
            f"{' '.join(key)} ({len(paths)} files)"
        )
        summaries[key] = summarize_dataset(
            paths, is_data=key[-1] in DATA_DATASETS
        )

    rows = []
    for key in sorted(summaries):
        row = dict(
            zip(["variant", "year", "task", "dataset"], key)
        )
        row.update(summaries[key])
        rows.append(row)

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    if args.merge_existing:
        json_path = output_prefix.with_suffix(".json")
        if not json_path.is_file():
            raise FileNotFoundError(json_path)
        with open(json_path) as stream:
            existing_rows = json.load(stream)
        replacement_keys = {
            (row["variant"], row["year"], row["task"], row["dataset"])
            for row in rows
        }
        rows.extend(
            row for row in existing_rows
            if (
                row["variant"],
                row["year"],
                row["task"],
                row["dataset"],
            ) not in replacement_keys
        )
        rows.sort(
            key=lambda row: (
                row["variant"],
                row["year"],
                row["task"],
                row["dataset"],
            )
        )
    with open(output_prefix.with_suffix(".json"), "w") as stream:
        json.dump(rows, stream, indent=2)
        stream.write("\n")
    with open(output_prefix.with_suffix(".csv"), "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_prefix}.csv and .json")


if __name__ == "__main__":
    main()
