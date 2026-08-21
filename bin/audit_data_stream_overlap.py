#!/usr/bin/env python3

"""Measure exact event overlap between primary data streams."""

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import uproot


DATASETS = {
    "2017": ["DoubleEG", "DoubleMuon", "SingleElectron", "SingleMuon", "MuonEG"],
    "2018": ["EGamma", "DoubleMuon", "SingleMuon", "MuonEG"],
}

REGIONS = {
    "SR": lambda arrays: (arrays["lepton_cat"] != 2)
    & (arrays["ptmiss"] >= 100),
    "DY_CR": lambda arrays: (arrays["lepton_cat"] != 2)
    & (arrays["ptmiss"] >= 30)
    & (arrays["ptmiss"] <= 90),
    "emu_CR": lambda arrays: (arrays["lepton_cat"] == 2)
    & (arrays["ptmiss"] >= 100),
}

EVENT_DTYPE = np.dtype(
    [("run", np.uint32), ("lumi", np.uint32), ("event", np.uint64)]
)


def event_ids(arrays, mask=None):
    if mask is None:
        mask = slice(None)
    result = np.empty(len(arrays["run"][mask]), dtype=EVENT_DTYPE)
    result["run"] = arrays["run"][mask]
    result["lumi"] = arrays["lumi"][mask]
    result["event"] = arrays["event"][mask]
    return result


def unique_event_ids(ids):
    if not len(ids):
        return ids
    order = np.lexsort((ids["event"], ids["lumi"], ids["run"]))
    sorted_ids = ids[order]
    keep = np.r_[
        True,
        (sorted_ids["run"][1:] != sorted_ids["run"][:-1])
        | (sorted_ids["lumi"][1:] != sorted_ids["lumi"][:-1])
        | (sorted_ids["event"][1:] != sorted_ids["event"][:-1]),
    ]
    return sorted_ids[keep]


def load_audit_ids(paths):
    chunks = []
    expressions = ["run", "lumi", "event"]
    for arrays in uproot.iterate(
        [f"{path}:TauVetoAudit" for path in paths],
        expressions=expressions,
        library="np",
        step_size="200 MB",
    ):
        chunks.append(event_ids(arrays))
    return np.concatenate(chunks) if chunks else np.empty(0, dtype=EVENT_DTYPE)


def load_vars_ids(paths, task):
    chunks = {"Vars": []}
    if task == "dilepton_jobs":
        chunks.update({region: [] for region in REGIONS})

    expressions = ["run", "lumi", "event"]
    if task == "dilepton_jobs":
        expressions.extend(["lepton_cat", "ptmiss"])

    for arrays in uproot.iterate(
        [f"{path}:Vars" for path in paths],
        expressions=expressions,
        library="np",
        step_size="100 MB",
    ):
        chunks["Vars"].append(event_ids(arrays))
        if task == "dilepton_jobs":
            for region, selection in REGIONS.items():
                chunks[region].append(event_ids(arrays, selection(arrays)))

    return {
        name: np.concatenate(region_chunks)
        if region_chunks
        else np.empty(0, dtype=EVENT_DTYPE)
        for name, region_chunks in chunks.items()
    }


def summarize_overlap(stream_ids):
    unique_by_stream = {}
    streams = {}
    for stream, ids in stream_ids.items():
        unique = unique_event_ids(ids)
        unique_by_stream[stream] = unique
        streams[stream] = {
            "entries": int(len(ids)),
            "unique_events": int(len(unique)),
            "duplicates_within_stream": int(len(ids) - len(unique)),
        }

    stream_names = list(unique_by_stream)
    all_ids = np.concatenate(list(unique_by_stream.values()))
    all_bits = np.concatenate(
        [
            np.full(len(unique_by_stream[stream]), 1 << index, dtype=np.uint16)
            for index, stream in enumerate(stream_names)
        ]
    )
    order = np.lexsort((all_ids["event"], all_ids["lumi"], all_ids["run"]))
    sorted_ids = all_ids[order]
    group_starts = np.r_[
        0,
        np.flatnonzero(
            (sorted_ids["run"][1:] != sorted_ids["run"][:-1])
            | (sorted_ids["lumi"][1:] != sorted_ids["lumi"][:-1])
            | (sorted_ids["event"][1:] != sorted_ids["event"][:-1])
        )
        + 1,
    ]
    group_ends = np.r_[group_starts[1:], len(sorted_ids)]
    membership = group_ends - group_starts
    stream_masks = np.bitwise_or.reduceat(all_bits[order], group_starts)

    pairs = []
    for first_index, first in enumerate(stream_names):
        for second_index, second in enumerate(
            stream_names[first_index + 1 :], start=first_index + 1
        ):
            pair_mask = (1 << first_index) | (1 << second_index)
            pairs.append(
                {
                    "first": first,
                    "second": second,
                    "overlap_events": int(
                        np.count_nonzero((stream_masks & pair_mask) == pair_mask)
                    ),
                }
            )

    multiplicities, counts = np.unique(membership, return_counts=True)
    return {
        "streams": streams,
        "sum_stream_unique_events": int(len(all_ids)),
        "union_unique_events": int(len(membership)),
        "duplicate_stream_entries": int(len(all_ids) - len(membership)),
        "events_in_multiple_streams": int(np.count_nonzero(membership > 1)),
        "max_stream_multiplicity": int(membership.max()) if len(membership) else 0,
        "stream_multiplicity": {
            str(int(multiplicity)): int(count)
            for multiplicity, count in zip(multiplicities, counts)
        },
        "pair_overlaps": pairs,
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="/eos/user/l/liwe/hadronic_tau_test/tau_none",
        help="Tau-variant production directory",
    )
    parser.add_argument(
        "--scope",
        choices=["all", "audit", "vars"],
        default="all",
        help="Trees to inspect",
    )
    parser.add_argument(
        "--years", nargs="+", choices=sorted(DATASETS), default=sorted(DATASETS)
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        choices=["dilepton_jobs", "trilepton_jobs"],
        default=["dilepton_jobs", "trilepton_jobs"],
    )
    parser.add_argument("--output", help="Optional output JSON path")
    return parser.parse_args()


def main():
    args = parse_args()
    result = {"base": str(Path(args.base).resolve()), "years": {}}

    for year in args.years:
        datasets = DATASETS[year]
        year_result = {}
        for task in args.tasks:
            print(f"Reading {year} {task} ({args.scope})", file=sys.stderr)
            audit_by_stream = {}
            vars_by_region = {}
            for dataset in datasets:
                pattern = str(
                    Path(args.base)
                    / year
                    / task
                    / "output"
                    / f"{dataset}_[0-9]*.root"
                )
                paths = sorted(glob.glob(pattern))
                if not paths:
                    raise RuntimeError(f"No files match {pattern}")

                if args.scope in ("all", "audit"):
                    audit_by_stream[dataset] = load_audit_ids(paths)
                if args.scope in ("all", "vars"):
                    for region, ids in load_vars_ids(paths, task).items():
                        vars_by_region.setdefault(region, {})[dataset] = ids

            task_result = {}
            if audit_by_stream:
                task_result["TauVetoAudit"] = summarize_overlap(audit_by_stream)
            if vars_by_region:
                task_result["Vars"] = {
                    region: summarize_overlap(stream_ids)
                    for region, stream_ids in vars_by_region.items()
                }
            year_result[task] = task_result
        result["years"][year] = year_result

    encoded = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded + "\n")
    print(encoded)


if __name__ == "__main__":
    main()
