#!/usr/bin/env python3

"""Validate the grouped Run-2 mono-Z trees and write a JSON report."""

import argparse
import json
import math
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


PROCESSES = [
    "DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"
]
YEARS = ["2016", "2017", "2018"]
REGIONS = {
    "SR": "lepton_cat != 2 && ptmiss >= 100",
    "DY_CR": "lepton_cat != 2 && ptmiss >= 30 && ptmiss <= 90",
    "emu_CR": "lepton_cat == 2 && ptmiss >= 100",
    "3l_CR": None,
}
REQUIRED_COMMON = {
    "lepton_cat", "jet_cat", "ll_pt", "ll_mass", "ptmiss",
    "run", "lumi", "event",
}


def open_root(path):
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {path}")
    return root_file


def duplicate_events(tree):
    seen = set()
    duplicates = 0
    for event in tree:
        key = (int(event.run), int(event.lumi), int(event.event))
        if key in seen:
            duplicates += 1
        else:
            seen.add(key)
    return duplicates, len(seen)


def weight_summary(tree):
    dataframe = ROOT.RDataFrame(tree)
    handles = [
        dataframe.Sum("weight"),
        dataframe.Min("weight"),
        dataframe.Max("weight"),
        dataframe.Filter("!std::isfinite(weight)").Count(),
        dataframe.Filter("weight < 0").Count(),
    ]
    ROOT.RDF.RunGraphs(handles)
    return {
        "sum": float(handles[0].GetValue()),
        "min": float(handles[1].GetValue()),
        "max": float(handles[2].GetValue()),
        "nonfinite": int(handles[3].GetValue()),
        "negative": int(handles[4].GetValue()),
    }


def validate_file(path, region, is_data):
    root_file = open_root(path)
    recovered = bool(root_file.TestBit(ROOT.TFile.kRecovered))
    tree = root_file.Get("Vars")
    if not tree or not tree.InheritsFrom("TTree"):
        root_file.Close()
        raise OSError(f"Vars tree is missing in {path}")
    branches = [branch.GetName() for branch in tree.GetListOfBranches()]
    branch_set = set(branches)
    entries = int(tree.GetEntries())
    result = {
        "entries": entries,
        "branches": len(branches),
        "recovered": recovered,
        "missing_required": sorted(REQUIRED_COMMON - branch_set),
        "has_dnn_score": "dnn_score" in branch_set,
        "has_weight": "weight" in branch_set,
        "selection_violations": 0,
    }

    selection = REGIONS[region]
    if selection is not None:
        result["selection_violations"] = int(
            tree.GetEntries(f"!({selection})")
        )

    if is_data:
        duplicates, unique = duplicate_events(tree)
        result["duplicate_events"] = duplicates
        result["unique_events"] = unique
    elif "weight" in branch_set:
        result["weight"] = weight_summary(tree)

    root_file.Close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base", default="/eos/user/l/liwe/monoz_full_dnn"
    )
    parser.add_argument(
        "--output",
        default=(
            "/eos/user/l/liwe/monoz_full_dnn/metadata/"
            "final_validation.json"
        ),
    )
    args = parser.parse_args()

    base = Path(args.base)
    expected_paths = []
    for year in YEARS:
        for region in REGIONS:
            expected_paths.extend(
                base / year / region / f"{process}.root"
                for process in PROCESSES
            )
            expected_paths.append(base / year / region / "Data.root")
        expected_paths.append(base / year / "SR" / "signal.root")

    actual_paths = sorted(
        path
        for year in YEARS
        for path in (base / year).glob("*/*.root")
    )
    expected_set = {str(path) for path in expected_paths}
    actual_set = {str(path) for path in actual_paths}

    report = {
        "base": str(base),
        "expected_files": len(expected_paths),
        "actual_files": len(actual_paths),
        "missing_files": sorted(expected_set - actual_set),
        "extra_files": sorted(actual_set - expected_set),
        "part_files": [str(path) for path in base.rglob("*.part")],
        "files": {},
        "summary": {},
    }

    for path in expected_paths:
        relative = path.relative_to(base)
        year, region, filename = relative.parts
        is_data = filename == "Data.root"
        report["files"][str(relative)] = validate_file(
            path, region, is_data
        )

    files = list(report["files"].values())
    report["summary"] = {
        "recovered_files": sum(item["recovered"] for item in files),
        "files_with_missing_branches": sum(
            bool(item["missing_required"]) for item in files
        ),
        "files_with_dnn_score": sum(item["has_dnn_score"] for item in files),
        "selection_violations": sum(
            item["selection_violations"] for item in files
        ),
        "data_duplicate_events": sum(
            item.get("duplicate_events", 0) for item in files
        ),
        "nonfinite_weights": sum(
            item.get("weight", {}).get("nonfinite", 0) for item in files
        ),
        "total_entries": sum(item["entries"] for item in files),
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
