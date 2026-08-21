#!/usr/bin/env python3

"""Validate all grouped ROOT inputs produced for the tau-veto Combine study."""

import argparse
import json
import time
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


VARIANTS = ["tau_medium", "tau_vloose", "tau_vvvloose", "tau_none"]
YEARS = ["2017", "2018"]
REGIONS = ["SR", "DY_CR", "emu_CR", "3l_CR"]
PROCESSES = [
    "DY",
    "Other",
    "ST",
    "VVV",
    "WW",
    "WZ",
    "ggZZ",
    "qqZZ",
    "ttbar",
    "Data",
]


def open_root(path, attempts=5):
    for attempt in range(attempts):
        root_file = ROOT.TFile.Open(str(path))
        if root_file and not root_file.IsZombie():
            return root_file
        if root_file:
            root_file.Close()
        if attempt + 1 < attempts:
            time.sleep(2)
    return None


def expected_paths(base):
    for variant in VARIANTS:
        for year in YEARS:
            for region in REGIONS:
                for process in PROCESSES:
                    yield (
                        variant,
                        year,
                        region,
                        process,
                        base / variant / "combine_input" / year
                        / region / f"{process}.root",
                    )
            yield (
                variant,
                year,
                "SR",
                "signal",
                base / variant / "combine_input" / year
                / "SR" / "signal.root",
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="/eos/user/l/liwe/hadronic_tau_test",
        type=Path,
    )
    parser.add_argument(
        "--report",
        default=(
            "/eos/user/l/liwe/hadronic_tau_test/metadata/"
            "region_integrity_report.json"
        ),
        type=Path,
    )
    args = parser.parse_args()

    rows = []
    problems = []
    for variant, year, region, process, path in expected_paths(args.base):
        row = {
            "variant": variant,
            "year": year,
            "region": region,
            "process": process,
            "path": str(path),
        }
        if not path.is_file():
            row["problem"] = "missing"
            problems.append(row)
            rows.append(row)
            continue

        root_file = open_root(path)
        if not root_file:
            row["problem"] = "zombie_or_unreadable"
            problems.append(row)
            rows.append(row)
            continue
        if root_file.TestBit(ROOT.TFile.kRecovered):
            row["problem"] = "recovered"
            problems.append(row)
            root_file.Close()
            rows.append(row)
            continue

        tree = root_file.Get("Vars")
        if not tree or not tree.InheritsFrom("TTree"):
            row["problem"] = "missing_Vars"
            problems.append(row)
            root_file.Close()
            rows.append(row)
            continue

        entries = int(tree.GetEntries())
        branches = tree.GetListOfBranches()
        row["entries"] = entries
        row["branches"] = int(branches.GetEntries())
        if not tree.GetBranch("ptmiss"):
            row["problem"] = "missing_ptmiss"
            problems.append(row)
        elif entries:
            # Force ROOT to read an actual basket, rather than only metadata.
            read_bytes = int(tree.GetEntry(entries // 2))
            row["probe_read_bytes"] = read_bytes
            if read_bytes <= 0:
                row["problem"] = "branch_read_failed"
                problems.append(row)
        root_file.Close()
        rows.append(row)

    report = {
        "expected_files": len(rows),
        "checked_files": sum("entries" in row for row in rows),
        "total_entries": sum(row.get("entries", 0) for row in rows),
        "num_problems": len(problems),
        "problems": problems,
        "files": rows,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({key: report[key] for key in (
        "expected_files",
        "checked_files",
        "total_entries",
        "num_problems",
    )}, indent=2))
    print(f"Report: {args.report}")
    if problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
