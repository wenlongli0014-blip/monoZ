#!/usr/bin/env python3

"""Validate every ROOT output listed in the tau-veto production manifest."""

import argparse
import datetime
import json
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

REQUIRED_BRANCHES = {
    "Vars": {"run", "event", "ptmiss", "lepton_cat"},
    "TauVetoAudit": {
        "run",
        "lumi",
        "event",
        "weight",
        "tau_veto_applied",
        "pass_tau_veto",
        "tau_vsjet_wp",
        "n_tau_preselected",
        "n_tau_selected",
        "n_tau_genuine",
        "tau_pt",
        "tau_id_vsjet",
        "tau_gen_part_flavour",
    },
}


def validate_file(path):
    problems = []
    entries = {}
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        return ["open_or_zombie"], entries
    if root_file.TestBit(ROOT.TFile.kRecovered):
        problems.append("recovered")

    for tree_name, required_branches in REQUIRED_BRANCHES.items():
        tree = root_file.Get(tree_name)
        if not tree or not tree.InheritsFrom("TTree"):
            problems.append(f"missing_{tree_name}")
            continue
        branches = {
            branch.GetName() for branch in tree.GetListOfBranches()
        }
        missing = sorted(required_branches - branches)
        if missing:
            problems.append(
                f"{tree_name}_missing_branches:{','.join(missing)}"
            )
        entries[tree_name] = int(tree.GetEntries())

    root_file.Close()
    return problems, entries


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
        "--output",
        default=(
            "/eos/user/l/liwe/hadronic_tau_test/metadata/"
            "root_integrity_report.json"
        ),
    )
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    with open(manifest_path) as stream:
        production = json.load(stream)

    checked_files = 0
    total_entries = {tree_name: 0 for tree_name in REQUIRED_BRANCHES}
    problems = []
    for index, task in enumerate(production["tasks"], start=1):
        task_dir = Path(task["task_dir"])
        output_dir = Path(task["output_dir"])
        with open(task_dir / "job_names.dat") as stream:
            job_names = [line.strip() for line in stream if line.strip()]

        for job_name in job_names:
            path = output_dir / f"{job_name}.root"
            file_problems, entries = validate_file(path)
            for problem in file_problems:
                problems.append({"path": str(path), "problem": problem})
            for tree_name, count in entries.items():
                total_entries[tree_name] += count
            checked_files += 1

        print(
            f"[{index}/{len(production['tasks'])}] "
            f"{task['variant']} {task['year']} {task['task']}: "
            f"{len(job_names)}"
        )

    report = {
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat(),
        "manifest": str(manifest_path),
        "checked_files": checked_files,
        "expected_files": sum(
            task["num_jobs"] for task in production["tasks"]
        ),
        "total_entries": total_entries,
        "num_problems": len(problems),
        "problems": problems,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")

    print(json.dumps({key: value for key, value in report.items()
                      if key != "problems"}, indent=2))
    print(f"Report: {output_path}")
    if problems:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
