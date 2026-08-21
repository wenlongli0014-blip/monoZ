#!/usr/bin/env python3

"""Validate the post-audit 2016 dataset repair jobs and ROOT outputs."""

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True

TASKS = {
    "2016HIPM/dilepton": 36,
    "2016HIPM/trilepton": 36,
    "2016noHIPM/dilepton": 20,
    "2016noHIPM/trilepton": 20,
}


def dataset_name(job_name):
    return job_name.rsplit("_", 1)[0]


def validate_root(path):
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        return {"error": "zombie_or_unopenable"}
    recovered = bool(root_file.TestBit(ROOT.TFile.kRecovered))
    tree = root_file.Get("Vars")
    if not tree or not tree.InheritsFrom("TTree"):
        root_file.Close()
        return {"error": "missing_Vars", "recovered": recovered}
    branches = tuple(
        (branch.GetName(), branch.GetClassName(), branch.GetTitle())
        for branch in tree.GetListOfBranches()
    )
    result = {
        "entries": int(tree.GetEntries()),
        "recovered": recovered,
        "branches": branches,
        "has_dnn_score": bool(tree.GetBranch("dnn_score")),
    }
    root_file.Close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-base", default="tasks/monoz_full_dnn_repair"
    )
    parser.add_argument(
        "--output-base",
        default="/eos/user/l/liwe/monoz_full_dnn/repairs",
    )
    parser.add_argument(
        "--report",
        default=(
            "/eos/user/l/liwe/monoz_full_dnn/metadata/"
            "repair_job_validation.json"
        ),
    )
    args = parser.parse_args()

    report = {"tasks": {}, "summary": {}}
    all_failures = []
    all_entries = defaultdict(int)
    schemas = {}
    total_expected = 0
    total_valid = 0

    for task, expected_count in TASKS.items():
        task_dir = Path(args.task_base) / task
        output_dir = Path(args.output_base) / task / "output"
        job_names = [
            line.strip()
            for line in (task_dir / "job_names.dat").read_text().splitlines()
            if line.strip()
        ]
        expected_files = {f"{name}.root" for name in job_names}
        actual_files = {path.name for path in output_dir.glob("*.root")}
        allowed_extra = {"MuonEG_0.root"} if "2016noHIPM" in task else set()
        missing = sorted(expected_files - actual_files)
        unexpected = sorted(actual_files - expected_files - allowed_extra)

        normal_logs = 0
        bad_logs = []
        for job_name in job_names:
            matches = list(
                (task_dir / "jobs" / "logs").glob(
                    f"runOnBatch_{job_name}.*.log"
                )
            )
            if len(matches) != 1:
                bad_logs.append(
                    {"job": job_name, "reason": f"log_count_{len(matches)}"}
                )
                continue
            log_text = matches[0].read_text(errors="replace")
            if "Normal termination (return value 0)" in log_text:
                normal_logs += 1
            else:
                bad_logs.append({"job": job_name, "reason": "not_return_0"})

            for suffix in ("outfile", "errors"):
                streams = list(
                    (task_dir / "jobs" / "logs").glob(
                        f"runOnBatch_{job_name}.{suffix}.*.txt"
                    )
                )
                for stream in streams:
                    text = stream.read_text(errors="replace")
                    if re.search(
                        r"\[ERROR\]|FATAL|Traceback|segmentation|"
                        r"Invalid address",
                        text,
                        re.IGNORECASE,
                    ):
                        bad_logs.append(
                            {
                                "job": job_name,
                                "reason": f"error_text_in_{stream.name}",
                            }
                        )

        root_failures = []
        dataset_entries = defaultdict(int)
        for filename in sorted(expected_files):
            path = output_dir / filename
            if not path.exists():
                continue
            result = validate_root(path)
            if result.get("error") or result.get("recovered"):
                root_failures.append({"file": filename, **result})
                continue
            if result["has_dnn_score"]:
                root_failures.append(
                    {"file": filename, "error": "contains_dnn_score"}
                )
                continue
            dataset = dataset_name(filename[:-5])
            schema_key = (task.split("/")[1], dataset)
            if schema_key in schemas and schemas[schema_key] != result["branches"]:
                root_failures.append(
                    {"file": filename, "error": "schema_mismatch"}
                )
                continue
            schemas[schema_key] = result["branches"]
            dataset_entries[dataset] += result["entries"]
            all_entries[f"{task}/{dataset}"] += result["entries"]
            total_valid += 1

        task_report = {
            "expected_jobs": expected_count,
            "job_names": len(job_names),
            "normal_exit_0": normal_logs,
            "missing_outputs": missing,
            "unexpected_outputs": unexpected,
            "allowed_preexisting_outputs": sorted(actual_files & allowed_extra),
            "bad_logs": bad_logs,
            "root_failures": root_failures,
            "dataset_entries": dict(sorted(dataset_entries.items())),
        }
        report["tasks"][task] = task_report
        total_expected += expected_count
        if (
            len(job_names) != expected_count or normal_logs != expected_count
            or missing or unexpected or bad_logs or root_failures
        ):
            all_failures.append(task)

    raw_scripts = sorted(
        (Path(args.task_base) / "2016HIPM" / "dilepton" / "jobs" / "scripts")
        .glob("runOnBatch_DYJetsToLL_PtZ-250To400_*.sh")
    )
    coverage = Counter()
    for script in raw_scripts:
        text = script.read_text()
        skip = int(re.search(r"--skip-files=(\d+)", text).group(1))
        maximum = int(re.search(r"--max-files=(\d+)", text).group(1))
        for index in range(skip, min(skip + maximum, 31)):
            coverage[index] += 1
    report["raw_UL16APV_coverage"] = {
        "DAS_files": 31,
        "jobs": len(raw_scripts),
        "missing_indices": [i for i in range(31) if coverage[i] == 0],
        "duplicate_indices": [i for i in range(31) if coverage[i] > 1],
    }
    if (
        report["raw_UL16APV_coverage"]["missing_indices"]
        or report["raw_UL16APV_coverage"]["duplicate_indices"]
    ):
        all_failures.append("raw_UL16APV_coverage")

    report["summary"] = {
        "expected_jobs": total_expected,
        "valid_root_outputs": total_valid,
        "failed_tasks": all_failures,
        "dataset_entries": dict(sorted(all_entries.items())),
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    print(json.dumps(report["raw_UL16APV_coverage"], indent=2, sort_keys=True))
    print(f"Wrote {report_path}")
    if all_failures or total_valid != total_expected:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
