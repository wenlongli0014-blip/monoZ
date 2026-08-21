#!/usr/bin/env python3

"""Prepare the full hadronic-tau veto study production."""

import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys


VARIANTS = {
    "tau_medium": ["--tau-vsjet-wp=16"],
    "tau_vloose": ["--tau-vsjet-wp=4"],
    "tau_vvvloose": ["--tau-vsjet-wp=1"],
    "tau_none": ["--disable-tau-veto"],
}


def git_output(repo, *args):
    return subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-base",
        default="/eos/user/l/liwe/hadronic_tau_test",
    )
    parser.add_argument(
        "--task-base",
        help="AFS directory for submit files and logs.",
    )
    parser.add_argument("--events-per-job", type=int, default=500000)
    parser.add_argument(
        "--years", nargs="+", choices=["2017", "2018"], default=["2017", "2018"]
    )
    parser.add_argument(
        "--variants", nargs="+", choices=VARIANTS, default=list(VARIANTS)
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse task directories that already contain job_names.dat.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate existing task scripts and submit files.",
    )
    args = parser.parse_args()

    repo = Path(os.environ["HZZ2L2NU_BASE"]).resolve()
    output_base = Path(os.path.abspath(args.output_base))
    if args.task_base:
        task_base = Path(args.task_base).resolve()
    else:
        task_base = repo / "tasks" / "hadronic_tau_test"
    prepare = repo / "bin" / "prepare_htcondor_jobs.py"
    tasks = []

    for variant in args.variants:
        for year in args.years:
            definitions = [
                (
                    "dilepton_jobs",
                    f"dilepton/{year}-ul.yaml",
                    f"config/samples_dilepton_{year}.txt",
                    "DileptonTrees",
                ),
                (
                    "trilepton_jobs",
                    f"trilepton/{year}-ul.yaml",
                    f"config/samples_trilepton_{year}.txt",
                    "TrileptonTrees",
                ),
                (
                    "signal_jobs",
                    f"dilepton/{year}-ul-signal-monoZ.yaml",
                    f"config/samples_zh_invisible_{year}_shears_skim.txt",
                    "DileptonTrees",
                ),
            ]

            for task_name, config, samples, analysis in definitions:
                task_dir = task_base / variant / year / task_name
                output_dir = output_base / variant / year / task_name / "output"
                reused = (task_dir / "job_names.dat").exists()
                if reused and not (args.resume or args.overwrite):
                    raise RuntimeError(f"Task is already prepared: {task_dir}")

                command = [
                    sys.executable,
                    str(prepare),
                    "--task-dir",
                    str(task_dir),
                    "--output-dir",
                    str(output_dir),
                    "--config",
                    config,
                    "--events-perjob",
                    str(args.events_per_job),
                    "--",
                    str(repo / samples),
                    "-a",
                    analysis,
                    "--more-vars",
                    *VARIANTS[variant],
                ]
                if not reused or args.overwrite:
                    subprocess.run(command, cwd=repo, check=True)

                with open(task_dir / "job_names.dat") as stream:
                    num_jobs = sum(1 for line in stream if line.strip())
                tasks.append(
                    {
                        "variant": variant,
                        "year": year,
                        "task": task_name,
                        "task_dir": str(task_dir),
                        "output_dir": str(output_dir),
                        "num_jobs": num_jobs,
                        "reused": reused and not args.overwrite,
                        "command": command,
                    }
                )

    metadata_dir = output_base / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    manifest = {
        "created_utc": timestamp,
        "repository": str(repo),
        "git_commit": git_output(repo, "rev-parse", "HEAD"),
        "git_status": git_output(repo, "status", "--short"),
        "events_per_job": args.events_per_job,
        "tasks": tasks,
        "total_jobs": sum(task["num_jobs"] for task in tasks),
    }
    with open(metadata_dir / "production_manifest.json", "w") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")

    task_base.mkdir(parents=True, exist_ok=True)
    submit_path = task_base / "submit_all.sh"
    with open(submit_path, "w") as stream:
        stream.write("#!/bin/bash\nset -euo pipefail\n")
        for task in tasks:
            stream.write(
                f"condor_submit {task['task_dir']}/jobs.sub\n"
            )
    submit_path.chmod(0o755)

    print(
        f"Prepared {len(tasks)} tasks with {manifest['total_jobs']} jobs. "
        f"Submit with {submit_path}"
    )


if __name__ == "__main__":
    main()
