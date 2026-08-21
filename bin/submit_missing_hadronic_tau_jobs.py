#!/usr/bin/env python3

"""Prepare or submit only hadronic-tau jobs whose EOS output is missing."""

import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess


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
        "--submit",
        action="store_true",
        help="Submit the generated missing-output queues.",
    )
    parser.add_argument(
        "--jobs-per-task",
        type=int,
        default=None,
        help=(
            "Submit at most this many missing jobs from each task. This "
            "creates ordinary small clusters without late materialization."
        ),
    )
    parser.add_argument(
        "--job-flavour",
        default=None,
        help="Optional CERN HTCondor JobFlavour for the recovery clusters.",
    )
    args = parser.parse_args()
    if args.jobs_per_task is not None and args.jobs_per_task <= 0:
        parser.error("--jobs-per-task must be positive")

    manifest_path = Path(args.manifest)
    with open(manifest_path) as stream:
        production = json.load(stream)

    submissions = []
    total_jobs = 0
    total_missing_jobs = 0
    for task in production["tasks"]:
        task_dir = Path(task["task_dir"])
        output_dir = Path(task["output_dir"])
        with open(task_dir / "job_names.dat") as stream:
            job_names = [line.strip() for line in stream if line.strip()]

        all_missing = [
            name for name in job_names
            if not os.path.isfile(output_dir / f"{name}.root")
        ]
        if not all_missing:
            continue
        total_missing_jobs += len(all_missing)
        missing = all_missing
        if args.jobs_per_task is not None:
            missing = missing[:args.jobs_per_task]

        missing_names_path = task_dir / "job_names_missing.dat"
        with open(missing_names_path, "w") as stream:
            stream.write("\n".join(missing))
            stream.write("\n")

        submit_path = task_dir / "jobs_missing.sub"
        with open(task_dir / "jobs.sub") as source, open(submit_path, "w") as target:
            for line in source:
                if line.startswith("queue Job_name from "):
                    if args.job_flavour:
                        target.write(
                            f'+JobFlavour = "{args.job_flavour}"\n'
                        )
                    target.write(f"queue Job_name from {missing_names_path}\n")
                else:
                    target.write(line)

        record = {
            "variant": task["variant"],
            "year": task["year"],
            "task": task["task"],
            "submit_file": str(submit_path),
            "num_jobs": len(missing),
            "job_names": missing,
            "total_missing_jobs": len(all_missing),
            "deferred_jobs": len(all_missing) - len(missing),
        }
        if args.submit:
            result = subprocess.run(
                ["condor_submit", "-terse", str(submit_path)],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
            )
            record["condor_result"] = result.stdout.strip()
            print(
                f"Submitted {len(missing):4d} jobs for "
                f"{task['variant']} {task['year']} {task['task']}: "
                f"{record['condor_result']}"
            )
        else:
            print(
                f"Would submit {len(missing):4d} jobs for "
                f"{task['variant']} {task['year']} {task['task']}"
            )

        submissions.append(record)
        total_jobs += len(missing)

    summary = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "production_manifest": str(manifest_path),
        "submitted": args.submit,
        "jobs_per_task": args.jobs_per_task,
        "job_flavour": args.job_flavour,
        "total_jobs": total_jobs,
        "total_missing_jobs": total_missing_jobs,
        "deferred_jobs": total_missing_jobs - total_jobs,
        "tasks": submissions,
    }
    summary_path = (
        manifest_path.parent
        / ("submission_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + ".json")
    )
    with open(summary_path, "w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    action = "Submitted" if args.submit else "Prepared"
    print(
        f"{action} {total_jobs} jobs from {len(submissions)} tasks; "
        f"{total_missing_jobs - total_jobs} deferred."
    )
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
