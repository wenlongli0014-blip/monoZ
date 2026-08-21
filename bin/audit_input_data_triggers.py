#!/usr/bin/env python3

"""Audit trigger decisions in the upstream dilepton data skims."""

import argparse
import json
from pathlib import Path

import numpy as np
import uproot
import yaml


DATASETS = {
    "2017": {
        "DoubleEG": "ee",
        "SingleElectron": "ee",
        "DoubleMuon": "mumu",
        "SingleMuon": "mumu",
        "MuonEG": "emu",
    },
    "2018": {
        "EGamma": "ee",
        "DoubleMuon": "mumu",
        "SingleMuon": "mumu",
        "MuonEG": "emu",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--years", nargs="+", choices=sorted(DATASETS), default=sorted(DATASETS)
    )
    parser.add_argument("--output", help="Optional output JSON path")
    return parser.parse_args()


def load_catalog(year, dataset):
    override = (
        Path("config")
        / "hadronic_tau_test"
        / f"ddf_{year}"
        / f"{dataset}.yaml"
    )
    central = Path(
        f"/eos/cms/store/group/phys_smp/ZZTo2L2Nu/HZZsample/"
        f"{year}/YAML/{dataset}.yaml"
    )
    catalog = override if override.exists() else central
    node = yaml.safe_load(catalog.read_text())
    return catalog, node["files"]


def trigger_names(blocks):
    return sorted(
        {
            f"HLT_{name}"
            for block in blocks
            for name in block["triggers"]
        }
    )


def audit_file(path, blocks, configured_triggers):
    tree = uproot.open(f"{path}:Events")
    available = [name for name in configured_triggers if name in tree]
    expressions = ["run", *available]
    arrays = tree.arrays(expressions, library="np")
    run = arrays["run"]
    decision = np.zeros(len(run), dtype=bool)

    for block in blocks:
        block_triggers = [
            f"HLT_{name}"
            for name in block["triggers"]
            if f"HLT_{name}" in arrays
        ]
        if not block_triggers:
            continue
        block_decision = np.logical_or.reduce(
            [arrays[name] for name in block_triggers]
        )
        if "run_range" in block:
            first_run, last_run = block["run_range"]
            block_decision &= (run >= first_run) & (run <= last_run)
        decision |= block_decision

    return {
        "entries": int(len(run)),
        "passed": int(np.count_nonzero(decision)),
        "failed": int(np.count_nonzero(~decision)),
        "available_triggers": available,
        "missing_triggers": sorted(set(configured_triggers) - set(available)),
    }


def main():
    args = parse_args()
    result = {"years": {}}

    for year in args.years:
        config = yaml.safe_load(
            Path(f"config/dilepton/{year}-ul.yaml").read_text()
        )
        year_result = {}
        for dataset, channel in DATASETS[year].items():
            catalog, files = load_catalog(year, dataset)
            blocks = config["trigger_filter"][channel]
            configured_triggers = trigger_names(blocks)
            dataset_result = {
                "channel": channel,
                "catalog": str(catalog),
                "num_files": len(files),
                "files": [],
            }
            for path in files:
                file_result = audit_file(path, blocks, configured_triggers)
                file_result["path"] = path
                dataset_result["files"].append(file_result)

            dataset_result["entries"] = sum(
                item["entries"] for item in dataset_result["files"]
            )
            dataset_result["passed"] = sum(
                item["passed"] for item in dataset_result["files"]
            )
            dataset_result["failed"] = sum(
                item["failed"] for item in dataset_result["files"]
            )
            year_result[dataset] = dataset_result
        result["years"][year] = year_result

    encoded = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded + "\n")
    print(encoded)


if __name__ == "__main__":
    main()
