#!/usr/bin/env python3

"""Build validated Run-2 Combine input trees without a DNN score."""

import argparse
import glob
import json
import os
import time
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


PROCESS_DATASETS = {
    "DY": [
        "DYJetsToLL",
        "DYJetsToLL_M50_HT100to200",
        "DYJetsToLL_M50_HT200to400",
        "DYJetsToLL_M50_HT400to600",
        "DYJetsToLL_M50_HT600to800",
        "DYJetsToLL_M50_HT800to1200",
        "DYJetsToLL_M50_HT1200to2500",
        "DYJetsToLL_M50_HT2500toInf",
    ],
    "Other": ["TTWJetsToLNu", "TTZToLLNuNu_M-10", "tZq_ll"],
    # Keep the historical grouping unchanged.  In particular, it contains
    # both the 5F and 4F t-channel samples and only the top tW sample.
    "ST": [
        "ST_s-channel",
        "ST_t-channel_antitop",
        "ST_t-channel_antitop_4f",
        "ST_t-channel_top",
        "ST_t-channel_top_4f",
        "ST_tW_top",
    ],
    "VVV": ["WWW", "WWZ", "WZZ", "ZZZ"],
    "WW": [
        "WWTo1L1Nu2Q",
        "WWTo2L2Nu",
        "GluGluToWWToENEN",
        "GluGluToWWToMNMN",
    ],
    "WZ": ["WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu"],
    "ggZZ": [
        "GluGluToContinToZZTo2e2nu",
        "GluGluToContinToZZTo2mu2nu",
    ],
    "qqZZ": ["ZZTo2L2Nu", "ZZTo2Q2L", "ZZTo4L"],
    "ttbar": ["TTTo2L2Nu"],
}

DATA_DATASETS = {
    "2016HIPM": [
        "DoubleEG", "DoubleMuon", "SingleElectron", "SingleMuon", "MuonEG"
    ],
    "2016noHIPM": [
        "DoubleEG", "DoubleMuon", "SingleElectron", "SingleMuon", "MuonEG"
    ],
    "2017": [
        "DoubleEG", "DoubleMuon", "SingleElectron", "SingleMuon", "MuonEG"
    ],
    "2018": ["EGamma", "DoubleMuon", "SingleMuon", "MuonEG"],
}

SIGNAL_DATASETS = {
    "2016HIPM": "ZH_ZToLL_HToInvisible_M125_UL16APV_full",
    "2016noHIPM": "ZH_ZToLL_HToInvisible_M125_UL16_full",
    "2017": "ZH_ZToLL_HToInvisible_M125_UL17_full",
    "2018": "ZH_ZToLL_HToInvisible_M125_UL18_full",
}

REGIONS = {
    "SR": ("dilepton", "lepton_cat != 2 && ptmiss >= 100"),
    "DY_CR": (
        "dilepton",
        "lepton_cat != 2 && ptmiss >= 30 && ptmiss <= 90",
    ),
    "emu_CR": ("dilepton", "lepton_cat == 2 && ptmiss >= 100"),
    "3l_CR": ("trilepton", "1"),
}

PERIODS = ["2016HIPM", "2016noHIPM", "2017", "2018"]
FINAL_YEARS = ["2016", "2017", "2018"]


def open_root(path, attempts=5):
    for attempt in range(attempts):
        root_file = ROOT.TFile.Open(str(path))
        if root_file and not root_file.IsZombie():
            return root_file
        if root_file:
            root_file.Close()
        if attempt + 1 < attempts:
            time.sleep(2)
    raise OSError(f"Could not open {path} after {attempts} attempts")


def files_for_datasets(
    output_dir, dataset_names, period=None, task=None, source_overrides=None
):
    files = []
    source_overrides = source_overrides or {}
    for dataset in dataset_names:
        dataset_dir = source_overrides.get(
            (period, task, dataset), output_dir
        )
        matches = sorted(
            glob.glob(os.path.join(dataset_dir, f"{dataset}_[0-9]*.root"))
        )
        if not matches:
            raise RuntimeError(
                f"No job outputs found for dataset {dataset} in {dataset_dir}"
            )
        files.extend(matches)
    return files


def parse_source_overrides(specifications):
    """Parse PERIOD:TASK:DATASET=OUTPUT_DIR source substitutions."""
    overrides = {}
    for specification in specifications:
        key, separator, directory = specification.partition("=")
        fields = key.split(":")
        if not separator or len(fields) != 3 or not directory:
            raise ValueError(
                "Invalid --source-override. Expected "
                "PERIOD:TASK:DATASET=OUTPUT_DIR, got "
                f"{specification!r}"
            )
        period, task, dataset = fields
        if period not in PERIODS:
            raise ValueError(f"Unknown override period {period!r}")
        if task not in {"dilepton", "trilepton", "signal"}:
            raise ValueError(f"Unknown override task {task!r}")
        override_key = (period, task, dataset)
        if override_key in overrides:
            raise ValueError(f"Duplicate source override for {override_key}")
        overrides[override_key] = directory
    return overrides


def inspect_sources(paths):
    expected_entries = 0
    schema = None
    for path in paths:
        root_file = open_root(path)
        if root_file.TestBit(ROOT.TFile.kRecovered):
            root_file.Close()
            raise OSError(f"Recovered ROOT file is not accepted: {path}")
        tree = root_file.Get("Vars")
        if not tree or not tree.InheritsFrom("TTree"):
            root_file.Close()
            raise OSError(f"Vars tree is missing in {path}")
        current_schema = tuple(
            (branch.GetName(), branch.GetClassName(), branch.GetTitle())
            for branch in tree.GetListOfBranches()
        )
        if schema is None:
            schema = current_schema
        elif current_schema != schema:
            root_file.Close()
            raise RuntimeError(f"Branch schema mismatch in {path}")
        expected_entries += int(tree.GetEntries())
        root_file.Close()
    return expected_entries


def copy_tree(paths, selection, output_path, overwrite):
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"{output_path} exists; pass --overwrite to replace it"
        )

    expected_entries = inspect_sources(paths)
    chain = ROOT.TChain("Vars")
    for path in paths:
        if chain.Add(path) != 1:
            raise OSError(f"Could not add {path} to TChain")
    num_input = int(chain.GetEntries())
    if chain.GetNtrees() != len(paths):
        raise RuntimeError(
            f"Expected {len(paths)} trees, TChain has {chain.GetNtrees()}"
        )
    if num_input != expected_entries:
        raise RuntimeError(
            f"TChain has {num_input} entries, expected {expected_entries}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    output = ROOT.TFile(str(temporary), "RECREATE")
    if not output or output.IsZombie():
        raise OSError(f"Could not create {temporary}")
    output.cd()
    source_file = None
    if num_input == 0:
        # TChain::CopyTree can segfault for an empty tree read through XRootD
        # (observed with ROOT 6.40).  Clone the schema from a concrete source
        # tree instead.  Any selection over zero entries has the same result.
        source_file = open_root(paths[0])
        selected = source_file.Get("Vars").CloneTree(0)
        output.cd()
        selected.SetDirectory(output)
    else:
        selected = chain.CopyTree(selection)
    if not selected:
        if source_file:
            source_file.Close()
        output.Close()
        raise RuntimeError(f"CopyTree failed for {output_path}")
    selected.Write("Vars")
    num_output = int(selected.GetEntries())
    output.Close()
    if source_file:
        source_file.Close()

    check_file = open_root(temporary)
    if check_file.TestBit(ROOT.TFile.kRecovered):
        check_file.Close()
        raise OSError(f"Recovered output ROOT file: {temporary}")
    check_tree = check_file.Get("Vars")
    check_entries = int(check_tree.GetEntries()) if check_tree else -1
    check_file.Close()
    if check_entries != num_output:
        raise RuntimeError(
            f"Output {temporary} has {check_entries} entries, "
            f"expected {num_output}"
        )
    os.replace(temporary, output_path)
    print(
        f"Wrote {output_path}: {num_output} / {num_input} entries "
        f"from {len(paths)} files",
        flush=True,
    )
    return {
        "path": str(output_path),
        "source_files": len(paths),
        "input_entries": num_input,
        "output_entries": num_output,
        "selection": selection,
    }


def build_period(
    period, input_base, staging_base, overwrite, report, source_overrides
):
    period_dir = Path(staging_base) / period
    for region, (task, selection) in REGIONS.items():
        source_dir = Path(input_base) / period / task / "output"
        for process, datasets in PROCESS_DATASETS.items():
            paths = files_for_datasets(
                str(source_dir), datasets, period, task, source_overrides
            )
            key = f"staging/{period}/{region}/{process}"
            report[key] = copy_tree(
                paths,
                selection,
                period_dir / region / f"{process}.root",
                overwrite,
            )

        data_paths = files_for_datasets(
            str(source_dir), DATA_DATASETS[period], period, task,
            source_overrides,
        )
        key = f"staging/{period}/{region}/Data"
        report[key] = copy_tree(
            data_paths,
            selection,
            period_dir / region / "Data.root",
            overwrite,
        )

    signal_dir = Path(input_base) / period / "signal" / "output"
    signal_paths = files_for_datasets(
        str(signal_dir), [SIGNAL_DATASETS[period]], period, "signal",
        source_overrides,
    )
    key = f"staging/{period}/SR/signal"
    report[key] = copy_tree(
        signal_paths,
        REGIONS["SR"][1],
        period_dir / "SR" / "signal.root",
        overwrite,
    )


def publish(staging_base, output_base, overwrite, report):
    filenames = [f"{name}.root" for name in PROCESS_DATASETS] + ["Data.root"]
    for year in FINAL_YEARS:
        periods = ["2016HIPM", "2016noHIPM"] if year == "2016" else [year]
        for region in REGIONS:
            for filename in filenames:
                paths = [str(Path(staging_base) / p / region / filename) for p in periods]
                key = f"final/{year}/{region}/{filename[:-5]}"
                report[key] = copy_tree(
                    paths,
                    "1",
                    Path(output_base) / year / region / filename,
                    overwrite,
                )

        signal_paths = [
            str(Path(staging_base) / p / "SR" / "signal.root")
            for p in periods
        ]
        key = f"final/{year}/SR/signal"
        report[key] = copy_tree(
            signal_paths,
            "1",
            Path(output_base) / year / "SR" / "signal.root",
            overwrite,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-base",
        default="/eos/user/l/liwe/monoz_full_dnn/production",
    )
    parser.add_argument(
        "--staging-base",
        default="/eos/user/l/liwe/monoz_full_dnn/staging/combine_input",
    )
    parser.add_argument(
        "--output-base", default="/eos/user/l/liwe/monoz_full_dnn"
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--publish-only",
        action="store_true",
        help="Publish final-year files from an existing staging directory.",
    )
    parser.add_argument(
        "--source-override",
        action="append",
        default=[],
        metavar="PERIOD:TASK:DATASET=OUTPUT_DIR",
        help=(
            "Use a replacement output directory for one dataset. May be "
            "specified more than once."
        ),
    )
    parser.add_argument(
        "--report",
        default=(
            "/eos/user/l/liwe/monoz_full_dnn/metadata/"
            "combine_build_report.json"
        ),
    )
    args = parser.parse_args()

    source_overrides = parse_source_overrides(args.source_override)
    report = {
        "_metadata": {
            "source_overrides": {
                ":".join(key): value
                for key, value in sorted(source_overrides.items())
            }
        }
    }
    if not args.publish_only:
        for period in PERIODS:
            build_period(
                period,
                args.input_base,
                args.staging_base,
                args.overwrite,
                report,
                source_overrides,
            )
    publish(
        args.staging_base,
        args.output_base,
        args.overwrite,
        report,
    )

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Wrote report {report_path}", flush=True)


if __name__ == "__main__":
    main()
