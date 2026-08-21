#!/usr/bin/env python3

"""Build the four grouped Combine input regions for the tau-veto study."""

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
    "ST": [
        "ST_s-channel",
        "ST_t-channel_antitop",
        "ST_t-channel_antitop_4f",
        "ST_t-channel_top",
        "ST_t-channel_top_4f",
        "ST_tW_top",
    ],
    "VVV": ["WWW", "WWZ", "WZZ", "ZZZ"],
    "WW": ["WWTo1L1Nu2Q", "WWTo2L2Nu", "GluGluToWWToENEN", "GluGluToWWToMNMN"],
    "WZ": ["WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu"],
    "ggZZ": ["GluGluToContinToZZTo2e2nu", "GluGluToContinToZZTo2mu2nu"],
    "qqZZ": ["ZZTo2L2Nu", "ZZTo2Q2L", "ZZTo4L"],
    "ttbar": ["TTTo2L2Nu"],
}

DATA_DATASETS = {
    "2017": ["DoubleEG", "DoubleMuon", "SingleElectron", "SingleMuon", "MuonEG"],
    "2018": ["EGamma", "DoubleMuon", "SingleMuon", "MuonEG"],
}

REGIONS = {
    "SR": ("dilepton_jobs", "lepton_cat != 2 && ptmiss >= 100"),
    "DY_CR": (
        "dilepton_jobs",
        "lepton_cat != 2 && ptmiss >= 30 && ptmiss <= 90",
    ),
    "emu_CR": ("dilepton_jobs", "lepton_cat == 2 && ptmiss >= 100"),
    "3l_CR": ("trilepton_jobs", "1"),
}


def open_root(path, attempts=5):
    for attempt in range(attempts):
        root_file = ROOT.TFile.Open(path)
        if root_file and not root_file.IsZombie():
            return root_file
        if root_file:
            root_file.Close()
        if attempt + 1 < attempts:
            time.sleep(2)
    raise OSError(f"Could not open {path} after {attempts} attempts")


def files_for_datasets(output_dir, dataset_names):
    files = []
    for dataset in dataset_names:
        matches = sorted(
            glob.glob(os.path.join(output_dir, f"{dataset}_[0-9]*.root"))
        )
        if not matches:
            raise RuntimeError(
                f"No job outputs found for dataset {dataset} in {output_dir}"
            )
        files.extend(matches)
    return files


def copy_tree(paths, selection, output_path, overwrite):
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"{output_path} exists; pass --overwrite to replace it."
        )

    expected_entries = 0
    for path in paths:
        root_file = open_root(path)
        if root_file.TestBit(ROOT.TFile.kRecovered):
            root_file.Close()
            raise OSError(f"Recovered ROOT file is not accepted: {path}")
        tree = root_file.Get("Vars")
        if not tree or not tree.InheritsFrom("TTree"):
            root_file.Close()
            raise OSError(f"Vars tree is missing in {path}")
        expected_entries += int(tree.GetEntries())
        root_file.Close()

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
            f"TChain for {output_path} has {num_input} entries, "
            f"but pre-opened files have {expected_entries}"
        )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    output = ROOT.TFile(output_path, "RECREATE")
    if not output or output.IsZombie():
        raise OSError(f"Could not create {output_path}")
    output.cd()
    selected = chain.CopyTree(selection)
    if not selected:
        output.Close()
        raise RuntimeError(f"CopyTree failed for {output_path}")
    selected.Write("Vars")
    num_output = int(selected.GetEntries())
    output.Close()

    check_file = open_root(output_path)
    if check_file.TestBit(ROOT.TFile.kRecovered):
        check_file.Close()
        raise OSError(f"Recovered output ROOT file: {output_path}")
    check_tree = check_file.Get("Vars")
    check_entries = (
        int(check_tree.GetEntries()) if check_tree else -1
    )
    check_file.Close()
    if check_entries != num_output:
        raise RuntimeError(
            f"Output {output_path} has {check_entries} entries after reopen, "
            f"expected {num_output}"
        )
    print(
        f"Wrote {output_path}: {num_output} / {num_input} entries "
        f"from {len(paths)} files"
    )


def ensure_complete(tasks, variant, year):
    selected = {
        task["task"]: task
        for task in tasks
        if task["variant"] == variant and task["year"] == year
    }
    required = {"dilepton_jobs", "trilepton_jobs", "signal_jobs"}
    if set(selected) != required:
        raise RuntimeError(
            f"Manifest tasks for {variant} {year}: {sorted(selected)}"
        )

    for task in selected.values():
        with open(Path(task["task_dir"]) / "job_names.dat") as stream:
            names = [line.strip() for line in stream if line.strip()]
        missing = [
            name for name in names
            if not os.path.isfile(Path(task["output_dir"]) / f"{name}.root")
        ]
        if missing:
            raise RuntimeError(
                f"{variant} {year} {task['task']} has {len(missing)} "
                f"missing outputs; first is {missing[0]}"
            )
    return selected


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
        "--variants",
        nargs="+",
        default=["tau_medium", "tau_vloose", "tau_vvvloose", "tau_none"],
    )
    parser.add_argument("--years", nargs="+", default=["2017", "2018"])
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--signal-only",
        action="store_true",
        help="Only rebuild SR/signal.root for the requested variants/years.",
    )
    args = parser.parse_args()

    with open(args.manifest) as stream:
        manifest = json.load(stream)
    study_base = Path(args.manifest).parent.parent

    for variant in args.variants:
        for year in args.years:
            tasks = ensure_complete(manifest["tasks"], variant, year)
            output_dirs = {
                name: task["output_dir"] for name, task in tasks.items()
            }
            region_base = study_base / variant / "combine_input" / year

            if not args.signal_only:
                for region, (task_name, selection) in REGIONS.items():
                    source_dir = output_dirs[task_name]
                    for process, datasets in PROCESS_DATASETS.items():
                        paths = files_for_datasets(source_dir, datasets)
                        copy_tree(
                            paths,
                            selection,
                            str(region_base / region / f"{process}.root"),
                            args.overwrite,
                        )

                    data_paths = files_for_datasets(
                        source_dir, DATA_DATASETS[year]
                    )
                    copy_tree(
                        data_paths,
                        selection,
                        str(region_base / region / "Data.root"),
                        args.overwrite,
                    )

            signal_dir = output_dirs["signal_jobs"]
            signal_dataset = (
                f"ZH_ZToLL_HToInvisible_M125_UL{year[-2:]}_shears_skim"
            )
            signal_paths = files_for_datasets(signal_dir, [signal_dataset])
            copy_tree(
                signal_paths,
                REGIONS["SR"][1],
                str(region_base / "SR" / "signal.root"),
                args.overwrite,
            )


if __name__ == "__main__":
    main()
