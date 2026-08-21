#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import ROOT


def parse_args():
    parser = argparse.ArgumentParser(
        description="Copy a Combine workspace and import a generated Asimov toy as asimovData_1."
    )
    parser.add_argument("--input-workspace", required=True)
    parser.add_argument("--toy-file", required=True)
    parser.add_argument("--output-workspace", required=True)
    parser.add_argument("--workspace-name", default="w")
    parser.add_argument("--toy-name", default="toys/toy_asimov")
    parser.add_argument("--dataset-name", default="asimovData_1")
    return parser.parse_args()


def main():
    args = parse_args()
    input_workspace = Path(args.input_workspace)
    toy_file = Path(args.toy_file)
    output_workspace = Path(args.output_workspace)

    if not input_workspace.exists():
        raise FileNotFoundError(input_workspace)
    if not toy_file.exists():
        raise FileNotFoundError(toy_file)

    output_workspace.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_workspace, output_workspace)

    toy_handle = ROOT.TFile.Open(str(toy_file))
    if not toy_handle or toy_handle.IsZombie():
        raise OSError(f"Could not open toy file: {toy_file}")
    toy = toy_handle.Get(args.toy_name)
    if not toy:
        raise KeyError(f"Could not find {args.toy_name} in {toy_file}")
    toy.SetName(args.dataset_name)

    out_handle = ROOT.TFile.Open(str(output_workspace), "UPDATE")
    if not out_handle or out_handle.IsZombie():
        raise OSError(f"Could not update workspace file: {output_workspace}")
    workspace = out_handle.Get(args.workspace_name)
    if not workspace:
        raise KeyError(f"Could not find workspace {args.workspace_name} in {output_workspace}")
    if workspace.data(args.dataset_name):
        raise RuntimeError(f"{args.dataset_name} already exists in {output_workspace}")

    getattr(workspace, "import")(toy)
    out_handle.cd()
    workspace.Write(args.workspace_name, ROOT.TObject.kOverwrite)
    out_handle.Close()
    toy_handle.Close()

    print(f"Imported {args.toy_name} from {toy_file}")
    print(f"Wrote {args.dataset_name} into {output_workspace}")


if __name__ == "__main__":
    main()
