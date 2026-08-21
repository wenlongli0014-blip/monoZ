#!/usr/bin/env python3

"""Rewrite the stale 2018 EGamma DDF with the current per-era paths."""

import argparse
import os
import re

import yaml


ERA_DIRECTORIES = {
    0: "EGamma_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1",
    1: "EGamma_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1",
    2: "EGamma_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1",
    3: "EGamma_Run2018D-UL2018_MiniAODv2_NanoAODv9-v3",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument(
        "--base",
        default="/eos/cms/store/group/phys_smp/ZZTo2L2Nu/HZZsample/2018",
    )
    args = parser.parse_args()

    with open(args.input) as stream:
        dataset = yaml.safe_load(stream)

    repaired = []
    for old_path in dataset["files"]:
        name = os.path.basename(old_path)
        match = re.match(r"EGamma_p([0-3])_", name)
        if not match:
            raise RuntimeError(f"Cannot infer Run era from {name}")

        era = int(match.group(1))
        new_path = os.path.join(args.base, ERA_DIRECTORIES[era], name)
        if not os.path.isfile(new_path):
            raise FileNotFoundError(new_path)
        repaired.append(new_path)

    if len(repaired) != len(set(repaired)):
        raise RuntimeError("The repaired DDF contains duplicate ROOT paths.")

    dataset["files"] = repaired
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as stream:
        yaml.safe_dump(dataset, stream, sort_keys=False)

    print(f"Wrote {len(repaired)} paths to {args.output}")


if __name__ == "__main__":
    main()
