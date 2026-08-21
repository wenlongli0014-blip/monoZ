#!/usr/bin/env python3

"""Rewrite DDF paths for ROOT files relocated below their dataset directory."""

import argparse
from collections import defaultdict
import os

import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument(
        "--replace",
        action="append",
        default=[],
        metavar="OLD_BASENAME=NEW_PATH",
        help="Explicit replacement when basename lookup is not sufficient.",
    )
    parser.add_argument(
        "--path-replace",
        action="append",
        default=[],
        metavar="OLD=NEW",
        help="Try a literal path substitution before basename lookup.",
    )
    args = parser.parse_args()

    explicit = defaultdict(list)
    for value in args.replace:
        old_name, separator, new_path = value.partition("=")
        if not separator:
            parser.error("--replace must have the form OLD_BASENAME=NEW_PATH")
        explicit[old_name].append(new_path)

    path_replacements = []
    for value in args.path_replace:
        old, separator, new = value.partition("=")
        if not separator:
            parser.error("--path-replace must have the form OLD=NEW")
        path_replacements.append((old, new))

    with open(args.input) as stream:
        dataset = yaml.safe_load(stream)

    repaired = []
    for old_path in dataset["files"]:
        if os.path.isfile(old_path):
            repaired.append(old_path)
            continue

        old_name = os.path.basename(old_path)
        if old_name in explicit:
            candidates = explicit[old_name]
        else:
            candidates = [
                old_path.replace(old, new)
                for old, new in path_replacements
                if old in old_path
            ]
            candidates = [path for path in candidates if os.path.isfile(path)]
            if not candidates:
                dataset_dir = os.path.dirname(old_path)
                for root, _, files in os.walk(dataset_dir):
                    if old_name in files:
                        candidates.append(os.path.join(root, old_name))

        candidates = [path for path in candidates if os.path.isfile(path)]
        if not candidates:
            raise RuntimeError(
                f"No replacement found for {old_path}"
            )
        repaired.extend(candidates)

    if len(repaired) != len(set(repaired)):
        raise RuntimeError("The repaired DDF contains duplicate ROOT paths.")

    dataset["files"] = repaired
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as stream:
        yaml.safe_dump(dataset, stream, sort_keys=False)

    print(f"Wrote {len(repaired)} paths to {args.output}")


if __name__ == "__main__":
    main()
