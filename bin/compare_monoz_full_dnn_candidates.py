#!/usr/bin/env python3

"""Compare two Combine-input directory trees branch by branch."""

import argparse
import json
from pathlib import Path

import awkward as ak
import uproot


def root_paths(base):
    return {
        str(path.relative_to(base)): path
        for path in Path(base).glob("20??/*/*.root")
    }


def compare_trees(reference_path, candidate_path, step_size):
    with uproot.open(f"{reference_path}:Vars") as reference, uproot.open(
        f"{candidate_path}:Vars"
    ) as candidate:
        reference_schema = reference.typenames()
        candidate_schema = candidate.typenames()
        result = {
            "reference_entries": int(reference.num_entries),
            "candidate_entries": int(candidate.num_entries),
            "schema_equal": reference_schema == candidate_schema,
            "common_prefix_equal": False,
        }
        if reference_schema != candidate_schema:
            result["reference_only_branches"] = sorted(
                set(reference_schema) - set(candidate_schema)
            )
            result["candidate_only_branches"] = sorted(
                set(candidate_schema) - set(reference_schema)
            )
            return result

        common_entries = min(reference.num_entries, candidate.num_entries)
        common_equal = True
        differing_branches = set()
        branches = list(reference_schema)
        for start in range(0, common_entries, step_size):
            stop = min(start + step_size, common_entries)
            reference_arrays = reference.arrays(
                branches, entry_start=start, entry_stop=stop, library="ak"
            )
            candidate_arrays = candidate.arrays(
                branches, entry_start=start, entry_stop=stop, library="ak"
            )
            if ak.almost_equal(
                reference_arrays, candidate_arrays, rtol=0, atol=0
            ):
                continue
            common_equal = False
            for branch in branches:
                if not ak.almost_equal(
                    reference_arrays[branch], candidate_arrays[branch],
                    rtol=0, atol=0,
                ):
                    differing_branches.add(branch)
        result["common_prefix_equal"] = common_equal
        if differing_branches:
            result["differing_branches"] = sorted(differing_branches)

        if candidate.num_entries > common_entries:
            identifiers = [
                name for name in ("run", "lumi", "event")
                if name in candidate_schema
            ]
            extras = candidate.arrays(
                identifiers,
                entry_start=common_entries,
                entry_stop=candidate.num_entries,
                library="np",
            )
            result["candidate_extra_event_ids"] = [
                {name: int(extras[name][index]) for name in identifiers}
                for index in range(candidate.num_entries - common_entries)
            ]
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--output", required=True)
    parser.add_argument("--step-size", type=int, default=100000)
    args = parser.parse_args()

    reference_paths = root_paths(args.reference)
    candidate_paths = root_paths(args.candidate)
    report = {
        "reference_base": args.reference,
        "candidate_base": args.candidate,
        "reference_files": len(reference_paths),
        "candidate_files": len(candidate_paths),
        "reference_only_files": sorted(
            set(reference_paths) - set(candidate_paths)
        ),
        "candidate_only_files": sorted(
            set(candidate_paths) - set(reference_paths)
        ),
        "files": {},
    }
    for relative_path in sorted(set(reference_paths) & set(candidate_paths)):
        comparison = compare_trees(
            reference_paths[relative_path], candidate_paths[relative_path],
            args.step_size,
        )
        if not (
            comparison["reference_entries"]
            == comparison["candidate_entries"]
            and comparison["schema_equal"]
            and comparison["common_prefix_equal"]
        ):
            report["files"][relative_path] = comparison

    report["different_files"] = len(report["files"])
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
