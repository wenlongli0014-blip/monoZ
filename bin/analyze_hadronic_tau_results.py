#!/usr/bin/env python3

"""Build compact physics and closure summaries for the tau-veto study."""

import argparse
import csv
import ctypes
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


VARIANTS = ["tau_medium", "tau_vloose", "tau_vvvloose", "tau_none"]
YEARS = ["2017", "2018"]
CHANNELS = ["SR", "DYCR", "EMUCR", "CR3L"]
BACKGROUND_PROCESSES = [
    "DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"
]
PROCESSES = [
    *BACKGROUND_PROCESSES,
    "data_obs",
    "signal",
]
FIXED_CANDIDATE_FIELDS = [
    "tau_preselected_candidates",
    "tau_vvvloose_candidates",
    "tau_vloose_candidates",
    "tau_medium_candidates",
]
FOCUS_DATASETS = {
    "WZ": {"WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu"},
    "ttbar": {"TTTo2L2Nu"},
    "signal": {
        "ZH_ZToLL_HToInvisible_M125_UL17_shears_skim",
        "ZH_ZToLL_HToInvisible_M125_UL18_shears_skim",
    },
}


def close_enough(left, right, rel_tol=1e-7, abs_tol=1e-6):
    return math.isclose(left, right, rel_tol=rel_tol, abs_tol=abs_tol)


def write_csv(path, rows):
    if not rows:
        raise RuntimeError(f"No rows to write to {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit_closure(audit_rows):
    problems = []
    grouped = defaultdict(dict)
    for row in audit_rows:
        key = (row["year"], row["task"], row["dataset"])
        grouped[key][row["variant"]] = row

    for key, rows in grouped.items():
        if set(rows) != set(VARIANTS):
            problems.append({
                "kind": "missing_audit_variant",
                "key": key,
                "variants": sorted(rows),
            })
            continue

        reference = rows["tau_none"]
        for variant, row in rows.items():
            if row["audit_entries"] != reference["audit_entries"]:
                problems.append({
                    "kind": "audit_entry_mismatch",
                    "key": key,
                    "variant": variant,
                })
            if not close_enough(
                row["audit_sum_weight"], reference["audit_sum_weight"]
            ):
                problems.append({
                    "kind": "audit_weight_mismatch",
                    "key": key,
                    "variant": variant,
                })
            for field in FIXED_CANDIDATE_FIELDS:
                if row[field] != reference[field]:
                    problems.append({
                        "kind": f"{field}_mismatch",
                        "key": key,
                        "variant": variant,
                    })

            if row["tau_genuine_candidates"] is not None:
                flavour_sum = sum(
                    row[field]
                    for field in [
                        "tau_genuine_candidates",
                        "tau_electron_candidates",
                        "tau_muon_candidates",
                        "tau_unmatched_candidates",
                    ]
                )
                if flavour_sum != row["tau_selected_candidates"]:
                    problems.append({
                        "kind": "gen_flavour_closure",
                        "key": key,
                        "variant": variant,
                        "selected": row["tau_selected_candidates"],
                        "flavour_sum": flavour_sum,
                    })

        expected_selected = {
            "tau_medium": reference["tau_medium_candidates"],
            "tau_none": reference["tau_medium_candidates"],
            "tau_vloose": reference["tau_vloose_candidates"],
            "tau_vvvloose": reference["tau_vvvloose_candidates"],
        }
        for variant, expected in expected_selected.items():
            if rows[variant]["tau_selected_candidates"] != expected:
                problems.append({
                    "kind": "selected_candidate_wp_mismatch",
                    "key": key,
                    "variant": variant,
                })

        rejected = {
            variant: rows[variant]["tau_rejected_events"]
            for variant in VARIANTS
        }
        if not (
            rejected["tau_none"] == 0
            and rejected["tau_medium"] <= rejected["tau_vloose"]
            and rejected["tau_vloose"] <= rejected["tau_vvvloose"]
        ):
            problems.append({
                "kind": "audit_rejection_nesting",
                "key": key,
                "values": rejected,
            })

        entries = {
            variant: rows[variant]["vars_entries"] for variant in VARIANTS
        }
        if not (
            entries["tau_vvvloose"] <= entries["tau_vloose"]
            and entries["tau_vloose"] <= entries["tau_medium"]
            and entries["tau_medium"] <= entries["tau_none"]
        ):
            problems.append({
                "kind": "dataset_vars_nesting",
                "key": key,
                "values": entries,
            })
    return len(grouped), problems


def region_closure(region_report):
    problems = []
    grouped = defaultdict(dict)
    for row in region_report["files"]:
        key = (row["year"], row["region"], row["process"])
        grouped[key][row["variant"]] = row["entries"]

    for key, entries in grouped.items():
        if set(entries) != set(VARIANTS):
            problems.append({
                "kind": "missing_region_variant",
                "key": key,
                "variants": sorted(entries),
            })
            continue
        if not (
            entries["tau_vvvloose"] <= entries["tau_vloose"]
            and entries["tau_vloose"] <= entries["tau_medium"]
            and entries["tau_medium"] <= entries["tau_none"]
        ):
            problems.append({
                "kind": "region_entry_nesting",
                "key": key,
                "values": entries,
            })
    return len(grouped), problems


def aggregate_focus_audit(audit_rows):
    aggregates = {}
    fields = [
        "audit_entries",
        "audit_sum_weight",
        "tau_rejected_events",
        "tau_rejected_sum_weight",
        "tau_preselected_candidates",
        "tau_selected_candidates",
        "tau_genuine_candidates",
        "tau_electron_candidates",
        "tau_muon_candidates",
        "tau_unmatched_candidates",
        "vars_entries",
        "vars_sum_weight",
    ]
    for row in audit_rows:
        for process, datasets in FOCUS_DATASETS.items():
            if row["dataset"] not in datasets:
                continue
            key = (row["variant"], row["task"], process)
            if key not in aggregates:
                aggregates[key] = {
                    "variant": row["variant"],
                    "task": row["task"],
                    "process": process,
                    **{field: 0 for field in fields},
                }
            target = aggregates[key]
            for field in fields:
                value = row[field]
                if value is not None:
                    target[field] += value

    rows = []
    for key in sorted(aggregates):
        row = aggregates[key]
        row["audit_rejected_fraction"] = (
            row["tau_rejected_events"] / row["audit_entries"]
            if row["audit_entries"] else 0.0
        )
        row["selected_genuine_fraction"] = (
            row["tau_genuine_candidates"] / row["tau_selected_candidates"]
            if row["tau_selected_candidates"] else 0.0
        )
        row["selected_unmatched_fraction"] = (
            row["tau_unmatched_candidates"] / row["tau_selected_candidates"]
            if row["tau_selected_candidates"] else 0.0
        )
        rows.append(row)
    return rows


def read_shape_yields(study_base):
    rows = []
    for variant in VARIANTS:
        path = (
            study_base / variant / "combine"
            / "shapes_ptmiss_2017_2018_statonly.root"
        )
        root_file = ROOT.TFile.Open(str(path))
        if not root_file or root_file.IsZombie():
            raise OSError(f"Could not open {path}")
        for channel in CHANNELS:
            for process in PROCESSES:
                hist = root_file.Get(f"{channel}/{process}")
                if not hist:
                    raise KeyError(f"Missing {channel}/{process} in {path}")
                error = ctypes.c_double(0.0)
                value = hist.IntegralAndError(
                    1, hist.GetNbinsX(), error
                )
                rows.append({
                    "variant": variant,
                    "channel": channel,
                    "process": process,
                    "yield": float(value),
                    "stat_error": error.value,
                    "raw_entries": float(hist.GetEntries()),
                })
        root_file.Close()

    lookup = {
        (row["variant"], row["channel"], row["process"]): row
        for row in rows
    }
    for row in rows:
        baseline = lookup[
            ("tau_none", row["channel"], row["process"])
        ]["yield"]
        medium = lookup[
            ("tau_medium", row["channel"], row["process"])
        ]["yield"]
        row["retention_vs_none"] = (
            row["yield"] / baseline if baseline else None
        )
        row["ratio_vs_medium"] = (
            row["yield"] / medium if medium else None
        )
    return rows


def read_yearly_sr_yields(study_base):
    rows = []
    for variant in VARIANTS:
        for year in YEARS:
            for process in [*BACKGROUND_PROCESSES, "signal"]:
                path = (
                    study_base / variant / "combine_input" / year / "SR"
                    / f"{process}.root"
                )
                frame = ROOT.RDataFrame("Vars", str(path))
                num_entries = int(frame.Count().GetValue())
                finite = frame.Filter("std::isfinite(weight)").Define(
                    "tau_study_weight2", "weight * weight"
                )
                num_finite = int(finite.Count().GetValue())
                rows.append({
                    "variant": variant,
                    "year": year,
                    "process": process,
                    "raw_entries": num_entries,
                    "nonfinite_weights": num_entries - num_finite,
                    "yield": float(finite.Sum("weight").GetValue()),
                    "sum_weight2": float(
                        finite.Sum("tau_study_weight2").GetValue()
                    ),
                })

    lookup = {
        (row["variant"], row["year"], row["process"]): row
        for row in rows
    }
    for row in rows:
        baseline = lookup[
            ("tau_none", row["year"], row["process"])
        ]["yield"]
        medium = lookup[
            ("tau_medium", row["year"], row["process"])
        ]["yield"]
        row["retention_vs_none"] = (
            row["yield"] / baseline if baseline else None
        )
        row["ratio_vs_medium"] = (
            row["yield"] / medium if medium else None
        )
    return rows


def read_sr_bin_diagnostics(study_base):
    rows = []
    for variant in VARIANTS:
        path = (
            study_base / variant / "combine"
            / "shapes_ptmiss_2017_2018_statonly.root"
        )
        root_file = ROOT.TFile.Open(str(path))
        if not root_file or root_file.IsZombie():
            raise OSError(f"Could not open {path}")
        signal = root_file.Get("SR/signal")
        backgrounds = {
            process: root_file.Get(f"SR/{process}")
            for process in BACKGROUND_PROCESSES
        }
        for process, hist in backgrounds.items():
            if not hist:
                raise KeyError(f"Missing SR/{process} in {path}")
        if not signal:
            raise KeyError(f"Missing SR/signal in {path}")

        for bin_index in range(1, signal.GetNbinsX() + 1):
            values = {
                process: hist.GetBinContent(bin_index)
                for process, hist in backgrounds.items()
            }
            total_background = sum(values.values())
            total_error = math.sqrt(sum(
                hist.GetBinError(bin_index) ** 2
                for hist in backgrounds.values()
            ))
            signal_yield = signal.GetBinContent(bin_index)
            floored = [
                process for process, value in values.items()
                if 0 < value <= 1.000001e-6
            ]
            floor_total = sum(values[process] for process in floored)
            rows.append({
                "variant": variant,
                "bin": bin_index,
                "ptmiss_low": signal.GetXaxis().GetBinLowEdge(bin_index),
                "ptmiss_high": signal.GetXaxis().GetBinUpEdge(bin_index),
                "signal": signal_yield,
                "background": total_background,
                "background_stat_error": total_error,
                "relative_background_stat_error": (
                    total_error / total_background
                    if total_background else None
                ),
                "s_over_b": (
                    signal_yield / total_background
                    if total_background else None
                ),
                "s_over_sqrt_b": (
                    signal_yield / math.sqrt(total_background)
                    if total_background > 0 else None
                ),
                "dy": values["DY"],
                "wz": values["WZ"],
                "ttbar": values["ttbar"],
                "num_floored_backgrounds": len(floored),
                "floored_backgrounds": ",".join(floored),
                "floor_yield": floor_total,
                "floor_fraction": (
                    floor_total / total_background
                    if total_background else None
                ),
            })
        root_file.Close()

    lookup = {
        (row["variant"], row["bin"]): row for row in rows
    }
    for row in rows:
        baseline = lookup[("tau_none", row["bin"])]
        medium = lookup[("tau_medium", row["bin"])]
        row["signal_retention_vs_none"] = (
            row["signal"] / baseline["signal"]
            if baseline["signal"] else None
        )
        row["background_retention_vs_none"] = (
            row["background"] / baseline["background"]
            if baseline["background"] else None
        )
        row["s_over_b_ratio_vs_none"] = (
            row["s_over_b"] / baseline["s_over_b"]
            if baseline["s_over_b"] else None
        )
        row["signal_ratio_vs_medium"] = (
            row["signal"] / medium["signal"]
            if medium["signal"] else None
        )
        row["background_ratio_vs_medium"] = (
            row["background"] / medium["background"]
            if medium["background"] else None
        )
    return rows


def read_limits(study_base):
    expected_pattern = re.compile(
        r"Expected\s+([0-9.]+)%:\s+r <\s+([0-9.]+)"
    )
    observed_pattern = re.compile(r"Observed Limit:\s+r <\s+([0-9.]+)")
    rows = []
    for variant in VARIANTS:
        path = (
            study_base / variant / "combine" / "logs"
            / "10_limit_expected.log"
        )
        text = path.read_text()
        expected = {
            percentile: float(value)
            for percentile, value in expected_pattern.findall(text)
        }
        observed_match = observed_pattern.search(text)
        if set(expected) != {"2.5", "16.0", "50.0", "84.0", "97.5"}:
            raise RuntimeError(f"Expected quantiles are incomplete in {path}")
        if not observed_match:
            raise RuntimeError(f"Observed Asimov limit is missing in {path}")
        rows.append({
            "variant": variant,
            "expected_2p5": expected["2.5"],
            "expected_16": expected["16.0"],
            "expected_50": expected["50.0"],
            "expected_84": expected["84.0"],
            "expected_97p5": expected["97.5"],
            # Combine labels this line "Observed" even for an explicit -t -1
            # Asimov run. It is not a limit from real SR data.
            "asimov_limit": float(observed_match.group(1)),
        })

    none_limit = next(
        row["expected_50"] for row in rows if row["variant"] == "tau_none"
    )
    medium_limit = next(
        row["expected_50"]
        for row in rows if row["variant"] == "tau_medium"
    )
    for row in rows:
        row["improvement_vs_none"] = 1.0 - row["expected_50"] / none_limit
        row["improvement_vs_medium"] = (
            1.0 - row["expected_50"] / medium_limit
        )
    return rows


def read_fit_summaries(study_base):
    summaries = {}
    for variant in VARIANTS:
        output_base = study_base / variant / "combine" / "fits"
        with (output_base / "cr_fit_summary.json").open() as stream:
            cr_fit = json.load(stream)
        with (
            output_base / "expected_closure_summary.json"
        ).open() as stream:
            closure = json.load(stream)
        summaries[variant] = {
            "cr_fit": cr_fit,
            "expected_closure": closure,
        }
    return summaries


def compact_physics_summary(shape_rows, limit_rows, fit_summaries):
    shape = {
        (row["variant"], row["channel"], row["process"]): row
        for row in shape_rows
    }
    limits = {row["variant"]: row for row in limit_rows}
    summary = {}
    for variant in VARIANTS:
        sr_background = sum(
            shape[(variant, "SR", process)]["yield"]
            for process in BACKGROUND_PROCESSES
        )
        summary[variant] = {
            "sr_background": sr_background,
            "sr_signal": shape[(variant, "SR", "signal")]["yield"],
            "sr_wz": shape[(variant, "SR", "WZ")]["yield"],
            "sr_ttbar": shape[(variant, "SR", "ttbar")]["yield"],
            "cr3l_wz": shape[(variant, "CR3L", "WZ")]["yield"],
            "expected_median_limit": limits[variant]["expected_50"],
            "limit_improvement_vs_none": (
                limits[variant]["improvement_vs_none"]
            ),
            "cr_fit": fit_summaries[variant]["cr_fit"],
            "expected_closure": (
                fit_summaries[variant]["expected_closure"]
            ),
        }

    none = summary["tau_none"]
    for values in summary.values():
        values["sr_background_retention_vs_none"] = (
            values["sr_background"] / none["sr_background"]
        )
        values["sr_signal_retention_vs_none"] = (
            values["sr_signal"] / none["sr_signal"]
        )
        values["sr_wz_retention_vs_none"] = (
            values["sr_wz"] / none["sr_wz"]
        )
        values["cr3l_wz_retention_vs_none"] = (
            values["cr3l_wz"] / none["cr3l_wz"]
        )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study-base",
        default="/eos/user/l/liwe/hadronic_tau_test",
        type=Path,
    )
    args = parser.parse_args()
    metadata = args.study_base / "metadata"

    with (metadata / "tau_audit_summary.json").open() as stream:
        audit_rows = json.load(stream)
    with (metadata / "region_integrity_report.json").open() as stream:
        region_report = json.load(stream)

    audit_groups, audit_problems = audit_closure(audit_rows)
    region_groups, region_problems = region_closure(region_report)
    focus_rows = aggregate_focus_audit(audit_rows)
    shape_rows = read_shape_yields(args.study_base)
    yearly_rows = read_yearly_sr_yields(args.study_base)
    sr_bin_rows = read_sr_bin_diagnostics(args.study_base)
    limit_rows = read_limits(args.study_base)
    fit_summaries = read_fit_summaries(args.study_base)

    closure = {
        "audit_groups_checked": audit_groups,
        "region_groups_checked": region_groups,
        "num_audit_problems": len(audit_problems),
        "num_region_problems": len(region_problems),
        "audit_problems": audit_problems,
        "region_problems": region_problems,
    }
    results = {
        "scope": {
            "years": YEARS,
            "observable": "ptmiss",
            "combine_model": (
                "three-CR-conditioned stat-only background-only Asimov"
            ),
            "sr_data": "blinded; SR data_obs is background-only Asimov",
            "rate_parameters": {
                "k_Zjet": ["DY"],
                "k_WZ": ["WZ"],
                "k_emu": ["ttbar", "WW", "ST"],
            },
            "autoMCStats": False,
            "variants": VARIANTS,
        },
        "closure": closure,
        "physics_summary": compact_physics_summary(
            shape_rows, limit_rows, fit_summaries
        ),
        "limits": limit_rows,
        "diagnostics": {
            "yearly_sr_nonfinite_weights": sum(
                row["nonfinite_weights"] for row in yearly_rows
            ),
            "max_sr_floor_fraction": max(
                row["floor_fraction"] for row in sr_bin_rows
                if row["floor_fraction"] is not None
            ),
        },
    }

    write_csv(metadata / "shape_yields.csv", shape_rows)
    write_csv(metadata / "yearly_sr_yields.csv", yearly_rows)
    write_csv(metadata / "sr_bin_diagnostics.csv", sr_bin_rows)
    write_csv(metadata / "combine_limits.csv", limit_rows)
    write_csv(metadata / "tau_audit_focus.csv", focus_rows)
    with (metadata / "study_results.json").open("w") as stream:
        json.dump(results, stream, indent=2)
        stream.write("\n")

    print(json.dumps(results, indent=2))
    print(f"Wrote result tables under {metadata}")
    if audit_problems or region_problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
