#!/usr/bin/env python3

"""Validate a Combine FitDiagnostics result and write a compact JSON summary."""

import argparse
import json
import math
from pathlib import Path

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


def get_fit_result(root_file):
    for name in ("fit_s", "fit_b"):
        result = root_file.Get(name)
        if result:
            return name, result
    raise KeyError("Could not find fit_s or fit_b in FitDiagnostics output")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", required=True, type=Path)
    parser.add_argument("--params", nargs="+", required=True)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    root_file = ROOT.TFile.Open(str(args.fit))
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {args.fit}")

    fit_name, result = get_fit_result(root_file)
    final_parameters = result.floatParsFinal()
    parameters = {}
    problems = []
    for name in args.params:
        parameter = final_parameters.find(name)
        if not parameter:
            problems.append(f"missing parameter {name}")
            continue
        value = float(parameter.getVal())
        error = float(parameter.getError())
        parameters[name] = {
            "value": value,
            "error": error,
            "error_lo": float(parameter.getAsymErrorLo()),
            "error_hi": float(parameter.getAsymErrorHi()),
        }
        if not math.isfinite(value) or not math.isfinite(error):
            problems.append(f"non-finite result for {name}")

    summary = {
        "fit_file": str(args.fit),
        "fit_result": fit_name,
        "status": int(result.status()),
        "covQual": int(result.covQual()),
        "edm": float(result.edm()),
        "minNll": float(result.minNll()),
        "parameters": parameters,
        "problems": problems,
    }
    if summary["status"] != 0:
        problems.append(f"fit status is {summary['status']}, expected 0")
    if summary["covQual"] < 3:
        problems.append(
            f"covariance quality is {summary['covQual']}, expected 3"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")
    root_file.Close()

    print(json.dumps(summary, indent=2))
    if problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
