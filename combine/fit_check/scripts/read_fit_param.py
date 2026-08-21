#!/usr/bin/env python3

import argparse
import sys

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


def get_fit_result(root_file):
    for name in ["fit_s", "fit_b"]:
        result = root_file.Get(name)
        if result:
            return result
    raise KeyError("Could not find fit_s or fit_b in fitDiagnostics file")


def parse_args():
    parser = argparse.ArgumentParser(description="Print one fitted parameter value from fitDiagnostics.")
    parser.add_argument("--fit", required=True)
    parser.add_argument("--param", required=True)
    parser.add_argument("--field", choices=["value", "error"], default="value")
    return parser.parse_args()


def main():
    args = parse_args()
    root_file = ROOT.TFile.Open(args.fit)
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open {args.fit}")
    result = get_fit_result(root_file)
    params = result.floatParsFinal()
    param = params.find(args.param)
    if not param:
        print(f"Parameter {args.param} not found in {args.fit}", file=sys.stderr)
        return 2
    if args.field == "value":
        print(f"{param.getVal():.10g}")
    else:
        print(f"{param.getError():.10g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
