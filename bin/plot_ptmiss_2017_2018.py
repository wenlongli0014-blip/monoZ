#!/usr/bin/env python3

"""Plot combined 2017+2018 ptmiss distributions in the mono-Z regions."""

import argparse
import os
import sys
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hzz import Hist1D, mpl_style
from plot_data_sim import (
    DecoratedSample,
    Sample,
    Selection,
    Variable,
    plot_data_sim,
)


REGIONS = {
    "3l_CR": {
        "label": r"$3\ell$ CR",
        "range": [30.0, 500.0],
        "bins": 24,
        "x_scale": "linear",
    },
    "DY_CR": {
        "label": "DY CR",
        "range": [30.0, 90.0],
        "bins": 20,
        "x_scale": "linear",
    },
    "emu_CR": {
        "label": r"$e\mu$ CR",
        "range": [90.0, 1200.0],
        "bins": 20,
        "x_scale": "log",
    },
    "SR": {
        "label": "SR",
        "range": [100.0, 1200.0],
        "bins": 20,
        "x_scale": "log",
        "blind": True,
    },
}

BACKGROUND_SAMPLES = [
    ("VVV", "VVV", "mediumpurple"),
    (r"$gg \to ZZ$", "ggZZ", "crimson"),
    (r"$q\bar q \to ZZ$", "qqZZ", "darkorange"),
    ("WW", "WW", "mediumseagreen"),
    ("WZ", "WZ", "lightseagreen"),
    (r"$t\bar t$", "ttbar", "greenyellow"),
    (r"Single $t$, $tW$", "ST", "mediumblue"),
    ("Drell-Yan", "DY", "steelblue"),
    ("Other", "Other", "slateblue"),
]


def sample_paths(input_dir, region, sample):
    """Return the extension-free paths for both data-taking years."""

    return [
        os.path.join(input_dir, year, region, sample)
        for year in ("2017", "2018")
    ]


def build_configuration(input_dir, region):
    """Build the configuration expected by the standard plotting helpers."""

    region_info = REGIONS[region]
    selection = Selection({
        "formula": "1",
        "weight": "weight",
        "tag": region,
        "label": region_info["label"],
        "blind": region_info.get("blind", False),
    })
    variable = Variable({
        "formula": "ptmiss",
        "tag": "ptmiss",
        "label": r"$p_\mathrm{T}^\mathrm{miss}$",
        "unit": "GeV",
        "y_scale": "log",
        "x_scale": region_info["x_scale"],
        "binning": {
            "default": {
                "range": region_info["range"],
                "bins": region_info["bins"],
            }
        },
    })

    data_sample = Sample(
        {"files": sample_paths(input_dir, region, "Data")},
        tree_name="Vars",
    )
    sim_samples = [
        DecoratedSample({
            "label": label,
            "tag": tag,
            "files": sample_paths(input_dir, region, tag),
            "color": color,
        }, tree_name="Vars")
        for label, tag, color in BACKGROUND_SAMPLES
    ]

    signal_samples = []
    if region == "SR":
        signal_samples.append(DecoratedSample({
            "label": "signal",
            "tag": "signal",
            "files": sample_paths(input_dir, region, "signal"),
            "color": "red",
            "linestyle": "--",
            "linewidth": 1.8,
            "scale": 1.0,
        }, tree_name="Vars"))

    return SimpleNamespace(
        selections=[selection],
        variables=[variable],
        data_sample=data_sample,
        sim_samples=sim_samples,
        signal_samples=signal_samples,
        entries_label="Events",
    )


def fill_histogram(sample, binning, weighted):
    """Fill a Hist1D from all files in a sample using uproot."""

    contents = np.zeros(len(binning) + 1, dtype=np.float64)
    sumw2 = np.zeros_like(contents)

    for path in sample.files:
        branches = ["ptmiss", "weight"] if weighted else ["ptmiss"]
        with uproot.open(path) as root_file:
            arrays = root_file[sample.tree_name].arrays(
                branches, library="np")

        values = np.asarray(arrays["ptmiss"], dtype=np.float64)
        weights = (
            np.asarray(arrays["weight"], dtype=np.float64)
            if weighted else np.ones_like(values)
        )
        finite = np.isfinite(values) & np.isfinite(weights)
        values = values[finite]
        weights = weights[finite]

        underflow = values < binning[0]
        overflow = values >= binning[-1]
        in_range = ~(underflow | overflow)

        contents[0] += weights[underflow].sum()
        sumw2[0] += np.square(weights[underflow]).sum()
        contents[-1] += weights[overflow].sum()
        sumw2[-1] += np.square(weights[overflow]).sum()
        contents[1:-1] += np.histogram(
            values[in_range], bins=binning, weights=weights[in_range])[0]
        sumw2[1:-1] += np.histogram(
            values[in_range], bins=binning,
            weights=np.square(weights[in_range]))[0]

    # Match HistogramBuilder._tidy_hist: fold tails into visible edge bins.
    contents[1] += contents[0]
    contents[-2] += contents[-1]
    sumw2[1] += sumw2[0]
    sumw2[-2] += sumw2[-1]
    contents[0] = contents[-1] = 0.0
    sumw2[0] = sumw2[-1] = 0.0
    np.clip(contents, 0.0, None, out=contents)

    return Hist1D(
        binning=binning,
        contents=contents,
        errors=np.sqrt(sumw2),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="/eos/user/l/liwe/monoZ_combine_17_18_dnn",
        help="Directory containing the 2017 and 2018 region directories.",
    )
    parser.add_argument(
        "--output",
        default="/afs/cern.ch/user/l/liwe/hzz2l2nu/plot",
        help="Output directory.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Plot formats to write.",
    )
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    plt.style.use(mpl_style)

    for region, region_info in REGIONS.items():
        config = build_configuration(args.input, region)

        selection = config.selections[0]
        variable = config.variables[0]
        binning = variable.binning(region)
        if selection.blind:
            # The plotting helper replaces this placeholder with the summed
            # background expectation.  Do not read the real SR data at all.
            data_hist = Hist1D(binning=binning)
        else:
            data_hist = fill_histogram(
                config.data_sample, binning, weighted=False)
        sim_hists = [
            fill_histogram(sample, binning, weighted=True)
            for sample in config.sim_samples
        ]
        signal_hists = [
            fill_histogram(sample, binning, weighted=True)
            for sample in config.signal_samples
        ]
        output_stem = os.path.join(args.output, f"{region}_ptmiss_2017_2018")
        plot_data_sim(
            variable,
            data_hist,
            list(zip(sim_hists, config.sim_samples)),
            selection,
            output_stem,
            signal_hists_infos=list(zip(
                signal_hists,
                config.signal_samples,
            )),
            formats=args.formats,
            entries_label=config.entries_label,
            info_label=f"{region_info['label']}, 2017+2018",
            root_path=f"{output_stem}.root",
        )
        print(f"Wrote {output_stem}")


if __name__ == "__main__":
    main()
