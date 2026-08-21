#!/usr/bin/env python3

"""Plot the WZ response to the hadronic-tau veto working point scan."""

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import uproot


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

VARIANTS = ["tau_none", "tau_medium", "tau_vloose", "tau_vvvloose"]
WP_VARIANTS = VARIANTS[1:]
LABELS = {
    "tau_none": "No veto",
    "tau_medium": "Medium (>=16)",
    "tau_vloose": "VLoose (>=4)",
    "tau_vvvloose": "VVVLoose (>=1)",
}
SHORT_LABELS = {
    "tau_none": "No veto",
    "tau_medium": "Medium",
    "tau_vloose": "VLoose",
    "tau_vvvloose": "VVVLoose",
}
COLORS = {
    "tau_none": "#303030",
    "tau_medium": "#2878B5",
    "tau_vloose": "#2B9A66",
    "tau_vvvloose": "#C53D3D",
}


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def as_float(row, field):
    return float(row[field])


def indexed(rows, *fields):
    return {tuple(row[field] for field in fields): row for row in rows}


def nested_ratio_error(numerator, numerator_error, denominator, denominator_error):
    """Error for A/B when selected events A are a subset of baseline B."""
    if denominator == 0:
        return 0.0
    variance_a = numerator_error ** 2
    variance_b = denominator_error ** 2
    covariance = variance_a
    variance = (
        variance_a / denominator ** 2
        + numerator ** 2 * variance_b / denominator ** 4
        - 2.0 * numerator * covariance / denominator ** 3
    )
    return math.sqrt(max(0.0, variance))


def independent_ratio_error(numerator, numerator_error, denominator, denominator_error):
    if numerator == 0 or denominator == 0:
        return 0.0
    ratio = numerator / denominator
    return abs(ratio) * math.sqrt(
        (numerator_error / numerator) ** 2
        + (denominator_error / denominator) ** 2
    )


def set_style():
    plt.rcParams.update({
        "figure.figsize": (7.4, 5.8),
        "figure.dpi": 130,
        "savefig.dpi": 180,
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 6.5,
    })


def cms_label(ax, qualifier="Simulation", right="2017+2018"):
    ax.text(
        0.0, 1.015, "CMS", transform=ax.transAxes, ha="left", va="bottom",
        fontsize=15, fontweight="bold",
    )
    ax.text(
        0.105, 1.015, qualifier, transform=ax.transAxes, ha="left", va="bottom",
        fontsize=11, style="italic",
    )
    ax.text(
        1.0, 1.015, right, transform=ax.transAxes, ha="right", va="bottom",
        fontsize=11,
    )


def save_figure(fig, output_dir, stem):
    fig.savefig(output_dir / (stem + ".png"), bbox_inches="tight")
    fig.savefig(output_dir / (stem + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_wz_retention(shape, output_dir, summary):
    x = np.arange(len(VARIANTS))
    fig, ax = plt.subplots()
    for channel, name, marker in [("SR", "SR", "o"), ("CR3L", "3l CR", "s")]:
        baseline = shape[("tau_none", channel, "WZ")]
        baseline_yield = as_float(baseline, "yield")
        baseline_error = as_float(baseline, "stat_error")
        values = []
        errors = []
        for variant in VARIANTS:
            row = shape[(variant, channel, "WZ")]
            value = as_float(row, "yield")
            error = as_float(row, "stat_error")
            values.append(100.0 * value / baseline_yield)
            errors.append(100.0 * nested_ratio_error(
                value, error, baseline_yield, baseline_error
            ))
        ax.errorbar(x, values, yerr=errors, marker=marker, capsize=3, label=name)
        summary["wz_retention"][channel] = dict(zip(VARIANTS, values))

    ax.axhline(100.0, color="0.65", linestyle="--", linewidth=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in VARIANTS])
    ax.set_ylabel("WZ retention relative to no veto [%]")
    ax.set_ylim(77.0, 103.0)
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    ax.legend(loc="lower left")
    cms_label(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "wz_retention_sr_cr3l")


def read_wz_histograms(study_base):
    result = {}
    for variant in VARIANTS:
        path = (
            study_base / variant / "combine"
            / "shapes_ptmiss_2017_2018_statonly.root"
        )
        root_file = ROOT.TFile.Open(str(path))
        if not root_file or root_file.IsZombie():
            raise OSError("Could not open {}".format(path))
        hist = root_file.Get("SR/WZ")
        if not hist:
            raise KeyError("Missing SR/WZ in {}".format(path))
        result[variant] = {
            "edges": [hist.GetXaxis().GetBinLowEdge(1)] + [
                hist.GetXaxis().GetBinUpEdge(i)
                for i in range(1, hist.GetNbinsX() + 1)
            ],
            "values": [
                hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)
            ],
            "errors": [
                hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1)
            ],
        }
        root_file.Close()
    return result


def plot_wz_ptmiss(histograms, output_dir, summary):
    edges = histograms["tau_none"]["edges"]
    bin_labels = [
        "{:.0f}-{:.0f}".format(edges[i], edges[i + 1])
        for i in range(len(edges) - 1)
    ]
    x = np.arange(len(bin_labels))
    fig, (ax, ratio_ax) = plt.subplots(
        2, 1, figsize=(7.4, 7.0), sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.2], "hspace": 0.04},
    )
    baseline = histograms["tau_none"]
    baseline_values = np.asarray(baseline["values"])
    baseline_errors = np.asarray(baseline["errors"])
    summary["wz_ptmiss"] = {}

    for variant in VARIANTS:
        values = np.asarray(histograms[variant]["values"])
        errors = np.asarray(histograms[variant]["errors"])
        ax.errorbar(
            x, values, yerr=errors, marker="o", capsize=2.5,
            color=COLORS[variant], label=LABELS[variant],
        )
        ratios = values / baseline_values
        ratio_errors = np.asarray([
            nested_ratio_error(a, ea, b, eb)
            for a, ea, b, eb in zip(values, errors, baseline_values, baseline_errors)
        ])
        ratio_ax.errorbar(
            x, ratios, yerr=ratio_errors, marker="o", capsize=2.5,
            color=COLORS[variant],
        )
        summary["wz_ptmiss"][variant] = {
            "yield": values.tolist(),
            "retention": ratios.tolist(),
        }

    ax.set_yscale("log")
    ax.set_ylabel("Pre-fit WZ yield")
    ax.set_ylim(0.8, 1000.0)
    ax.grid(axis="y", which="both", color="0.91", linewidth=0.8)
    ax.legend(ncol=2, loc="upper right")
    cms_label(ax)

    ratio_ax.axhline(1.0, color="0.55", linestyle="--", linewidth=1.2)
    ratio_ax.set_ylabel("Ratio to\nno veto")
    ratio_ax.set_xlabel(r"$p_{T}^{miss}$ bin [GeV]")
    ratio_ax.set_ylim(0.53, 1.05)
    ratio_ax.set_xticks(x)
    ratio_ax.set_xticklabels(bin_labels)
    ratio_ax.grid(axis="y", color="0.9", linewidth=0.8)
    fig.align_ylabels([ax, ratio_ax])
    save_figure(fig, output_dir, "wz_sr_ptmiss_shape_ratio")


def plot_transfer_factor(shape, output_dir, summary):
    x = np.arange(len(VARIANTS))
    factors = []
    errors = []
    for variant in VARIANTS:
        sr = shape[(variant, "SR", "WZ")]
        cr = shape[(variant, "CR3L", "WZ")]
        sr_yield = as_float(sr, "yield")
        cr_yield = as_float(cr, "yield")
        factors.append(sr_yield / cr_yield)
        errors.append(independent_ratio_error(
            sr_yield, as_float(sr, "stat_error"),
            cr_yield, as_float(cr, "stat_error"),
        ))
    baseline = factors[0]
    relative = [value / baseline for value in factors]
    summary["wz_transfer_factor"] = {
        variant: {"value": factors[i], "relative_to_none": relative[i]}
        for i, variant in enumerate(VARIANTS)
    }

    fig, (ax, ratio_ax) = plt.subplots(
        2, 1, figsize=(7.4, 6.8), sharex=True,
        gridspec_kw={"height_ratios": [2.7, 1.1], "hspace": 0.04},
    )
    ax.errorbar(x, factors, yerr=errors, marker="o", capsize=3, color="#2878B5")
    for index, value in enumerate(factors):
        ax.text(index, value + 0.0018, "{:.4f}".format(value), ha="center", fontsize=10)
    ax.set_ylabel(r"$T_{WZ}=N_{WZ}^{SR}/N_{WZ}^{3l\,CR}$")
    ax.set_ylim(0.098, 0.130)
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    cms_label(ax)

    ratio_ax.plot(x, relative, marker="o", color="#2878B5")
    ratio_ax.axhline(1.0, color="0.55", linestyle="--", linewidth=1.2)
    ratio_ax.set_ylabel("Ratio to\nno veto")
    ratio_ax.set_ylim(0.81, 1.03)
    ratio_ax.set_xticks(x)
    ratio_ax.set_xticklabels([LABELS[v] for v in VARIANTS])
    ratio_ax.grid(axis="y", color="0.9", linewidth=0.8)
    save_figure(fig, output_dir, "wz_sr_cr3l_transfer_factor")


def plot_signal_wz_tradeoff(shape, output_dir, summary):
    baseline_signal = as_float(shape[("tau_none", "SR", "signal")], "yield")
    baseline_wz = as_float(shape[("tau_none", "SR", "WZ")], "yield")
    fig, ax = plt.subplots()
    values = {}
    for variant in VARIANTS:
        signal = as_float(shape[(variant, "SR", "signal")], "yield")
        wz = as_float(shape[(variant, "SR", "WZ")], "yield")
        signal_loss = 100.0 * (1.0 - signal / baseline_signal)
        wz_rejection = 100.0 * (1.0 - wz / baseline_wz)
        values[variant] = {
            "signal_loss_percent": signal_loss,
            "wz_rejection_percent": wz_rejection,
        }
        ax.scatter(
            signal_loss, wz_rejection, s=75, color=COLORS[variant], zorder=3,
        )
        offset = (5, 7) if variant != "tau_none" else (6, 6)
        ax.annotate(
            SHORT_LABELS[variant], (signal_loss, wz_rejection),
            xytext=offset, textcoords="offset points", fontsize=10,
        )
    summary["signal_wz_tradeoff"] = values
    ax.plot(
        [values[v]["signal_loss_percent"] for v in VARIANTS],
        [values[v]["wz_rejection_percent"] for v in VARIANTS],
        color="0.7", linewidth=1.2, zorder=1,
    )
    ax.set_xlabel("Signal yield loss relative to no veto [%]")
    ax.set_ylabel("SR WZ rejection relative to no veto [%]")
    ax.set_xlim(-0.2, 3.7)
    ax.set_ylim(-1.0, 21.0)
    ax.grid(color="0.9", linewidth=0.8)
    cms_label(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "signal_efficiency_wz_rejection_tradeoff")


def plot_wz_prefit_postfit(shape, study_base, output_dir, summary):
    x = np.arange(len(VARIANTS))
    prefit = []
    prefit_errors = []
    postfit = []
    postfit_errors = []
    summary["wz_sr_prefit_postfit"] = {}
    for variant in VARIANTS:
        row = shape[(variant, "SR", "WZ")]
        value = as_float(row, "yield")
        error = as_float(row, "stat_error")
        fit_path = study_base / variant / "combine" / "fits" / "cr_fit_summary.json"
        with fit_path.open() as stream:
            fit = json.load(stream)
        k_wz = float(fit["parameters"]["k_WZ"]["value"])
        k_wz_error = float(fit["parameters"]["k_WZ"]["error"])
        fitted_value = k_wz * value
        fitted_error = math.sqrt((k_wz * error) ** 2 + (value * k_wz_error) ** 2)
        prefit.append(value)
        prefit_errors.append(error)
        postfit.append(fitted_value)
        postfit_errors.append(fitted_error)
        summary["wz_sr_prefit_postfit"][variant] = {
            "prefit_yield": value,
            "k_WZ": k_wz,
            "k_WZ_error": k_wz_error,
            "cr_normalized_yield": fitted_value,
        }

    fig, ax = plt.subplots()
    offset = 0.08
    ax.errorbar(
        x - offset, prefit, yerr=prefit_errors, marker="o", capsize=3,
        linestyle="none", color="#303030", label="Pre-fit",
    )
    ax.errorbar(
        x + offset, postfit, yerr=postfit_errors, marker="s", capsize=3,
        linestyle="none", color="#2878B5", label=r"CR-normalized ($k_{WZ}$)",
    )
    ax.plot(x, prefit, color="#303030", alpha=0.45, linewidth=1.3)
    ax.plot(x, postfit, color="#2878B5", alpha=0.45, linewidth=1.3)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in VARIANTS])
    ax.set_ylabel("SR WZ yield")
    ax.set_ylim(840.0, 1160.0)
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    ax.legend(loc="upper right")
    ax.text(
        0.02, 0.03, r"CR-normalized error: MC stat. $\oplus$ $k_{WZ}$ fit error",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=9,
    )
    cms_label(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "wz_sr_prefit_crnormalized_yield")


def plot_yearly_wz_retention(yearly_rows, output_dir, summary):
    data = indexed(yearly_rows, "variant", "year", "process")
    x = np.arange(len(VARIANTS))
    fig, ax = plt.subplots()
    summary["yearly_wz_retention"] = {}
    for year, marker in [("2017", "o"), ("2018", "s")]:
        baseline = data[("tau_none", year, "WZ")]
        baseline_yield = as_float(baseline, "yield")
        baseline_error = math.sqrt(as_float(baseline, "sum_weight2"))
        values = []
        errors = []
        for variant in VARIANTS:
            row = data[(variant, year, "WZ")]
            value = as_float(row, "yield")
            error = math.sqrt(as_float(row, "sum_weight2"))
            values.append(100.0 * value / baseline_yield)
            errors.append(100.0 * nested_ratio_error(
                value, error, baseline_yield, baseline_error
            ))
        summary["yearly_wz_retention"][year] = dict(zip(VARIANTS, values))
        ax.errorbar(x, values, yerr=errors, marker=marker, capsize=3, label=year)
    ax.axhline(100.0, color="0.65", linestyle="--", linewidth=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in VARIANTS])
    ax.set_ylabel("SR WZ retention relative to no veto [%]")
    ax.set_ylim(77.0, 103.0)
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    ax.legend(loc="lower left")
    cms_label(ax, right="2017 / 2018")
    fig.tight_layout()
    save_figure(fig, output_dir, "wz_sr_retention_by_year")


def plot_wz_candidate_composition(audit_rows, output_dir, summary):
    rows = indexed(audit_rows, "variant", "task", "process")
    fields = [
        ("tau_genuine_candidates", r"Genuine $\tau_h$", "#2878B5"),
        ("tau_electron_candidates", "Electron-related", "#E3B23C"),
        ("tau_muon_candidates", "Muon-related", "#C53D3D"),
        ("tau_unmatched_candidates", "Jet/unmatched", "#7A7A7A"),
    ]
    x = np.arange(len(WP_VARIANTS))
    bottoms = np.zeros(len(WP_VARIANTS))
    totals = np.asarray([
        as_float(rows[(variant, "dilepton_jobs", "WZ")], "tau_selected_candidates")
        for variant in WP_VARIANTS
    ])
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    summary["wz_candidate_composition"] = {}
    for field, label, color in fields:
        counts = np.asarray([
            as_float(rows[(variant, "dilepton_jobs", "WZ")], field)
            for variant in WP_VARIANTS
        ])
        fractions = 100.0 * counts / totals
        bars = ax.bar(x, fractions, bottom=bottoms, color=color, label=label, width=0.68)
        for index, (bar, fraction) in enumerate(zip(bars, fractions)):
            if fraction >= 3.0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bottoms[index] + fraction / 2.0,
                    "{:.1f}%".format(fraction), ha="center", va="center",
                    color="white" if color != "#E3B23C" else "black", fontsize=9,
                )
        for index, variant in enumerate(WP_VARIANTS):
            summary["wz_candidate_composition"].setdefault(variant, {})[field] = fractions[index]
        bottoms += fractions
    for index, total in enumerate(totals):
        ax.text(
            index, 101.5, r"$N_{\rm cand}$=" + "{:,.0f}".format(total),
            ha="center", fontsize=9,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in WP_VARIANTS])
    ax.set_ylabel("DileptonTrees WZ tau-candidate composition [%]")
    ax.set_ylim(0.0, 108.0)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.04), ncol=2)
    cms_label(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "wz_tau_candidate_genmatch_composition")


def read_weighted_wz_candidate_composition(study_base):
    datasets = ["WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu"]
    branches = {
        "genuine": "n_tau_genuine",
        "electron": "n_tau_electron",
        "muon": "n_tau_muon",
        "unmatched": "n_tau_unmatched",
        "selected": "n_tau_selected",
    }
    result = {}
    for variant in WP_VARIANTS:
        chain = ROOT.TChain("TauVetoAudit")
        paths = []
        for year in ["2017", "2018"]:
            audit_dir = study_base / variant / year / "dilepton_jobs" / "output"
            for dataset in datasets:
                paths.extend(sorted(audit_dir.glob(dataset + "_*.root")))
        if not paths:
            raise FileNotFoundError("No WZ audit files found for {}".format(variant))
        for path in paths:
            if chain.Add(str(path)) != 1:
                raise OSError("Could not add TauVetoAudit from {}".format(path))

        frame = ROOT.RDataFrame(chain)
        actions = {}
        for category, branch in branches.items():
            column = "weighted_tau_{}".format(category)
            frame = frame.Define(column, "weight * {}".format(branch))
            actions[category] = frame.Sum(column)
        values = {
            category: float(action.GetValue())
            for category, action in actions.items()
        }
        category_sum = sum(
            values[key] for key in ["genuine", "electron", "muon", "unmatched"]
        )
        if not math.isclose(
            category_sum, values["selected"], rel_tol=1e-10, abs_tol=1e-7
        ):
            raise RuntimeError(
                "Weighted candidate composition does not close for {}: {} != {}".format(
                    variant, category_sum, values["selected"]
                )
            )
        result[variant] = values
    return result


def plot_weighted_wz_candidate_composition(weighted, output_dir, summary):
    fields = [
        ("genuine", r"Genuine $\tau_h$", "#2878B5"),
        ("electron", "Electron-related", "#E3B23C"),
        ("muon", "Muon-related", "#C53D3D"),
        ("unmatched", "Jet/unmatched", "#7A7A7A"),
    ]
    x = np.arange(len(WP_VARIANTS))
    totals = np.asarray([
        weighted[variant]["selected"] for variant in WP_VARIANTS
    ])
    bottoms = np.zeros(len(WP_VARIANTS))
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    summary["wz_candidate_composition_weighted"] = {}
    for field, label, color in fields:
        values = np.asarray([
            weighted[variant][field] for variant in WP_VARIANTS
        ])
        if np.any(values < 0.0):
            raise RuntimeError(
                "Weighted WZ candidate category {} is negative".format(field)
            )
        fractions = 100.0 * values / totals
        bars = ax.bar(
            x, fractions, bottom=bottoms, color=color, label=label, width=0.68
        )
        for index, (bar, fraction) in enumerate(zip(bars, fractions)):
            if fraction >= 3.0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bottoms[index] + fraction / 2.0,
                    "{:.1f}%".format(fraction), ha="center", va="center",
                    color="white" if color != "#E3B23C" else "black", fontsize=9,
                )
        for index, variant in enumerate(WP_VARIANTS):
            summary["wz_candidate_composition_weighted"].setdefault(
                variant, {}
            )[field] = {
                "sum_weight": values[index],
                "fraction_percent": fractions[index],
            }
        bottoms += fractions
    for index, (variant, total) in enumerate(zip(WP_VARIANTS, totals)):
        summary["wz_candidate_composition_weighted"][variant][
            "selected_sum_weight"
        ] = total
        ax.text(
            index, 101.5,
            r"$\sum (w\,n_{\rm cand})$=" + "{:,.1f}".format(total),
            ha="center", fontsize=9,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in WP_VARIANTS])
    ax.set_ylabel("Weighted DileptonTrees WZ tau-candidate composition [%]")
    ax.set_ylim(0.0, 108.0)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.04), ncol=2)
    cms_label(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "wz_tau_candidate_genmatch_composition_weighted")


def event_keys(arrays):
    keys = np.empty(
        len(arrays["event"]),
        dtype=[("run", "<u4"), ("lumi", "<u4"), ("event", "<u8")],
    )
    keys["run"] = arrays["run"]
    keys["lumi"] = arrays["lumi"]
    keys["event"] = arrays["event"]
    return keys


def read_weighted_wz_candidates_in_sr(study_base):
    datasets = ["WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu"]
    categories = ["genuine", "electron", "muon", "unmatched", "selected"]
    branches = [
        "run", "lumi", "event", "weight",
        *["n_tau_{}".format(category) for category in categories],
    ]
    result = {}
    for variant in WP_VARIANTS:
        values = {category: 0.0 for category in categories}
        values["no_veto_sr_events"] = 0
        values["no_veto_sr_sum_weight"] = 0.0
        values["selected_event_count"] = 0
        values["selected_event_sum_weight"] = 0.0
        values["genuine_event_count"] = 0
        values["genuine_event_sum_weight"] = 0.0
        values["fake_only_event_count"] = 0
        values["fake_only_event_sum_weight"] = 0.0
        for category in ["electron", "muon", "unmatched"]:
            values[category + "_event_count"] = 0
            values[category + "_event_sum_weight"] = 0.0
        for category in ["electron_only", "muon_only", "unmatched_only", "mixed_fake"]:
            values[category + "_event_count"] = 0
            values[category + "_event_sum_weight"] = 0.0
        for year in ["2017", "2018"]:
            audit_dir = study_base / variant / year / "dilepton_jobs" / "output"
            none_dir = study_base / "tau_none" / year / "dilepton_jobs" / "output"
            for dataset in datasets:
                for audit_path in sorted(audit_dir.glob(dataset + "_*.root")):
                    none_path = none_dir / audit_path.name
                    if not none_path.is_file():
                        raise FileNotFoundError(none_path)
                    with uproot.open(none_path) as none_file:
                        final_arrays = none_file["Vars"].arrays(
                            [
                                "run", "lumi", "event", "weight",
                                "lepton_cat", "ptmiss",
                            ],
                            library="np",
                        )
                    sr_mask = (
                        (final_arrays["lepton_cat"] != 2)
                        & (final_arrays["ptmiss"] >= 100.0)
                    )
                    if not np.any(sr_mask):
                        continue
                    final_keys = event_keys(final_arrays)[sr_mask]
                    if len(np.unique(final_keys)) != len(final_keys):
                        raise RuntimeError(
                            "Duplicate no-veto SR event key in {}".format(none_path)
                        )
                    with uproot.open(audit_path) as audit_file:
                        audit_arrays = audit_file["TauVetoAudit"].arrays(
                            branches, library="np"
                        )
                    audit_keys = event_keys(audit_arrays)
                    mask = np.isin(audit_keys, final_keys)
                    if np.count_nonzero(mask) != len(final_keys):
                        raise RuntimeError(
                            "Could not match all {} final Vars events in {} TauVetoAudit".format(
                                len(final_keys), audit_path
                            )
                        )
                    weights = audit_arrays["weight"][mask]
                    for category in categories:
                        counts = audit_arrays["n_tau_{}".format(category)][mask]
                        values[category] += float(np.sum(weights * counts))
                    selected_mask = audit_arrays["n_tau_selected"][mask] > 0
                    genuine_mask = (
                        selected_mask
                        & (audit_arrays["n_tau_genuine"][mask] > 0)
                    )
                    fake_only_mask = selected_mask & ~genuine_mask
                    has_electron = audit_arrays["n_tau_electron"][mask] > 0
                    has_muon = audit_arrays["n_tau_muon"][mask] > 0
                    has_unmatched = audit_arrays["n_tau_unmatched"][mask] > 0
                    event_origin_masks = {
                        "electron": fake_only_mask & has_electron,
                        "muon": (
                            fake_only_mask & ~has_electron & has_muon
                        ),
                        "unmatched": (
                            fake_only_mask & ~has_electron & ~has_muon
                            & has_unmatched
                        ),
                    }
                    num_fake_origins = (
                        has_electron.astype(np.int8)
                        + has_muon.astype(np.int8)
                        + has_unmatched.astype(np.int8)
                    )
                    fake_category_masks = {
                        "electron_only": (
                            fake_only_mask & has_electron
                            & (num_fake_origins == 1)
                        ),
                        "muon_only": (
                            fake_only_mask & has_muon
                            & (num_fake_origins == 1)
                        ),
                        "unmatched_only": (
                            fake_only_mask & has_unmatched
                            & (num_fake_origins == 1)
                        ),
                        "mixed_fake": fake_only_mask & (num_fake_origins > 1),
                    }
                    values["selected_event_count"] += int(
                        np.count_nonzero(selected_mask)
                    )
                    values["selected_event_sum_weight"] += float(
                        np.sum(weights[selected_mask])
                    )
                    values["genuine_event_count"] += int(
                        np.count_nonzero(genuine_mask)
                    )
                    values["genuine_event_sum_weight"] += float(
                        np.sum(weights[genuine_mask])
                    )
                    values["fake_only_event_count"] += int(
                        np.count_nonzero(fake_only_mask)
                    )
                    values["fake_only_event_sum_weight"] += float(
                        np.sum(weights[fake_only_mask])
                    )
                    for category, category_mask in event_origin_masks.items():
                        values[category + "_event_count"] += int(
                            np.count_nonzero(category_mask)
                        )
                        values[category + "_event_sum_weight"] += float(
                            np.sum(weights[category_mask])
                        )
                    for category, category_mask in fake_category_masks.items():
                        values[category + "_event_count"] += int(
                            np.count_nonzero(category_mask)
                        )
                        values[category + "_event_sum_weight"] += float(
                            np.sum(weights[category_mask])
                        )
                    values["no_veto_sr_events"] += len(final_keys)
                    values["no_veto_sr_sum_weight"] += float(
                        np.sum(final_arrays["weight"][sr_mask])
                    )

        category_sum = sum(
            values[key] for key in ["genuine", "electron", "muon", "unmatched"]
        )
        if not math.isclose(
            category_sum, values["selected"], rel_tol=1e-10, abs_tol=1e-7
        ):
            raise RuntimeError(
                "SR weighted candidate composition does not close for {}".format(
                    variant
                )
            )
        event_category_sum = (
            values["genuine_event_sum_weight"]
            + values["fake_only_event_sum_weight"]
        )
        if not math.isclose(
            event_category_sum,
            values["selected_event_sum_weight"],
            rel_tol=1e-10,
            abs_tol=1e-5,
        ):
            raise RuntimeError(
                "SR weighted event composition does not close for {}".format(
                    variant
                )
            )
        fake_category_sum = sum(
            values[category + "_event_sum_weight"]
            for category in [
                "electron_only", "muon_only", "unmatched_only", "mixed_fake"
            ]
        )
        if not math.isclose(
            fake_category_sum,
            values["fake_only_event_sum_weight"],
            rel_tol=1e-10,
            abs_tol=1e-5,
        ):
            raise RuntimeError(
                "SR fake-only event composition does not close for {}".format(
                    variant
                )
            )
        four_class_sum = (
            values["genuine_event_sum_weight"]
            + values["electron_event_sum_weight"]
            + values["muon_event_sum_weight"]
            + values["unmatched_event_sum_weight"]
        )
        if not math.isclose(
            four_class_sum,
            values["selected_event_sum_weight"],
            rel_tol=1e-10,
            abs_tol=1e-5,
        ):
            raise RuntimeError(
                "SR four-class event composition does not close for {}".format(
                    variant
                )
            )
        result[variant] = values
    return result


def plot_weighted_wz_veto_events_in_sr(weighted, output_dir, summary):
    fields = [
        (
            "genuine_event_sum_weight",
            r"Contains genuine $\tau_h$",
            "#2878B5",
        ),
        (
            "electron_event_sum_weight",
            "Electron-related",
            "#E3B23C",
        ),
        (
            "muon_event_sum_weight",
            "Muon-related",
            "#C53D3D",
        ),
        (
            "unmatched_event_sum_weight",
            "Jet/unmatched",
            "#7A7A7A",
        ),
    ]
    x = np.arange(len(WP_VARIANTS))
    totals = np.asarray([
        weighted[variant]["selected_event_sum_weight"]
        for variant in WP_VARIANTS
    ])
    bottoms = np.zeros(len(WP_VARIANTS))
    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    key = "wz_tau_veto_event_composition_weighted_sr"
    summary[key] = {}
    for field, label, color in fields:
        values = np.asarray([
            weighted[variant][field] for variant in WP_VARIANTS
        ])
        if np.any(values < 0.0):
            raise RuntimeError(
                "SR veto-event category {} is negative".format(field)
            )
        fractions = 100.0 * values / totals
        bars = ax.bar(
            x, fractions, bottom=bottoms, color=color, label=label, width=0.68
        )
        for index, (bar, fraction) in enumerate(zip(bars, fractions)):
            if fraction >= 3.0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bottoms[index] + fraction / 2.0,
                    "{:.1f}%".format(fraction), ha="center", va="center",
                    color="black" if color == "#E3B23C" else "white",
                    fontsize=9,
                )
        for index, variant in enumerate(WP_VARIANTS):
            summary[key].setdefault(variant, {})[field] = {
                "sum_weight": values[index],
                "fraction_percent": fractions[index],
            }
        bottoms += fractions
    for index, (variant, total) in enumerate(zip(WP_VARIANTS, totals)):
        summary[key][variant]["veto_event_sum_weight"] = total
        summary[key][variant]["no_veto_sr_events"] = weighted[variant][
            "no_veto_sr_events"
        ]
        summary[key][variant]["no_veto_sr_sum_weight"] = weighted[variant][
            "no_veto_sr_sum_weight"
        ]
        summary[key][variant]["selected_event_count"] = weighted[variant][
            "selected_event_count"
        ]
        summary[key][variant]["selected_event_sum_weight"] = weighted[variant][
            "selected_event_sum_weight"
        ]
        summary[key][variant]["genuine_event_count"] = weighted[variant][
            "genuine_event_count"
        ]
        summary[key][variant]["fake_only_event_count"] = weighted[variant][
            "fake_only_event_count"
        ]
        summary[key][variant]["fake_only_event_sum_weight"] = weighted[variant][
            "fake_only_event_sum_weight"
        ]
        for category in ["electron", "muon", "unmatched"]:
            summary[key][variant][category + "_event_count"] = weighted[variant][
                category + "_event_count"
            ]
        for category in [
            "electron_only", "muon_only", "unmatched_only", "mixed_fake"
        ]:
            summary[key][variant][category + "_event_count"] = weighted[variant][
                category + "_event_count"
            ]
        ax.text(
            index, 101.5,
            r"$\sum w_{\rm events}$=" + "{:,.2f}".format(total),
            ha="center", fontsize=9,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in WP_VARIANTS])
    ax.set_xlabel(
        r"No-veto SR: lepton category $\ne e\mu$, "
        r"$p_T^{miss} \geq 100$ GeV"
    )
    ax.set_ylabel("Weighted SR WZ tau-veto event composition [%]")
    ax.set_ylim(0.0, 108.0)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.04), ncol=2)
    cms_label(ax)
    fig.tight_layout()
    save_figure(
        fig, output_dir,
        "wz_tau_veto_event_genmatch_composition_weighted_sr",
    )


def plot_wz_subsample_rejection(audit_summary_rows, output_dir, summary):
    datasets = ["WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu"]
    labels = {
        "WZTo1L1Nu2Q": "WZTo1L1Nu2Q",
        "WZTo1L3Nu": "WZTo1L3Nu",
        "WZTo2Q2L": "WZTo2Q2L",
        "WZTo3LNu": "WZTo3LNu",
    }
    markers = ["o", "s", "^", "D"]
    colors = ["#2878B5", "#E07A24", "#2B9A66", "#C53D3D"]
    x = np.arange(len(WP_VARIANTS))
    fig, ax = plt.subplots()
    summary["wz_subsample_tau_stage_rejection"] = {}
    for dataset, marker, color in zip(datasets, markers, colors):
        fractions = []
        errors = []
        for variant in WP_VARIANTS:
            selected_rows = [
                row for row in audit_summary_rows
                if row["variant"] == variant
                and row["task"] == "dilepton_jobs"
                and row["dataset"] == dataset
            ]
            entries = sum(as_float(row, "audit_entries") for row in selected_rows)
            rejected = sum(as_float(row, "tau_rejected_events") for row in selected_rows)
            fraction = rejected / entries if entries else 0.0
            fractions.append(100.0 * fraction)
            errors.append(
                100.0 * math.sqrt(fraction * (1.0 - fraction) / entries)
                if entries else 0.0
            )
        summary["wz_subsample_tau_stage_rejection"][dataset] = dict(
            zip(WP_VARIANTS, fractions)
        )
        ax.errorbar(
            x, fractions, yerr=errors, marker=marker, color=color,
            capsize=3, label=labels[dataset],
        )
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[v] for v in WP_VARIANTS])
    ax.set_ylabel("Events rejected at tau-veto stage [%]")
    ax.set_ylim(-0.5, 22.0)
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    ax.legend(loc="upper left", ncol=2)
    ax.text(
        0.99, 0.03, "DileptonTrees audit, unweighted events",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
    )
    cms_label(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "wz_subsample_tau_stage_rejection")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study-base", type=Path,
        default=Path("/eos/user/l/liwe/hadronic_tau_test"),
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("/eos/user/l/liwe/hadronic_tau_test/plot"),
    )
    args = parser.parse_args()
    metadata = args.study_base / "metadata"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    shape_rows = read_csv(metadata / "shape_yields.csv")
    shape = indexed(shape_rows, "variant", "channel", "process")
    yearly_rows = read_csv(metadata / "yearly_sr_yields.csv")
    audit_rows = read_csv(metadata / "tau_audit_focus.csv")
    audit_summary_rows = read_csv(metadata / "tau_audit_summary.csv")
    histograms = read_wz_histograms(args.study_base)
    weighted_candidates = read_weighted_wz_candidate_composition(
        args.study_base
    )
    weighted_candidates_sr = (
        read_weighted_wz_candidates_in_sr(args.study_base)
    )

    set_style()
    summary = {
        "inputs": {
            "study_base": str(args.study_base),
            "shape_yields": str(metadata / "shape_yields.csv"),
            "yearly_sr_yields": str(metadata / "yearly_sr_yields.csv"),
            "tau_audit_focus": str(metadata / "tau_audit_focus.csv"),
            "tau_audit_summary": str(metadata / "tau_audit_summary.csv"),
        },
        "notes": [
            "All yields and shapes are pre-fit MC.",
            "The gen-match composition is candidate-level before the tau veto, not a final-SR event classification.",
            "The SR composition uses tau_none Vars event IDs passing lepton_cat != 2 and ptmiss >= 100, joined to each active working point's TauVetoAudit.",
            "The SR veto-event plot counts each event once using the exclusive priority genuine tau_h, electron-related, muon-related, then jet/unmatched.",
            "Nested retention uncertainties include the subset covariance approximation.",
        ],
        "wz_retention": {},
    }
    plot_wz_retention(shape, args.output_dir, summary)
    plot_wz_ptmiss(histograms, args.output_dir, summary)
    plot_transfer_factor(shape, args.output_dir, summary)
    plot_signal_wz_tradeoff(shape, args.output_dir, summary)
    plot_wz_prefit_postfit(shape, args.study_base, args.output_dir, summary)
    plot_yearly_wz_retention(yearly_rows, args.output_dir, summary)
    plot_wz_candidate_composition(audit_rows, args.output_dir, summary)
    plot_weighted_wz_candidate_composition(
        weighted_candidates, args.output_dir, summary
    )
    plot_weighted_wz_veto_events_in_sr(
        weighted_candidates_sr,
        args.output_dir,
        summary,
    )
    plot_wz_subsample_rejection(audit_summary_rows, args.output_dir, summary)

    with (args.output_dir / "wz_impact_plot_summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print("Wrote plots and summary to {}".format(args.output_dir))


if __name__ == "__main__":
    main()
