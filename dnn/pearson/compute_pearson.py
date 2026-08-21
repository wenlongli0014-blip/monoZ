#!/usr/bin/env python3

"""Compute Tiwari-style Pearson correlation diagnostics for the HZZ DNN.

The calculation deliberately mirrors the current training inputs:
  * the explicit candidate feature registry from each DNN config;
  * the SR selection implemented by train_dnn_score.load_training_arrays;
  * the 50% class-stratified training split with seed 12345;
  * non-negative physical weights clipped at the configured weight_clip.

Weighted correlations are produced separately for all events, signal, and
background. Pairwise missing values use the feature sentinel and are excluded
independently for each feature pair.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path

# Keep matplotlib caches out of the read-only home configuration directory.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/liwe/matplotlib-pearson")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT


THIS_DIR = Path(__file__).resolve().parent
DNN_DIR = THIS_DIR.parent

BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]


DEFAULT_CONFIGS = {
    "2017": DNN_DIR / "2017" / "dnn_config_2017.json",
    "2018": DNN_DIR / "2018" / "dnn_config_2018.json",
    "combined_2017_2018": DNN_DIR / "combined_2017_2018" / "dnn_config_combined_2017_2018.json",
}


@dataclass
class PearsonResult:
    correlation: np.ndarray
    pair_weight_sum: np.ndarray
    effective_events: np.ndarray
    pair_event_count: np.ndarray


def _index_alias(branch: str, index: int) -> str:
    return f"__pearson_{branch}_{index}"


def _required_base_branches(feature_specs: list[dict]) -> set[str]:
    branches = {"weight", "lepton_cat", "ptmiss"}
    for spec in feature_specs:
        if spec.get("kind") == "derived":
            branches.update(spec.get("branches", []))
        elif spec.get("branch"):
            branches.add(str(spec["branch"]))
    return branches


def _required_indices(feature_specs: list[dict]) -> set[tuple[str, int]]:
    indices: set[tuple[str, int]] = set()
    for spec in feature_specs:
        kind = spec.get("kind")
        if kind in {"fixed_index", "jagged_index"}:
            indices.add((str(spec["branch"]), int(spec["index"])))
        if kind == "derived":
            for branch in spec.get("branches", []):
                if branch.startswith("lepton_"):
                    indices.update({(branch, 0), (branch, 1)})
                elif branch.startswith("jet_") and branch != "jet_size":
                    indices.update({(branch, 0), (branch, 1)})
    return indices


def _sanitize(values: np.ndarray, sentinel: float) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    out[~np.isfinite(out)] = sentinel
    return out


def _valid(*values: np.ndarray, sentinel: float) -> np.ndarray:
    mask = np.ones(len(values[0]), dtype=bool)
    for value in values:
        mask &= np.isfinite(value) & (value != sentinel)
    return mask


def _wrapped_abs_dphi(phi1: np.ndarray, phi2: np.ndarray, sentinel: float) -> np.ndarray:
    out = np.full(len(phi1), sentinel, dtype=np.float64)
    mask = _valid(phi1, phi2, sentinel=sentinel)
    difference = phi1[mask] - phi2[mask]
    out[mask] = np.abs(np.arctan2(np.sin(difference), np.cos(difference)))
    return out


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray, sentinel: float) -> np.ndarray:
    out = np.full(len(numerator), sentinel, dtype=np.float64)
    mask = _valid(numerator, denominator, sentinel=sentinel) & (np.abs(denominator) > 1e-12)
    out[mask] = numerator[mask] / denominator[mask]
    return out


def _extract_feature_matrix(
    arrays: dict[str, np.ndarray],
    feature_specs: list[dict],
    sentinel: float,
) -> np.ndarray:
    """Reproduce the current dnn_common feature formulas from numeric arrays."""
    length = len(arrays["weight"])
    cache: dict[tuple, np.ndarray] = {}

    def scalar(branch: str) -> np.ndarray:
        key = ("scalar", branch)
        if key not in cache:
            cache[key] = _sanitize(arrays.get(branch, np.full(length, sentinel)), sentinel)
        return cache[key]

    def indexed(branch: str, index: int) -> np.ndarray:
        key = ("indexed", branch, index)
        if key not in cache:
            cache[key] = _sanitize(arrays.get(_index_alias(branch, index), np.full(length, sentinel)), sentinel)
        return cache[key]

    def derived(name: str) -> np.ndarray:
        key = ("derived", name)
        if key in cache:
            return cache[key]

        ll_pt = scalar("ll_pt")
        met = scalar("ptmiss")
        if name == "dphi_ll_ptmiss":
            out = _wrapped_abs_dphi(scalar("ll_phi"), scalar("ptmiss_phi"), sentinel)
        elif name == "mt_ll_ptmiss":
            mass_ll = scalar("ll_mass")
            dphi = _wrapped_abs_dphi(scalar("ll_phi"), scalar("ptmiss_phi"), sentinel)
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(ll_pt, met, mass_ll, dphi, sentinel=sentinel) & (ll_pt >= 0.0) & (met >= 0.0)
            invisible_z_mass = 91.1876
            et_ll = np.sqrt(np.maximum(ll_pt[mask] ** 2 + mass_ll[mask] ** 2, 0.0))
            et_met = np.sqrt(np.maximum(met[mask] ** 2 + invisible_z_mass**2, 0.0))
            mt2 = mass_ll[mask] ** 2 + invisible_z_mass**2 + 2.0 * (
                et_ll * et_met - ll_pt[mask] * met[mask] * np.cos(dphi[mask])
            )
            out[mask] = np.sqrt(np.maximum(mt2, 0.0))
        elif name == "ptmiss_over_ll_pt":
            out = _safe_ratio(met, ll_pt, sentinel)
        elif name == "pt_balance_ll_ptmiss":
            out = _safe_ratio(np.abs(met - ll_pt), met + ll_pt, sentinel)
        elif name == "deta_leptons":
            eta0, eta1 = indexed("lepton_eta", 0), indexed("lepton_eta", 1)
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(eta0, eta1, sentinel=sentinel)
            out[mask] = np.abs(eta0[mask] - eta1[mask])
        elif name == "dphi_leptons":
            out = _wrapped_abs_dphi(indexed("lepton_phi", 0), indexed("lepton_phi", 1), sentinel)
        elif name == "dr_leptons":
            deta, dphi = derived("deta_leptons"), derived("dphi_leptons")
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(deta, dphi, sentinel=sentinel)
            out[mask] = np.hypot(deta[mask], dphi[mask])
        elif name == "lepton_pt_ratio":
            out = _safe_ratio(indexed("lepton_pt", 1), indexed("lepton_pt", 0), sentinel)
        elif name == "dphi_lepton0_ptmiss":
            out = _wrapped_abs_dphi(indexed("lepton_phi", 0), scalar("ptmiss_phi"), sentinel)
        elif name == "dphi_lepton1_ptmiss":
            out = _wrapped_abs_dphi(indexed("lepton_phi", 1), scalar("ptmiss_phi"), sentinel)
        elif name == "lepton_ht":
            p0, p1 = indexed("lepton_pt", 0), indexed("lepton_pt", 1)
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(p0, p1, sentinel=sentinel)
            out[mask] = p0[mask] + p1[mask]
        elif name == "deta_jets":
            eta0, eta1 = indexed("jet_eta", 0), indexed("jet_eta", 1)
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(eta0, eta1, sentinel=sentinel)
            out[mask] = np.abs(eta0[mask] - eta1[mask])
        elif name == "dphi_jets":
            out = _wrapped_abs_dphi(indexed("jet_phi", 0), indexed("jet_phi", 1), sentinel)
        elif name == "dr_jets":
            deta, dphi = derived("deta_jets"), derived("dphi_jets")
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(deta, dphi, sentinel=sentinel)
            out[mask] = np.hypot(deta[mask], dphi[mask])
        elif name in {"dijet_mass", "dijet_pt"}:
            p0, p1 = indexed("jet_pt", 0), indexed("jet_pt", 1)
            phi0, phi1 = indexed("jet_phi", 0), indexed("jet_phi", 1)
            out = np.full(length, sentinel, dtype=np.float64)
            if name == "dijet_pt":
                mask = _valid(p0, p1, phi0, phi1, sentinel=sentinel)
                px = p0[mask] * np.cos(phi0[mask]) + p1[mask] * np.cos(phi1[mask])
                py = p0[mask] * np.sin(phi0[mask]) + p1[mask] * np.sin(phi1[mask])
                out[mask] = np.hypot(px, py)
            else:
                eta0, eta1 = indexed("jet_eta", 0), indexed("jet_eta", 1)
                mass0, mass1 = indexed("jet_mass", 0), indexed("jet_mass", 1)
                mask = _valid(p0, p1, eta0, eta1, phi0, phi1, mass0, mass1, sentinel=sentinel)
                e0 = np.sqrt(np.maximum((p0[mask] * np.cosh(eta0[mask])) ** 2 + mass0[mask] ** 2, 0.0))
                e1 = np.sqrt(np.maximum((p1[mask] * np.cosh(eta1[mask])) ** 2 + mass1[mask] ** 2, 0.0))
                px = p0[mask] * np.cos(phi0[mask]) + p1[mask] * np.cos(phi1[mask])
                py = p0[mask] * np.sin(phi0[mask]) + p1[mask] * np.sin(phi1[mask])
                pz = p0[mask] * np.sinh(eta0[mask]) + p1[mask] * np.sinh(eta1[mask])
                out[mask] = np.sqrt(np.maximum((e0 + e1) ** 2 - px**2 - py**2 - pz**2, 0.0))
        elif name == "jet_ht":
            p0, p1 = indexed("jet_pt", 0), indexed("jet_pt", 1)
            out = np.zeros(length, dtype=np.float64)
            valid0, valid1 = _valid(p0, sentinel=sentinel), _valid(p1, sentinel=sentinel)
            out[valid0] += p0[valid0]
            out[valid1] += p1[valid1]
            out[~(valid0 | valid1)] = sentinel
        elif name == "min_dphi_jet_ptmiss":
            dphi0 = _wrapped_abs_dphi(indexed("jet_phi", 0), scalar("ptmiss_phi"), sentinel)
            dphi1 = _wrapped_abs_dphi(indexed("jet_phi", 1), scalar("ptmiss_phi"), sentinel)
            out = np.full(length, sentinel, dtype=np.float64)
            valid0, valid1 = _valid(dphi0, sentinel=sentinel), _valid(dphi1, sentinel=sentinel)
            out[valid0 & ~valid1] = dphi0[valid0 & ~valid1]
            out[valid1 & ~valid0] = dphi1[valid1 & ~valid0]
            both = valid0 & valid1
            out[both] = np.minimum(dphi0[both], dphi1[both])
        else:
            raise ValueError(f"Unsupported derived feature: {name}")

        cache[key] = out
        return out

    columns = []
    for spec in feature_specs:
        kind = spec["kind"]
        if kind == "scalar":
            column = scalar(str(spec["branch"]))
        elif kind in {"fixed_index", "jagged_index"}:
            column = indexed(str(spec["branch"]), int(spec["index"]))
        elif kind == "derived":
            column = derived(str(spec["name"]))
        else:
            raise ValueError(f"Unknown feature kind: {kind}")
        columns.append(_sanitize(column, sentinel).astype(np.float32))
    return np.column_stack(columns).astype(np.float32)


def _read_process(
    path: Path,
    label: int,
    process: str,
    feature_specs: list[dict],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read one process with the same branches and SR selection as training."""
    sentinel = float(feature_specs[0].get("missing_sentinel", -9999.0))
    dataframe = ROOT.RDataFrame("Vars", str(path)).Filter(
        "lepton_cat != 2 && ptmiss >= 100 && std::isfinite(ptmiss)"
    )
    available = {str(name) for name in dataframe.GetColumnNames()}
    missing = sorted(_required_base_branches(feature_specs) - available)
    if missing:
        raise KeyError(f"Missing branches in {path}: {missing}")

    aliases = []
    for branch, index in sorted(_required_indices(feature_specs)):
        alias = _index_alias(branch, index)
        expression = f"{branch}.size() > {index} ? static_cast<double>({branch}[{index}]) : {sentinel}"
        dataframe = dataframe.Define(alias, expression)
        aliases.append(alias)

    scalar_branches = {
        str(spec["branch"])
        for spec in feature_specs
        if spec.get("kind") == "scalar"
    }
    for spec in feature_specs:
        if spec.get("kind") == "derived":
            scalar_branches.update(
                branch
                for branch in spec.get("branches", [])
                if not branch.startswith("lepton_") and not branch.startswith("jet_")
            )
    requested = sorted(scalar_branches | {"weight"} | set(aliases))
    arrays = {name: np.asarray(values) for name, values in dataframe.AsNumpy(requested).items()}

    raw_features = _extract_feature_matrix(arrays, feature_specs, sentinel)
    signed_weights = np.asarray(arrays["weight"], dtype=np.float32)
    signed_weights = np.where(np.isfinite(signed_weights), signed_weights, 0.0).astype(np.float32)
    n_selected = len(signed_weights)
    return (
        raw_features.astype(np.float32),
        np.full(n_selected, label, dtype=np.int64),
        signed_weights,
        np.full(n_selected, process),
    )


def _load_training_arrays(
    input_bases: list[str],
    signal_file: str,
    feature_specs: list[dict],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pieces = []
    for source_index, input_base in enumerate(input_bases):
        source = Path(input_base).name or f"source{source_index}"
        inputs = [(signal_file, 1, f"{source}:signal")]
        inputs.extend((f"{process}.root", 0, f"{source}:{process}") for process in BKG_PROCESSES)
        for filename, label, process_name in inputs:
            piece = _read_process(Path(input_base) / "SR" / filename, label, process_name, feature_specs)
            print(f"Loaded {process_name}: {len(piece[1])} events", flush=True)
            pieces.append(piece)
    return tuple(np.concatenate([piece[index] for piece in pieces]) for index in range(4))


def _nonnegative_clipped_weights(signed_weights: np.ndarray, clip: float) -> np.ndarray:
    weights = np.where(np.isfinite(signed_weights), signed_weights, 0.0).astype(np.float64)
    weights = np.maximum(weights, 0.0)
    if clip > 0.0:
        weights = np.minimum(weights, clip)
    return weights


def _stratified_train_indices(labels: np.ndarray, fraction: float, seed: int) -> np.ndarray:
    """Pure-NumPy reproducible class-stratified training subset."""
    rng = np.random.RandomState(seed)
    selected = []
    for label in np.unique(labels):
        indices = np.flatnonzero(labels == label)
        rng.shuffle(indices)
        selected.append(indices[: int(np.floor(fraction * len(indices)))])
    train_indices = np.concatenate(selected)
    rng.shuffle(train_indices)
    return train_indices


def _accumulate_pairwise_moments(
    x: np.ndarray,
    weights: np.ndarray,
    *,
    sentinel: float,
    chunk_size: int,
) -> PearsonResult:
    """Return a pairwise-missing weighted Pearson matrix.

    Sufficient statistics are accumulated with small matrix products so the
    full event-by-feature array is never duplicated as float64.
    """
    x = np.asarray(x)
    weights = np.asarray(weights, dtype=np.float64)
    n_features = x.shape[1]

    pair_weight_sum = np.zeros((n_features, n_features), dtype=np.float64)
    pair_weight2_sum = np.zeros_like(pair_weight_sum)
    pair_event_count = np.zeros_like(pair_weight_sum)
    sum_x = np.zeros_like(pair_weight_sum)
    sum_x2 = np.zeros_like(pair_weight_sum)
    sum_xy = np.zeros_like(pair_weight_sum)

    for start in range(0, len(x), chunk_size):
        stop = min(start + chunk_size, len(x))
        values = np.asarray(x[start:stop], dtype=np.float64)
        w = weights[start:stop]

        valid = np.isfinite(values) & (values != float(sentinel))
        valid &= np.isfinite(w)[:, None] & (w > 0.0)[:, None]
        valid_f = valid.astype(np.float64)
        clean = np.where(valid, values, 0.0)
        w_col = np.where(np.isfinite(w) & (w > 0.0), w, 0.0)[:, None]

        weighted_valid = valid_f * w_col
        weighted_clean = clean * w_col

        pair_event_count += valid_f.T @ valid_f
        pair_weight_sum += weighted_valid.T @ valid_f
        pair_weight2_sum += (valid_f * (w_col * w_col)).T @ valid_f
        sum_x += weighted_clean.T @ valid_f
        sum_x2 += ((clean * clean) * w_col).T @ valid_f
        sum_xy += weighted_clean.T @ clean

    with np.errstate(divide="ignore", invalid="ignore"):
        mean_x = sum_x / pair_weight_sum
        mean_y = mean_x.T
        covariance = sum_xy / pair_weight_sum - mean_x * mean_y
        variance_x = sum_x2 / pair_weight_sum - mean_x * mean_x
        variance_y = variance_x.T
        denominator = np.sqrt(np.maximum(variance_x, 0.0) * np.maximum(variance_y, 0.0))
        correlation = covariance / denominator
        effective_events = pair_weight_sum * pair_weight_sum / pair_weight2_sum

    bad = (pair_weight_sum <= 0.0) | (denominator <= 0.0) | ~np.isfinite(correlation)
    correlation[bad] = np.nan
    correlation = np.clip(correlation, -1.0, 1.0)
    effective_events[~np.isfinite(effective_events)] = 0.0

    return PearsonResult(
        correlation=correlation,
        pair_weight_sum=pair_weight_sum,
        effective_events=effective_events,
        pair_event_count=pair_event_count,
    )


def _write_matrix_csv(path: Path, names: list[str], matrix: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature"] + names)
        for name, row in zip(names, matrix):
            writer.writerow([name] + ["" if not np.isfinite(v) else f"{v:.10g}" for v in row])


def _plot_matrix(path: Path, names: list[str], matrix: np.ndarray, title: str) -> None:
    n_features = len(names)
    size = max(9.0, 0.47 * n_features)
    fig, ax = plt.subplots(figsize=(size + 2.0, size))
    masked = np.ma.masked_invalid(matrix)
    image = ax.imshow(masked, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="equal")

    ticks = np.arange(n_features)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticklabels(names, fontsize=7)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=13, pad=14)

    cell_font = 4.0 if n_features > 25 else 7.0
    for row in range(n_features):
        for col in range(n_features):
            value = matrix[row, col]
            if not np.isfinite(value):
                continue
            color = "white" if abs(value) >= 0.58 else "black"
            ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=cell_font, color=color)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    colorbar.set_label("Pearson correlation coefficient")
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def _write_pair_summary(
    path: Path,
    names: list[str],
    results: dict[tuple[str, str], PearsonResult],
) -> None:
    columns = [f"{scope}_weighted" for scope in ("all", "signal", "background")]

    rows = []
    for i, name_i in enumerate(names):
        for j in range(i + 1, len(names)):
            values = {column: results[tuple(column.rsplit("_", 1))].correlation[i, j] for column in columns}
            finite = [(column, value) for column, value in values.items() if np.isfinite(value)]
            max_scope, max_value = max(finite, key=lambda item: abs(item[1])) if finite else ("", np.nan)
            rows.append({
                "feature_1": name_i,
                "feature_2": names[j],
                **values,
                "max_abs_correlation": abs(max_value) if np.isfinite(max_value) else np.nan,
                "max_correlation_scope": max_scope,
                "max_signed_correlation": max_value,
                "signal_weighted_effective_events": results[("signal", "weighted")].effective_events[i, j],
                "background_weighted_effective_events": results[("background", "weighted")].effective_events[i, j],
            })

    rows.sort(key=lambda row: (-np.nan_to_num(row["max_abs_correlation"], nan=-1.0), row["feature_1"], row["feature_2"]))
    fieldnames = list(rows[0]) if rows else ["feature_1", "feature_2"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_strong_pairs(path: Path, pair_summary_path: Path, threshold: float) -> int:
    with pair_summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    strong = [row for row in rows if float(row["max_abs_correlation"]) >= threshold]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["feature_1", "feature_2"])
        writer.writeheader()
        writer.writerows(strong)
    return len(strong)


def _run_feature_set(
    outdir: Path,
    dataset_name: str,
    set_name: str,
    x_train: np.ndarray,
    y_train: np.ndarray,
    physical_weights: np.ndarray,
    names: list[str],
    sentinel: float,
    chunk_size: int,
    strong_threshold: float,
) -> dict:
    set_dir = outdir / dataset_name / set_name
    set_dir.mkdir(parents=True, exist_ok=True)
    results: dict[tuple[str, str], PearsonResult] = {}

    scopes = {
        "all": np.ones(len(y_train), dtype=bool),
        "signal": y_train == 1,
        "background": y_train == 0,
    }
    for scope, mask in scopes.items():
        weighting = "weighted"
        print(f"[{dataset_name}/{set_name}] {scope} {weighting}: {np.count_nonzero(mask)} events", flush=True)
        result = _accumulate_pairwise_moments(
            x_train[mask], physical_weights[mask], sentinel=sentinel, chunk_size=chunk_size,
        )
        results[(scope, weighting)] = result
        stem = f"pearson_{scope}_{weighting}"
        _write_matrix_csv(set_dir / f"{stem}.csv", names, result.correlation)
        _write_matrix_csv(set_dir / f"{stem}_effective_events.csv", names, result.effective_events)
        _write_matrix_csv(set_dir / f"{stem}_pair_event_count.csv", names, result.pair_event_count)
        _plot_matrix(
            set_dir / f"{stem}.png",
            names,
            result.correlation,
            f"{dataset_name} {set_name}: {scope}, weighted Pearson (train only)",
        )

    np.savez_compressed(
        set_dir / "pearson_matrices.npz",
        feature_names=np.asarray(names),
        **{
            f"{scope}_{weighting}_{quantity}": getattr(result, quantity)
            for (scope, weighting), result in results.items()
            for quantity in ("correlation", "pair_weight_sum", "effective_events", "pair_event_count")
        },
    )
    pair_summary = set_dir / "pearson_pairs.csv"
    _write_pair_summary(pair_summary, names, results)
    n_strong = _write_strong_pairs(set_dir / "strong_pairs.csv", pair_summary, strong_threshold)
    return {"n_features": len(names), "n_strong_pairs": n_strong, "directory": str(set_dir)}


def run_dataset(
    dataset_name: str,
    config_path: Path,
    outdir: Path,
    *,
    seed: int,
    chunk_size: int,
    strong_threshold: float,
) -> dict:
    config = json.loads(config_path.read_text())
    configured_inputs = config.get("input_bases", config.get("input_base"))
    input_bases = configured_inputs if isinstance(configured_inputs, list) else [configured_inputs]
    signal_file = config["signal_file"]
    sentinel = float(config.get("missing_sentinel", -9999.0))
    weight_clip = float(config.get("weight_clip", 100.0))

    metadata_path = DNN_DIR / dataset_name / "output" / "dnn_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    feature_specs = list(metadata["candidate_feature_specs"])
    candidate_names = [str(spec["name"]) for spec in feature_specs]
    skipped = list(metadata.get("skipped_branches", []))

    x_all, y_all, signed_weights, process = _load_training_arrays(
        input_bases, signal_file, feature_specs,
    )
    physical_weights = _nonnegative_clipped_weights(signed_weights, weight_clip)

    train_idx = _stratified_train_indices(y_all, fraction=0.50, seed=seed)
    x_train = x_all[train_idx]
    y_train = y_all[train_idx]
    w_train = physical_weights[train_idx]

    selected_names = list(metadata["feature_names"])
    selected_indices = [candidate_names.index(name) for name in selected_names]

    feature_sets = {
        "candidate_features": (x_train, candidate_names),
        "selected_top10": (x_train[:, selected_indices], selected_names),
    }
    summaries = {}
    for set_name, (features, names) in feature_sets.items():
        summaries[set_name] = _run_feature_set(
            outdir,
            dataset_name,
            set_name,
            features,
            y_train,
            w_train,
            names,
            sentinel,
            chunk_size,
            strong_threshold,
        )

    dataset_summary = {
        "dataset": dataset_name,
        "config": str(config_path),
        "input_bases": input_bases,
        "seed": seed,
        "selection": "lepton_cat != 2 and ptmiss >= 100 (via current training loader)",
        "correlation_scope": "50% class-stratified training split only",
        "n_all_candidates": int(len(y_all)),
        "n_train": int(len(train_idx)),
        "n_train_signal": int(np.count_nonzero(y_train == 1)),
        "n_train_background": int(np.count_nonzero(y_train == 0)),
        "negative_weight_events_in_train": int(np.count_nonzero(signed_weights[train_idx] < 0.0)),
        "zero_physical_weight_events_in_train": int(np.count_nonzero(w_train <= 0.0)),
        "weight_clip": weight_clip,
        "strong_pair_threshold": strong_threshold,
        "candidate_feature_names": candidate_names,
        "selected_feature_names": selected_names,
        "skipped_branches": skipped,
        "processes_in_train": {str(name): int(np.count_nonzero(process[train_idx] == name)) for name in np.unique(process[train_idx])},
        "feature_sets": summaries,
    }
    dataset_dir = outdir / dataset_name
    (dataset_dir / "summary.json").write_text(json.dumps(dataset_summary, indent=2) + "\n")

    with (dataset_dir / "feature_provenance.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["feature", "category", "kind", "root_dependencies", "selected_top10"],
        )
        writer.writeheader()
        for spec in feature_specs:
            kind = str(spec["kind"])
            if kind == "scalar":
                category = "root_scalar"
                dependencies = [str(spec["branch"])]
            elif kind in {"fixed_index", "jagged_index"}:
                category = "root_vector_element"
                dependencies = [f"{spec['branch']}[{spec['index']}]"]
            else:
                category = "derived"
                dependencies = [str(branch) for branch in spec.get("branches", [])]
            writer.writerow({
                "feature": spec["name"],
                "category": category,
                "kind": kind,
                "root_dependencies": ";".join(dependencies),
                "selected_top10": str(spec["name"] in selected_names).lower(),
            })
    return dataset_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=list(DEFAULT_CONFIGS),
        default=list(DEFAULT_CONFIGS),
        help="Datasets to process (default: all).",
    )
    parser.add_argument("--outdir", type=Path, default=THIS_DIR)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--strong-threshold", type=float, default=0.80)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for dataset in args.datasets:
        summaries.append(
            run_dataset(
                dataset,
                DEFAULT_CONFIGS[dataset],
                args.outdir,
                seed=args.seed,
                chunk_size=args.chunk_size,
                strong_threshold=args.strong_threshold,
            )
        )

    run_summary = {
        "method": "Tiwari-style pairwise-missing Pearson correlation",
        "weighted_definition": "non-negative physical event weights clipped identically to current DNN training",
        "datasets": summaries,
    }
    (args.outdir / "run_summary.json").write_text(json.dumps(run_summary, indent=2) + "\n")
    print(f"Pearson outputs written to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
