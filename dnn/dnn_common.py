#!/usr/bin/env python3

from __future__ import annotations

import copy
import json
from pathlib import Path

import awkward as ak
import numpy as np
import torch
from torch import nn


BKG_PROCESSES = ["DY", "Other", "ST", "VVV", "WW", "WZ", "ggZZ", "qqZZ", "ttbar"]

CHANNELS = [
    ("SR", "SR", "lepton_cat!=2 && ptmiss>=100", True),
    ("DYCR", "DY_CR", "lepton_cat!=2 && ptmiss>=30 && ptmiss<=90", False),
    ("EMUCR", "emu_CR", "lepton_cat==2 && ptmiss>=100", False),
    ("CR3L", "3l_CR", "ptmiss>=30", False),
]

DEFAULT_CONFIG_PATH = Path("/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/dnn_config.json")
MISSING_SENTINEL = -9999.0


# The model inputs are intentionally explicit. Absolute phi values are only
# dependencies of wrapped angular differences; they are not model features.
DEFAULT_FEATURE_NAMES = [
    "lepton_cat",
    "jet_cat",
    "ll_pt",
    "ll_eta",
    "ll_mass",
    "ptmiss",
    "ptmiss_significance_corrected",
    "lepton_pt_0",
    "lepton_pt_1",
    "lepton_eta_0",
    "lepton_eta_1",
    "jet_size",
    "jet_pt_0",
    "jet_pt_1",
    "jet_eta_0",
    "jet_eta_1",
    "jet_mass_0",
    "jet_mass_1",
    "dphi_ll_ptmiss",
    "mt_ll_ptmiss",
    "ptmiss_over_ll_pt",
    "pt_balance_ll_ptmiss",
    "deta_leptons",
    "dphi_leptons",
    "dr_leptons",
    "lepton_pt_ratio",
    "dphi_lepton0_ptmiss",
    "dphi_lepton1_ptmiss",
    "lepton_ht",
    "deta_jets",
    "dphi_jets",
    "dr_jets",
    "dijet_mass",
    "dijet_pt",
    "jet_ht",
    "min_dphi_jet_ptmiss",
]


RAW_FEATURE_SPECS = {
    "lepton_cat": {"kind": "scalar", "branch": "lepton_cat"},
    "jet_cat": {"kind": "scalar", "branch": "jet_cat"},
    "ll_pt": {"kind": "scalar", "branch": "ll_pt"},
    "ll_eta": {"kind": "scalar", "branch": "ll_eta"},
    "ll_mass": {"kind": "scalar", "branch": "ll_mass"},
    "ptmiss": {"kind": "scalar", "branch": "ptmiss"},
    "ptmiss_significance_corrected": {"kind": "scalar", "branch": "ptmiss_significance_corrected"},
    "lepton_pt_0": {"kind": "fixed_index", "branch": "lepton_pt", "index": 0},
    "lepton_pt_1": {"kind": "fixed_index", "branch": "lepton_pt", "index": 1},
    "lepton_eta_0": {"kind": "fixed_index", "branch": "lepton_eta", "index": 0},
    "lepton_eta_1": {"kind": "fixed_index", "branch": "lepton_eta", "index": 1},
    "jet_size": {"kind": "scalar", "branch": "jet_size"},
    "jet_pt_0": {"kind": "jagged_index", "branch": "jet_pt", "index": 0},
    "jet_pt_1": {"kind": "jagged_index", "branch": "jet_pt", "index": 1},
    "jet_eta_0": {"kind": "jagged_index", "branch": "jet_eta", "index": 0},
    "jet_eta_1": {"kind": "jagged_index", "branch": "jet_eta", "index": 1},
    "jet_mass_0": {"kind": "jagged_index", "branch": "jet_mass", "index": 0},
    "jet_mass_1": {"kind": "jagged_index", "branch": "jet_mass", "index": 1},
}


DERIVED_FEATURE_BRANCHES = {
    "dphi_ll_ptmiss": ["ll_phi", "ptmiss_phi"],
    "mt_ll_ptmiss": ["ll_pt", "ll_mass", "ll_phi", "ptmiss", "ptmiss_phi"],
    "ptmiss_over_ll_pt": ["ptmiss", "ll_pt"],
    "pt_balance_ll_ptmiss": ["ptmiss", "ll_pt"],
    "deta_leptons": ["lepton_eta"],
    "dphi_leptons": ["lepton_phi"],
    "dr_leptons": ["lepton_eta", "lepton_phi"],
    "lepton_pt_ratio": ["lepton_pt"],
    "dphi_lepton0_ptmiss": ["lepton_phi", "ptmiss_phi"],
    "dphi_lepton1_ptmiss": ["lepton_phi", "ptmiss_phi"],
    "lepton_ht": ["lepton_pt"],
    "deta_jets": ["jet_eta"],
    "dphi_jets": ["jet_phi"],
    "dr_jets": ["jet_eta", "jet_phi"],
    "dijet_mass": ["jet_pt", "jet_eta", "jet_phi", "jet_mass"],
    "dijet_pt": ["jet_pt", "jet_phi"],
    "jet_ht": ["jet_pt"],
    "min_dphi_jet_ptmiss": ["jet_phi", "ptmiss_phi"],
}


DEFAULT_CONFIG = {
    "input_base": "/eos/user/l/liwe/monoZ_combine",
    "scored_output_base": "/eos/user/l/liwe/monoZ_combine_dnn",
    "signal_file": "signal.root",
    "tree_name": "Vars",
    "training_region": "SR",
    "training_selection": "lepton_cat != 2 && ptmiss >= 100",
    "feature_names": DEFAULT_FEATURE_NAMES,
    "missing_sentinel": MISSING_SENTINEL,
    "preprocessing": "standardize_nonmissing_then_map_missing_to_zero",
    "top_k_significance": 10,
    "feature_significance_bins": 40,
    "feature_significance_background_syst": 0.20,
    "feature_selection_scope": "full_dataset_before_split",
    "feature_selection_known_issue": "test-set leakage by explicit user choice; revisit later",
    "weight_clip": 100.0,
    "split": {"train": 0.50, "validation": 0.20, "test": 0.30},
    "binning": {
        "SR": [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.875, 1.0],
        "DYCR": [0.0, 0.15, 0.30, 1.0],
        "EMUCR": [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 1.0],
        "CR3L": [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 1.0],
    },
    "rate_parameters": {
        "k_Zjet": ["DY"],
        "k_WZ": ["WZ"],
        "k_emu": ["ttbar", "WW", "ST"],
    },
    "notes": [
        "Data.root is scored but never used for supervised training.",
        "Candidate inputs are an explicit HZZ2l2nu physics-feature allowlist.",
        "Absolute phi values are used only to derive wrapped delta-phi features.",
        "Missing and nonfinite feature values use -9999 before scaling and map to zero after scaling.",
        "Training/ranking weights clamp negative values to zero and clip at 100; signed weights remain available for templates.",
    ],
}


class SimpleDNN(nn.Module):
    def __init__(
        self,
        input_size: int,
        hidden_size: int = 64,
        dropout: float = 0.08,
        hidden_layers: list[int] | tuple[int, ...] | None = None,
    ):
        super().__init__()
        if hidden_layers is None:
            # Legacy architecture used by the archived two-logit models.
            hidden_layers = [hidden_size, hidden_size, max(hidden_size // 2, 8)]
        layers: list[nn.Module] = []
        previous = int(input_size)
        for width in hidden_layers:
            layers.extend([nn.Linear(previous, int(width)), nn.ReLU(), nn.Dropout(float(dropout))])
            previous = int(width)
        layers.append(nn.Linear(previous, 2))
        self.net = nn.Sequential(*layers)
        self.hidden_layers = tuple(int(x) for x in hidden_layers)
        self.score_activation = "softmax_signal"

    def forward(self, x):
        return self.net(x)


def load_config(path: str | Path = DEFAULT_CONFIG_PATH):
    path = Path(path)
    if not path.exists():
        path.write_text(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
    config = copy.deepcopy(DEFAULT_CONFIG)
    config.update(json.loads(path.read_text()))
    return config


def build_feature_specs(config) -> list[dict]:
    requested = list(config.get("feature_names", DEFAULT_FEATURE_NAMES))
    specs: list[dict] = []
    unknown: list[str] = []
    sentinel = float(config.get("missing_sentinel", MISSING_SENTINEL))
    for name in requested:
        if name in RAW_FEATURE_SPECS:
            spec = {"name": name, **copy.deepcopy(RAW_FEATURE_SPECS[name]), "missing_sentinel": sentinel}
        elif name in DERIVED_FEATURE_BRANCHES:
            branches = list(DERIVED_FEATURE_BRANCHES[name])
            spec = {
                "name": name,
                "kind": "derived",
                "branches": branches,
                "source": f"derived({','.join(branches)})",
                "missing_sentinel": sentinel,
            }
        else:
            unknown.append(str(name))
            continue
        specs.append(spec)
    if unknown:
        raise ValueError(f"Unknown explicit DNN feature names: {unknown}")
    if not specs:
        raise ValueError("The explicit DNN feature list is empty")
    return specs


def discover_feature_specs(root_path: str | Path, config):
    """Backward-compatible name for the new explicit allowlist builder."""
    specs = build_feature_specs(config)
    import uproot

    with uproot.open(root_path)[config.get("tree_name", "Vars")] as tree:
        available = set(tree.keys())
    missing = sorted(set(branches_for_specs(specs)) - available)
    return specs, [(name, "missing_in_reference; will use sentinel") for name in missing]


def feature_names(feature_specs):
    return [spec["name"] for spec in feature_specs]


def branches_for_specs(feature_specs, extra=()):
    names: set[str] = set(extra)
    for spec in feature_specs:
        if spec.get("kind") == "derived":
            names.update(spec.get("branches", []))
        elif spec.get("branch"):
            names.add(spec["branch"])
    return sorted(names)


def available_branches(tree, requested) -> list[str]:
    present = set(tree.keys())
    return [name for name in requested if name in present]


def _fields(arrays) -> set[str]:
    if hasattr(arrays, "fields"):
        return set(arrays.fields)
    if hasattr(arrays, "keys"):
        return set(arrays.keys())
    return set()


def _sanitize(values, sentinel: float) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    out[~np.isfinite(out)] = sentinel
    return out


def _scalar(arrays, branch: str, length: int, sentinel: float) -> np.ndarray:
    if branch not in _fields(arrays):
        return np.full(length, sentinel, dtype=np.float64)
    return _sanitize(ak.to_numpy(arrays[branch]), sentinel)


def _fixed(arrays, branch: str, index: int, length: int, sentinel: float) -> np.ndarray:
    if branch not in _fields(arrays):
        return np.full(length, sentinel, dtype=np.float64)
    try:
        return _sanitize(ak.to_numpy(arrays[branch][:, index]), sentinel)
    except (IndexError, ValueError):
        return np.full(length, sentinel, dtype=np.float64)


def _jagged(arrays, branch: str, index: int, length: int, sentinel: float) -> np.ndarray:
    if branch not in _fields(arrays):
        return np.full(length, sentinel, dtype=np.float64)
    padded = ak.pad_none(arrays[branch], index + 1, clip=True)
    values = ak.fill_none(padded[:, index], sentinel)
    return _sanitize(ak.to_numpy(values), sentinel)


def _valid(*values: np.ndarray, sentinel: float) -> np.ndarray:
    mask = np.ones(len(values[0]), dtype=bool)
    for value in values:
        mask &= np.isfinite(value) & (value != sentinel)
    return mask


def _wrapped_abs_dphi(phi1: np.ndarray, phi2: np.ndarray, sentinel: float) -> np.ndarray:
    out = np.full(len(phi1), sentinel, dtype=np.float64)
    mask = _valid(phi1, phi2, sentinel=sentinel)
    out[mask] = np.abs(np.arctan2(np.sin(phi1[mask] - phi2[mask]), np.cos(phi1[mask] - phi2[mask])))
    return out


def _safe_ratio(num: np.ndarray, den: np.ndarray, sentinel: float) -> np.ndarray:
    out = np.full(len(num), sentinel, dtype=np.float64)
    mask = _valid(num, den, sentinel=sentinel) & (np.abs(den) > 1e-12)
    out[mask] = num[mask] / den[mask]
    return out


def _derived_feature(name: str, arrays, length: int, sentinel: float, cache=None) -> np.ndarray:
    cache = {} if cache is None else cache

    def cached(key, factory):
        if key not in cache:
            cache[key] = factory()
        return cache[key]

    ll_pt = lambda: cached(("scalar", "ll_pt"), lambda: _scalar(arrays, "ll_pt", length, sentinel))
    ll_phi = lambda: cached(("scalar", "ll_phi"), lambda: _scalar(arrays, "ll_phi", length, sentinel))
    met = lambda: cached(("scalar", "ptmiss"), lambda: _scalar(arrays, "ptmiss", length, sentinel))
    met_phi = lambda: cached(("scalar", "ptmiss_phi"), lambda: _scalar(arrays, "ptmiss_phi", length, sentinel))
    lep_pt = lambda index: cached(("fixed", "lepton_pt", index), lambda: _fixed(arrays, "lepton_pt", index, length, sentinel))
    lep_eta = lambda index: cached(("fixed", "lepton_eta", index), lambda: _fixed(arrays, "lepton_eta", index, length, sentinel))
    lep_phi = lambda index: cached(("fixed", "lepton_phi", index), lambda: _fixed(arrays, "lepton_phi", index, length, sentinel))
    jet_pt = lambda index: cached(("jagged", "jet_pt", index), lambda: _jagged(arrays, "jet_pt", index, length, sentinel))
    jet_eta = lambda index: cached(("jagged", "jet_eta", index), lambda: _jagged(arrays, "jet_eta", index, length, sentinel))
    jet_phi = lambda index: cached(("jagged", "jet_phi", index), lambda: _jagged(arrays, "jet_phi", index, length, sentinel))
    jet_mass = lambda index: cached(("jagged", "jet_mass", index), lambda: _jagged(arrays, "jet_mass", index, length, sentinel))

    if name == "dphi_ll_ptmiss":
        return _wrapped_abs_dphi(ll_phi(), met_phi(), sentinel)
    if name == "mt_ll_ptmiss":
        p_ll, p_met = ll_pt(), met()
        mass_ll = _scalar(arrays, "ll_mass", length, sentinel)
        dphi = _wrapped_abs_dphi(ll_phi(), met_phi(), sentinel)
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(p_ll, p_met, mass_ll, dphi, sentinel=sentinel) & (p_ll >= 0.0) & (p_met >= 0.0)
        invisible_z_mass = 91.1876
        et_ll = np.sqrt(np.maximum(p_ll[mask] ** 2 + mass_ll[mask] ** 2, 0.0))
        et_met = np.sqrt(np.maximum(p_met[mask] ** 2 + invisible_z_mass**2, 0.0))
        mt2 = (
            mass_ll[mask] ** 2
            + invisible_z_mass**2
            + 2.0 * (et_ll * et_met - p_ll[mask] * p_met[mask] * np.cos(dphi[mask]))
        )
        out[mask] = np.sqrt(np.maximum(mt2, 0.0))
        return out
    if name == "ptmiss_over_ll_pt":
        return _safe_ratio(met(), ll_pt(), sentinel)
    if name == "pt_balance_ll_ptmiss":
        p_ll, p_met = ll_pt(), met()
        return _safe_ratio(np.abs(p_met - p_ll), p_met + p_ll, sentinel)
    if name == "deta_leptons":
        eta0, eta1 = lep_eta(0), lep_eta(1)
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(eta0, eta1, sentinel=sentinel)
        out[mask] = np.abs(eta0[mask] - eta1[mask])
        return out
    if name == "dphi_leptons":
        return _wrapped_abs_dphi(lep_phi(0), lep_phi(1), sentinel)
    if name == "dr_leptons":
        deta = cached(("derived", "deta_leptons"), lambda: _derived_feature("deta_leptons", arrays, length, sentinel, cache))
        dphi = cached(("derived", "dphi_leptons"), lambda: _derived_feature("dphi_leptons", arrays, length, sentinel, cache))
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(deta, dphi, sentinel=sentinel)
        out[mask] = np.hypot(deta[mask], dphi[mask])
        return out
    if name == "lepton_pt_ratio":
        return _safe_ratio(lep_pt(1), lep_pt(0), sentinel)
    if name == "dphi_lepton0_ptmiss":
        return _wrapped_abs_dphi(lep_phi(0), met_phi(), sentinel)
    if name == "dphi_lepton1_ptmiss":
        return _wrapped_abs_dphi(lep_phi(1), met_phi(), sentinel)
    if name == "lepton_ht":
        p0, p1 = lep_pt(0), lep_pt(1)
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(p0, p1, sentinel=sentinel)
        out[mask] = p0[mask] + p1[mask]
        return out
    if name == "deta_jets":
        eta0, eta1 = jet_eta(0), jet_eta(1)
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(eta0, eta1, sentinel=sentinel)
        out[mask] = np.abs(eta0[mask] - eta1[mask])
        return out
    if name == "dphi_jets":
        return _wrapped_abs_dphi(jet_phi(0), jet_phi(1), sentinel)
    if name == "dr_jets":
        deta = cached(("derived", "deta_jets"), lambda: _derived_feature("deta_jets", arrays, length, sentinel, cache))
        dphi = cached(("derived", "dphi_jets"), lambda: _derived_feature("dphi_jets", arrays, length, sentinel, cache))
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(deta, dphi, sentinel=sentinel)
        out[mask] = np.hypot(deta[mask], dphi[mask])
        return out
    if name in {"dijet_mass", "dijet_pt"}:
        p0, p1 = jet_pt(0), jet_pt(1)
        phi0, phi1 = jet_phi(0), jet_phi(1)
        if name == "dijet_pt":
            out = np.full(length, sentinel, dtype=np.float64)
            mask = _valid(p0, p1, phi0, phi1, sentinel=sentinel)
            px = p0[mask] * np.cos(phi0[mask]) + p1[mask] * np.cos(phi1[mask])
            py = p0[mask] * np.sin(phi0[mask]) + p1[mask] * np.sin(phi1[mask])
            out[mask] = np.hypot(px, py)
            return out
        eta0, eta1 = jet_eta(0), jet_eta(1)
        mass0, mass1 = jet_mass(0), jet_mass(1)
        out = np.full(length, sentinel, dtype=np.float64)
        mask = _valid(p0, p1, eta0, eta1, phi0, phi1, mass0, mass1, sentinel=sentinel)
        e0 = np.sqrt(np.maximum((p0[mask] * np.cosh(eta0[mask])) ** 2 + mass0[mask] ** 2, 0.0))
        e1 = np.sqrt(np.maximum((p1[mask] * np.cosh(eta1[mask])) ** 2 + mass1[mask] ** 2, 0.0))
        px = p0[mask] * np.cos(phi0[mask]) + p1[mask] * np.cos(phi1[mask])
        py = p0[mask] * np.sin(phi0[mask]) + p1[mask] * np.sin(phi1[mask])
        pz = p0[mask] * np.sinh(eta0[mask]) + p1[mask] * np.sinh(eta1[mask])
        out[mask] = np.sqrt(np.maximum((e0 + e1) ** 2 - px**2 - py**2 - pz**2, 0.0))
        return out
    if name == "jet_ht":
        p0, p1 = jet_pt(0), jet_pt(1)
        out = np.zeros(length, dtype=np.float64)
        valid0 = _valid(p0, sentinel=sentinel)
        valid1 = _valid(p1, sentinel=sentinel)
        out[valid0] += p0[valid0]
        out[valid1] += p1[valid1]
        out[~(valid0 | valid1)] = sentinel
        return out
    if name == "min_dphi_jet_ptmiss":
        dphi0 = _wrapped_abs_dphi(jet_phi(0), met_phi(), sentinel)
        dphi1 = _wrapped_abs_dphi(jet_phi(1), met_phi(), sentinel)
        out = np.full(length, sentinel, dtype=np.float64)
        valid0 = _valid(dphi0, sentinel=sentinel)
        valid1 = _valid(dphi1, sentinel=sentinel)
        out[valid0 & ~valid1] = dphi0[valid0 & ~valid1]
        out[valid1 & ~valid0] = dphi1[valid1 & ~valid0]
        both = valid0 & valid1
        out[both] = np.minimum(dphi0[both], dphi1[both])
        return out
    raise ValueError(f"Unsupported derived feature: {name}")


def extract_features_from_arrays(arrays, feature_specs, length: int | None = None):
    if length is None:
        fields = list(_fields(arrays))
        length = len(arrays[fields[0]]) if fields else 0
    cols: list[np.ndarray] = []
    cache: dict[tuple, np.ndarray] = {}
    for spec in feature_specs:
        sentinel = float(spec.get("missing_sentinel", MISSING_SENTINEL))
        kind = spec["kind"]
        if kind == "scalar":
            key = ("scalar", spec["branch"])
            col = cache.setdefault(key, _scalar(arrays, spec["branch"], length, sentinel))
        elif kind == "fixed_index":
            key = ("fixed", spec["branch"], int(spec["index"]))
            col = cache.setdefault(key, _fixed(arrays, spec["branch"], int(spec["index"]), length, sentinel))
        elif kind == "jagged_index":
            key = ("jagged", spec["branch"], int(spec["index"]))
            col = cache.setdefault(key, _jagged(arrays, spec["branch"], int(spec["index"]), length, sentinel))
        elif kind == "derived":
            key = ("derived", spec["name"])
            if key not in cache:
                cache[key] = _derived_feature(spec["name"], arrays, length, sentinel, cache)
            col = cache[key]
        else:
            raise ValueError(f"Unknown feature kind: {kind}")
        cols.append(_sanitize(col, sentinel).astype(np.float32))
    if not cols:
        return np.empty((length, 0), dtype=np.float32)
    return np.column_stack(cols).astype(np.float32)


def _legacy_signed_log(raw_features: np.ndarray) -> np.ndarray:
    x = np.asarray(raw_features, dtype=np.float32)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return (np.sign(x) * np.log1p(np.abs(x))).astype(np.float32)


def fit_scaler(raw_features: np.ndarray, missing_sentinel: float = MISSING_SENTINEL) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(raw_features, dtype=np.float64)
    mean = np.zeros(x.shape[1], dtype=np.float32)
    scale = np.ones(x.shape[1], dtype=np.float32)
    for index in range(x.shape[1]):
        valid = np.isfinite(x[:, index]) & (x[:, index] != float(missing_sentinel))
        if np.any(valid):
            mean[index] = float(np.mean(x[valid, index]))
            std = float(np.std(x[valid, index]))
            scale[index] = std if std > 1e-12 else 1.0
    return mean, scale


def transform_features(
    raw_features: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
    *,
    missing_sentinel: float = MISSING_SENTINEL,
    preprocessing: str = "standardize_nonmissing_then_map_missing_to_zero",
) -> np.ndarray:
    if "log1p" in str(preprocessing) or "signed_log" in str(preprocessing):
        x = _legacy_signed_log(raw_features)
        return ((x - mean) / scale).astype(np.float32)
    x = np.asarray(raw_features, dtype=np.float32)
    out = np.zeros_like(x, dtype=np.float32)
    valid = np.isfinite(x) & (x != float(missing_sentinel))
    standardized = (x - np.asarray(mean, dtype=np.float32)) / np.asarray(scale, dtype=np.float32)
    out[valid] = standardized[valid]
    return out


def save_metadata(
    path: str | Path,
    *,
    mean,
    scale,
    feature_specs,
    skipped_branches,
    config,
    training_summary,
    candidate_feature_specs=None,
):
    payload = {
        "feature_names": feature_names(feature_specs),
        "feature_specs": feature_specs,
        "candidate_feature_names": feature_names(candidate_feature_specs or feature_specs),
        "candidate_feature_specs": candidate_feature_specs or feature_specs,
        "skipped_branches": skipped_branches,
        "preprocessing": "standardize nonmissing values using train-only mean/std; map missing to zero",
        "missing_sentinel": float(config.get("missing_sentinel", MISSING_SENTINEL)),
        "mean": [float(x) for x in mean],
        "scale": [float(x) for x in scale],
        "config": config,
        "training_summary": training_summary,
    }
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True))


def load_metadata(path: str | Path):
    payload = json.loads(Path(path).read_text())
    mean = np.asarray(payload["mean"], dtype=np.float32)
    scale = np.asarray(payload["scale"], dtype=np.float32)
    return payload, mean, scale


def load_model(model_path: str | Path, device: torch.device):
    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    features = checkpoint.get("features") or checkpoint.get("feature_names")
    if not features:
        raise KeyError(f"Missing features in model checkpoint: {model_path}")
    hidden_layers = checkpoint.get("hidden_layers")
    model = SimpleDNN(
        input_size=len(features),
        hidden_size=int(checkpoint.get("hidden_size", 64)),
        dropout=float(checkpoint.get("dropout", 0.08)),
        hidden_layers=hidden_layers,
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.score_activation = str(checkpoint.get("score_activation", "softmax_signal"))
    model.to(device)
    model.eval()
    return model


def predict_scores(model, x: np.ndarray, device: torch.device, batch_size: int = 65536) -> np.ndarray:
    scores: list[np.ndarray] = []
    activation = getattr(model, "score_activation", "softmax_signal")
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            xb = torch.from_numpy(x[start : start + batch_size].astype(np.float32)).to(device)
            logits = model(xb)
            if activation == "sigmoid_signal":
                prob = torch.sigmoid(logits[:, 1])
            else:
                prob = torch.softmax(logits, dim=1)[:, 1]
            scores.append(prob.cpu().numpy())
    if not scores:
        return np.asarray([], dtype=np.float32)
    return np.concatenate(scores).astype(np.float32)
