#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
import uproot
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset

from dnn_common import (
    BKG_PROCESSES,
    DEFAULT_CONFIG_PATH,
    MISSING_SENTINEL,
    SimpleDNN,
    available_branches,
    branches_for_specs,
    discover_feature_specs,
    extract_features_from_arrays,
    feature_names,
    fit_scaler,
    load_config,
    predict_scores,
    save_metadata,
    transform_features,
)


def read_process(path, label, process, feature_specs, max_events=None, rng=None):
    requested = branches_for_specs(feature_specs, extra=["weight", "lepton_cat", "ptmiss"])
    with uproot.open(path)["Vars"] as tree:
        # A cap is primarily a debugging option. Read a modest oversample so
        # the SR selection can still return the requested number of events.
        entry_stop = None if max_events is None else min(int(tree.num_entries), max(5 * int(max_events), int(max_events)))
        arrays = tree.arrays(available_branches(tree, requested), entry_stop=entry_stop, library="ak")
        n_entries = len(arrays["weight"])

    raw_features = extract_features_from_arrays(arrays, feature_specs, length=n_entries)
    lepton_cat = np.asarray(arrays["lepton_cat"])
    ptmiss = np.asarray(arrays["ptmiss"], dtype=np.float32)
    w_signed = np.asarray(arrays["weight"], dtype=np.float32)
    w_signed = np.where(np.isfinite(w_signed), w_signed, 0.0).astype(np.float32)
    mask = (lepton_cat != 2) & (ptmiss >= 100.0) & np.isfinite(ptmiss)
    x = raw_features[mask].astype(np.float32)
    y = np.full(len(x), label, dtype=np.int64)
    w_signed = w_signed[mask]
    proc = np.full(len(x), process)
    if max_events is not None and len(x) > max_events:
        if rng is None:
            rng = np.random.default_rng(12345)
        keep = rng.choice(len(x), size=max_events, replace=False)
        x, y, w_signed, proc = x[keep], y[keep], w_signed[keep], proc[keep]
    return x, y, w_signed, proc


def load_training_arrays(input_bases, signal_file, feature_specs, max_events_per_process, rng):
    pieces = []
    for source_index, input_base in enumerate(input_bases):
        source = Path(input_base).name or f"source{source_index}"
        sig = read_process(
            os.path.join(input_base, "SR", signal_file), 1, f"{source}:signal",
            feature_specs, max_events_per_process, rng,
        )
        print(f"Loaded {source}:signal: {len(sig[1])} events", flush=True)
        pieces.append(sig)
        for process in BKG_PROCESSES:
            piece = read_process(
                os.path.join(input_base, "SR", f"{process}.root"), 0, f"{source}:{process}",
                feature_specs, max_events_per_process, rng,
            )
            print(f"Loaded {source}:{process}: {len(piece[1])} events", flush=True)
            pieces.append(piece)
    x = np.concatenate([piece[0] for piece in pieces])
    y = np.concatenate([piece[1] for piece in pieces])
    w_signed = np.concatenate([piece[2] for piece in pieces])
    process = np.concatenate([piece[3] for piece in pieces])
    return x, y, w_signed, process


def nonnegative_clipped_weights(w_signed, clip):
    w = np.where(np.isfinite(w_signed), w_signed, 0.0).astype(np.float64)
    w = np.maximum(w, 0.0)
    if clip is not None and clip > 0:
        w = np.minimum(w, float(clip))
    return w.astype(np.float32)


def _weighted_percentile(values, weights, quantiles):
    sorter = np.argsort(values)
    values = values[sorter]
    weights = weights[sorter]
    cdf = np.cumsum(weights)
    if not len(cdf) or cdf[-1] <= 0:
        return np.percentile(values, quantiles * 100.0)
    return np.interp(quantiles, cdf / cdf[-1], values)


def _asimov_z(sig, bkg):
    s = np.maximum(np.asarray(sig, dtype=np.float64), 0.0)
    b = np.maximum(np.asarray(bkg, dtype=np.float64), 0.0)
    term = np.where(
        b > 1e-12,
        (s + b) * np.log1p(np.divide(s, b, out=np.zeros_like(s), where=b > 1e-12)) - s,
        0.0,
    )
    return float(np.sqrt(2.0 * np.sum(np.maximum(term, 0.0))))


def _asimov_z_syst(sig, bkg, sigma_rel):
    if sigma_rel <= 0.0:
        return _asimov_z(sig, bkg)
    s = np.maximum(np.asarray(sig, dtype=np.float64), 0.0)
    b = np.maximum(np.asarray(bkg, dtype=np.float64), 0.0)
    sigma_b2 = (float(sigma_rel) * b) ** 2
    term = np.zeros_like(s)
    mask = (b > 1e-12) & (sigma_b2 > 1e-12)
    if np.any(mask):
        sm, bm, sb2 = s[mask], b[mask], sigma_b2[mask]
        numerator = (sm + bm) * (bm + sb2)
        denominator = bm * bm + (sm + bm) * sb2
        ratio1 = np.divide(numerator, denominator, out=np.ones_like(sm), where=denominator > 1e-12)
        term1 = np.where(ratio1 > 1.0 + 1e-12, (sm + bm) * np.log(ratio1), 0.0)
        ratio2 = sb2 * sm / np.maximum(bm * (bm + sb2), 1e-12)
        term2 = (bm * bm / np.maximum(sb2, 1e-12)) * np.log1p(ratio2)
        term[mask] = term1 - term2
    pure = ~mask & (b > 1e-12)
    if np.any(pure):
        term[pure] = (s[pure] + b[pure]) * np.log1p(s[pure] / b[pure]) - s[pure]
    return float(np.sqrt(2.0 * np.sum(np.maximum(term, 0.0))))


def compute_feature_significance(raw_x, y, weights, names, n_bins, sig_syst, outdir):
    """Tiwari-style full-dataset ranking; intentionally runs before the split."""
    rows = []
    y = np.asarray(y, dtype=np.int64)
    weights = np.maximum(np.asarray(weights, dtype=np.float64), 0.0)
    total_s = float(weights[y == 1].sum())
    total_b = float(weights[y == 0].sum())
    baseline_z = _asimov_z([total_s], [total_b])
    baseline_z_syst = _asimov_z_syst([total_s], [total_b], sig_syst)

    for index, name in enumerate(names):
        values = np.asarray(raw_x[:, index], dtype=np.float64)
        finite = np.isfinite(values) & np.isfinite(weights)
        x, yy, ww = values[finite], y[finite], weights[finite]
        if not len(x) or ww[yy == 0].sum() <= 0 or ww[yy == 1].sum() <= 0:
            rows.append({"feature": name, "auc": None, "asimov_z": 0.0, "asimov_z_syst": 0.0})
            continue
        lo, hi = _weighted_percentile(x, ww, np.asarray([0.01, 0.99]))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = float(np.min(x)), float(np.max(x))
        if hi <= lo:
            hi = lo + 1.0
        edges = np.linspace(float(lo), float(hi), int(n_bins) + 1)
        hs, _ = np.histogram(x[yy == 1], bins=edges, weights=ww[yy == 1])
        hb, _ = np.histogram(x[yy == 0], bins=edges, weights=ww[yy == 0])
        auc = float(roc_auc_score(yy, x, sample_weight=ww))
        z = _asimov_z(hs, hb)
        z_syst = _asimov_z_syst(hs, hb, sig_syst)
        rows.append({
            "feature": name,
            "auc": auc,
            "asimov_z": z,
            "asimov_z_syst": z_syst,
            "delta_z": z - baseline_z,
            "delta_z_syst": z_syst - baseline_z_syst,
            "range_low_weighted_p01": float(lo),
            "range_high_weighted_p99": float(hi),
            "n_bins": int(n_bins),
            "background_systematic_fraction": float(sig_syst),
        })
    rows.sort(
        key=lambda row: (
            float(row["asimov_z_syst"]),
            abs(float(row["auc"] if row["auc"] is not None else 0.5) - 0.5),
        ),
        reverse=True,
    )
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        row["selection_scope"] = "full_dataset_before_split"
        row["known_leakage_issue"] = True

    (outdir / "feature_significance.json").write_text(json.dumps(rows, indent=2) + "\n")
    with open(outdir / "feature_significance.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def class_balance_factors(y_train, w_train):
    sum_b = float(w_train[y_train == 0].sum())
    sum_s = float(w_train[y_train == 1].sum())
    if sum_b <= 0 or sum_s <= 0:
        raise ValueError("Cannot balance classes: one training class has zero weight")
    factor_b = 0.5 / sum_b
    factor_s = 0.5 / sum_s
    scale = (sum_b + sum_s) / (sum_b * factor_b + sum_s * factor_s)
    return float(factor_b * scale), float(factor_s * scale)


def apply_balance(y, weights, factor_b, factor_s):
    return (weights * np.where(y == 1, factor_s, factor_b)).astype(np.float32)


def make_loader(x, y, weights, batch_size, shuffle):
    dataset = TensorDataset(
        torch.from_numpy(x.astype(np.float32)),
        torch.from_numpy(y.astype(np.int64)),
        torch.from_numpy(weights.astype(np.float32)),
    )
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=0)


def weighted_two_logit_bce(logits, labels, weights):
    targets = F.one_hot(labels, num_classes=2).to(dtype=logits.dtype)
    per_event = F.binary_cross_entropy_with_logits(logits, targets, reduction="none").mean(dim=1)
    return (per_event * weights).sum() / torch.clamp(weights.sum(), min=1e-12)


def evaluate(model, x, y, loss_weights, metric_weights, device, batch_size):
    model.eval()
    losses, logits_all = [], []
    with torch.no_grad():
        for xb, yb, wb in make_loader(x, y, loss_weights, batch_size, False):
            xb, yb, wb = xb.to(device), yb.to(device), wb.to(device)
            logits = model(xb)
            losses.append(float(weighted_two_logit_bce(logits, yb, wb).cpu()))
            logits_all.append(logits.cpu())
    logits = torch.cat(logits_all, dim=0)
    score = torch.sigmoid(logits[:, 1]).numpy()
    auc = roc_auc_score(y, score, sample_weight=metric_weights)
    return float(np.mean(losses)), float(auc), score


def parse_hidden_layers(text):
    layers = [int(value.strip()) for value in str(text).split(",") if value.strip()]
    if not layers or any(value <= 0 for value in layers):
        raise ValueError(f"Invalid hidden layers: {text}")
    return layers


def parse_args():
    parser = argparse.ArgumentParser(description="Train the weighted HZZ2l2nu two-logit DNN score.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--input-base", action="append", default=None, help="Repeat for combined-year training")
    parser.add_argument("--signal-file", default=None)
    parser.add_argument("--outdir", default="/afs/cern.ch/user/l/liwe/hzz2l2nu/dnn/output")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--min-delta", type=float, default=1e-6)
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--hidden-layers", default="128,64,32")
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-clip", type=float, default=None)
    parser.add_argument("--top-k", type=int, default=None)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--max-events-per-process", type=int, default=0)
    parser.add_argument("--torch-threads", type=int, default=4)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    configured_inputs = config.get("input_bases", config.get("input_base"))
    if args.input_base:
        input_bases = args.input_base
    elif isinstance(configured_inputs, list):
        input_bases = configured_inputs
    else:
        input_bases = [configured_inputs]
    signal_file = args.signal_file or config["signal_file"]
    weight_clip = float(args.weight_clip if args.weight_clip is not None else config.get("weight_clip", 100.0))
    top_k = int(args.top_k if args.top_k is not None else config.get("top_k_significance", 10))
    missing_sentinel = float(config.get("missing_sentinel", MISSING_SENTINEL))
    n_significance_bins = int(config.get("feature_significance_bins", 40))
    sig_syst = float(config.get("feature_significance_background_syst", 0.20))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    torch.set_num_threads(args.torch_threads)
    rng = np.random.default_rng(args.seed)

    reference_file = os.path.join(input_bases[0], "SR", signal_file)
    candidate_specs, skipped_branches = discover_feature_specs(reference_file, config)
    candidate_names = feature_names(candidate_specs)
    print(f"Explicit HZZ candidate features ({len(candidate_names)}):", flush=True)
    for index, name in enumerate(candidate_names, start=1):
        print(f"  {index:02d}. {name}", flush=True)

    max_events = args.max_events_per_process if args.max_events_per_process > 0 else None
    raw_x_all, y, w_signed, process = load_training_arrays(
        input_bases, signal_file, candidate_specs, max_events, rng,
    )
    w_physical = nonnegative_clipped_weights(w_signed, weight_clip)
    if float(w_physical.sum()) <= 0:
        raise ValueError("All training/ranking weights are zero after sanitation and clipping")

    # Intentionally match Tiwari: rank/select on the full dataset before splitting.
    significance_rows = compute_feature_significance(
        raw_x_all, y, w_physical, candidate_names,
        n_significance_bins, sig_syst, outdir,
    )
    ranked_names = [row["feature"] for row in significance_rows]
    selected_names = ranked_names[: min(top_k, len(ranked_names))] if top_k > 0 else ranked_names
    selected_indices = [candidate_names.index(name) for name in selected_names]
    selected_specs = [candidate_specs[index] for index in selected_indices]
    raw_x = raw_x_all[:, selected_indices]
    print(f"Selected top-{len(selected_names)} features: {selected_names}", flush=True)
    print("[KNOWN ISSUE] Feature selection used the full dataset before splitting (test leakage).", flush=True)

    all_indices = np.arange(len(y))
    train_idx, temp_idx = train_test_split(
        all_indices, test_size=0.50, random_state=args.seed, stratify=y,
    )
    val_idx, test_idx = train_test_split(
        temp_idx, test_size=0.60, random_state=args.seed, stratify=y[temp_idx],
    )
    print(f"Split: train={len(train_idx)} (50%), val={len(val_idx)} (20%), test={len(test_idx)} (30%)", flush=True)

    mean, scale = fit_scaler(raw_x[train_idx], missing_sentinel=missing_sentinel)
    x = transform_features(
        raw_x, mean, scale, missing_sentinel=missing_sentinel,
        preprocessing="standardize_nonmissing_then_map_missing_to_zero",
    )

    factor_b, factor_s = class_balance_factors(y[train_idx], w_physical[train_idx])
    w_train_eff = apply_balance(y[train_idx], w_physical[train_idx], factor_b, factor_s)
    w_val_eff = apply_balance(y[val_idx], w_physical[val_idx], factor_b, factor_s)
    w_test_eff = apply_balance(y[test_idx], w_physical[test_idx], factor_b, factor_s)
    print(f"Train-only class balance factors: background={factor_b:.6g}, signal={factor_s:.6g}", flush=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    hidden_layers = parse_hidden_layers(args.hidden_layers)
    model = SimpleDNN(
        input_size=len(selected_names), hidden_layers=hidden_layers, dropout=args.dropout,
    ).to(device)
    model.score_activation = "sigmoid_signal"
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.learning_rate, weight_decay=1e-4)
    train_loader = make_loader(x[train_idx], y[train_idx], w_train_eff, args.batch_size, True)

    best_val_auc = -np.inf
    best_state = None
    bad_epochs = 0
    history = []
    for epoch in range(1, args.epochs + 1):
        model.train()
        train_losses = []
        for xb, yb, wb in train_loader:
            xb, yb, wb = xb.to(device), yb.to(device), wb.to(device)
            optimizer.zero_grad(set_to_none=True)
            loss = weighted_two_logit_bce(model(xb), yb, wb)
            loss.backward()
            optimizer.step()
            train_losses.append(float(loss.detach().cpu()))
        val_loss, val_auc, _ = evaluate(
            model, x[val_idx], y[val_idx], w_val_eff, w_physical[val_idx], device, args.batch_size,
        )
        train_loss = float(np.mean(train_losses))
        history.append((epoch, train_loss, val_loss, val_auc))
        print(f"epoch {epoch:03d}: train_loss={train_loss:.6f} val_loss={val_loss:.6f} val_auc={val_auc:.6f}", flush=True)
        if val_auc > best_val_auc + args.min_delta:
            best_val_auc = val_auc
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            bad_epochs = 0
        else:
            bad_epochs += 1
            if args.patience > 0 and bad_epochs >= args.patience:
                print(f"Early stopping: validation AUC did not improve for {bad_epochs} epochs", flush=True)
                break

    if best_state is None:
        raise RuntimeError("Training did not produce a valid checkpoint")
    model.load_state_dict(best_state)
    test_loss, test_auc, test_score = evaluate(
        model, x[test_idx], y[test_idx], w_test_eff, w_physical[test_idx], device, args.batch_size,
    )
    _, train_auc, train_score = evaluate(
        model, x[train_idx], y[train_idx], w_train_eff, w_physical[train_idx], device, args.batch_size,
    )
    _, val_auc, val_score = evaluate(
        model, x[val_idx], y[val_idx], w_val_eff, w_physical[val_idx], device, args.batch_size,
    )
    fpr, tpr, thresholds = roc_curve(y[test_idx], test_score, sample_weight=w_physical[test_idx])

    np.savetxt(
        outdir / "training_history.csv", np.asarray(history), delimiter=",",
        header="epoch,train_loss,val_loss,val_auc", comments="",
    )
    np.savetxt(
        outdir / "test_roc.csv", np.column_stack([fpr, tpr, thresholds]), delimiter=",",
        header="fpr,tpr,threshold", comments="",
    )
    np.savez(
        outdir / "test_scores.npz",
        score=test_score,
        label=y[test_idx],
        train_weight=w_test_eff,
        physical_nonnegative_weight=w_physical[test_idx],
        signed_weight=w_signed[test_idx],
        raw_features=raw_x[test_idx],
        transformed_features=x[test_idx],
        process=process[test_idx],
    )
    np.savez(
        outdir / "split_scores.npz",
        train_score=train_score, train_label=y[train_idx], train_weight=w_physical[train_idx],
        val_score=val_score, val_label=y[val_idx], val_weight=w_physical[val_idx],
        test_score=test_score, test_label=y[test_idx], test_weight=w_physical[test_idx],
    )

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "hidden_layers": hidden_layers,
            "dropout": args.dropout,
            "features": selected_names,
            "feature_specs": selected_specs,
            "candidate_features": candidate_names,
            "output_nodes": 2,
            "loss_function": "weighted_BCEWithLogitsLoss_two_output_one_hot",
            "score_activation": "sigmoid_signal",
            "best_val_auc": float(best_val_auc),
            "test_auc": float(test_auc),
        },
        outdir / "dnn_model.pt",
    )

    missing_fraction = {
        name: float(np.mean(raw_x_all[:, index] == missing_sentinel))
        for index, name in enumerate(candidate_names)
    }
    training_summary = {
        "best_val_auc": float(best_val_auc),
        "train_auc": float(train_auc),
        "validation_auc": float(val_auc),
        "test_auc": float(test_auc),
        "test_loss": float(test_loss),
        "candidate_features": len(candidate_names),
        "selected_features": len(selected_names),
        "selected_feature_names": selected_names,
        "num_training_candidates": int(len(y)),
        "n_train": int(len(train_idx)),
        "n_validation": int(len(val_idx)),
        "n_test": int(len(test_idx)),
        "epochs_completed": int(len(history)),
        "early_stopped": bool(len(history) < args.epochs),
        "class_balance_factors_train_only": {"background": factor_b, "signal": factor_s},
        "negative_weight_events_ignored_for_training": int(np.sum(w_signed < 0)),
        "weights_clipped_at_max": int(np.sum(w_signed > weight_clip)),
        "weight_clip": float(weight_clip),
        "candidate_missing_fractions": missing_fraction,
        "feature_selection_scope": "full_dataset_before_split",
        "feature_selection_known_test_leakage": True,
        "max_events_per_process": max_events,
    }
    resolved_config = {
        **config,
        "input_bases": input_bases,
        "signal_file": signal_file,
        "epochs": args.epochs,
        "patience": args.patience,
        "min_delta": args.min_delta,
        "batch_size": args.batch_size,
        "hidden_layers": hidden_layers,
        "dropout": args.dropout,
        "learning_rate": args.learning_rate,
        "weight_clip": weight_clip,
        "top_k_significance": top_k,
        "loss_function": "weighted BCEWithLogitsLoss over two one-hot output targets",
        "score_activation": "sigmoid of signal logit; outputs are not softmax-normalized",
        "template_weight": "original signed ROOT weight remains unchanged for Combine templates",
    }
    save_metadata(
        outdir / "dnn_metadata.json",
        mean=mean,
        scale=scale,
        feature_specs=selected_specs,
        candidate_feature_specs=candidate_specs,
        skipped_branches=skipped_branches,
        config=resolved_config,
        training_summary=training_summary,
    )

    with open(outdir / "training_sample_summary.txt", "w", encoding="utf-8") as handle:
        handle.write(f"candidate features ({len(candidate_names)}): {', '.join(candidate_names)}\n")
        handle.write(f"selected top-{len(selected_names)}: {', '.join(selected_names)}\n")
        handle.write("preprocessing: train-only standardization; missing sentinel -> zero after scaling\n")
        handle.write("loss: weighted two-output BCEWithLogitsLoss with one-hot labels\n")
        handle.write("KNOWN ISSUE: feature significance/top-K used all events before split (test leakage)\n")
        handle.write(f"train_auc: {train_auc:.6f}\nvalidation_auc: {val_auc:.6f}\ntest_auc: {test_auc:.6f}\n\n")
        for proc_name in sorted(set(process)):
            mask = process == proc_name
            handle.write(
                f"{proc_name:20s} entries={mask.sum():9d} sumw={w_signed[mask].sum():14.6g} "
                f"sum_trainw={w_physical[mask].sum():14.6g} negative={(w_signed[mask] < 0).sum():8d}\n"
            )

    print(f"Saved model and diagnostics in {outdir}", flush=True)
    print(
        f"AUC train={train_auc:.6f}, validation={val_auc:.6f}, test={test_auc:.6f}; "
        f"epochs={len(history)}; device={device}",
        flush=True,
    )


if __name__ == "__main__":
    main()
