"""Benchmark serial versus length-aware native-context extraction on a subset.

This is a pilot benchmark only. It never discovers or processes the full
production universe unless the caller explicitly supplies such an input.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from protein_embeddings.esm import ESMConfig, load_model
from protein_embeddings.parti import pool_parti
from protein_embeddings.pooling import mean_representation, mean_sd_representation
from scripts.extract_esm_features import extract_native_batch, extract_one, make_length_batches


def select_stratified(df, count):
    bins = [(0, 350), (450, 550), (650, 750), (850, 950), (980, 1022)]
    selected = []
    per_bin = max(1, count // len(bins))
    for low, high in bins:
        part = df[(df.protein_length >= low) & (df.protein_length <= high)]
        selected.append(part.sort_values("sequence_hash").iloc[:per_bin])
    chosen = pd.concat(selected).drop_duplicates("sequence_hash")
    if len(chosen) < count:
        rest = df[~df.sequence_hash.isin(chosen.sequence_hash)].sort_values("sequence_hash")
        chosen = pd.concat([chosen, rest.iloc[:count - len(chosen)]])
    return chosen.sort_values("sequence_hash").iloc[:count].reset_index(drop=True)


def timed(fn, device):
    if device.startswith("cuda"):
        torch.cuda.synchronize()
        torch.cuda.reset_peak_memory_stats()
    start = time.perf_counter()
    result = fn()
    if device.startswith("cuda"):
        torch.cuda.synchronize()
        peak = torch.cuda.max_memory_allocated() / 1024**3
    else:
        peak = None
    return result, time.perf_counter() - start, peak


def compare(serial, batched):
    output = {}
    for feature_index, feature in enumerate(("mean", "mean_sd", "parti")):
        values = []
        for key in serial:
            a, b = serial[key][feature_index], batched[key][feature_index]
            values.append({"sequence_hash": key,
                           "cosine": float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))),
                           "max_abs_error": float(np.max(np.abs(a - b))),
                           "finite": bool(np.isfinite(b).all())})
        output[feature] = {"min_cosine": min(x["cosine"] for x in values),
                           "mean_cosine": float(np.mean([x["cosine"] for x in values])),
                           "max_abs_error": max(x["max_abs_error"] for x in values),
                           "all_finite": all(x["finite"] for x in values)}
    output["identical_sequence_hash_order"] = list(serial) == list(batched)
    return output


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--count", type=int, default=200)
    ap.add_argument("--max-tokens", type=int, action="append", default=None)
    ap.add_argument("--max-batch-size", type=int, default=16)
    args = ap.parse_args()
    budgets = args.max_tokens or [4096, 8192, 16384]

    df = pd.read_parquet(args.input)
    if "protein_length" not in df:
        df["protein_length"] = df.amino_acid_sequence.str.len()
    df = df.drop_duplicates("sequence_hash").sort_values("sequence_hash")
    subset = select_stratified(df, min(args.count, len(df)))
    model, alphabet, device = load_model(args.model, args.device)
    cfg = ESMConfig()
    native = [r for r in subset.itertuples(index=False) if len(r.amino_acid_sequence) <= cfg.window_size]

    def serial_fn():
        return {r.sequence_hash: extract_one(model, alphabet, r.amino_acid_sequence, device, cfg)
                for r in native}
    serial, serial_seconds, serial_peak = timed(serial_fn, device)
    results = {"subset_count": len(subset), "native_count": len(native),
               "lengths": subset.protein_length.tolist(),
               "serial": {"seconds": serial_seconds,
                          "proteins_per_hour": len(native) / serial_seconds * 3600,
                          "residues_per_second": float(subset.protein_length.sum() / serial_seconds),
                          "peak_vram_gib": serial_peak}}

    for budget in budgets:
        batches = make_length_batches(native, budget, args.max_batch_size)
        def batch_fn():
            output = {}
            for batch in batches:
                values = extract_native_batch(model, alphabet,
                                              [r.amino_acid_sequence for r in batch], device, cfg)
                for row, value in zip(batch, values):
                    output[row.sequence_hash] = value
            # The execution order is length-sorted, but the comparison/output
            # contract is deterministic sequence_hash order.
            return {key: output[key] for key in sorted(output)}
        batched, seconds, peak = timed(batch_fn, device)
        results[f"batched_{budget}"] = {
            "seconds": seconds, "proteins_per_hour": len(native) / seconds * 3600,
            "residues_per_second": float(subset.protein_length.sum() / seconds),
            "peak_vram_gib": peak, "batch_count": len(batches),
            "max_batch_size": max(len(b) for b in batches),
            "max_padded_tokens_per_batch": max(
                (max(len(r.amino_acid_sequence) for r in b) + 2) * len(b)
                for b in batches
            ),
            "speedup_vs_serial": serial_seconds / seconds,
            "comparison": compare(serial, batched),
        }
    Path(args.output).write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
