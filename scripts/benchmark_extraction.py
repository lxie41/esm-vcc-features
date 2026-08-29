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

from protein_embeddings.esm import (ESMConfig, load_model, residue_embeddings,
                                     residue_embeddings_and_parti_attention_streaming)
from protein_embeddings.parti import (parti_attention_matrix, pagerank_weights,
                                       pagerank_weights_tensor, pool_parti)
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


def compare_weights(reference, optimized):
    """Compare normalized residue weights and their final weighted vectors."""
    rows = []
    for sequence_hash in reference:
        weights_a, vector_a = reference[sequence_hash]
        weights_b, vector_b = optimized[sequence_hash]
        midpoint = (weights_a + weights_b) / 2
        js = 0.5 * np.sum(weights_a * np.log(weights_a / midpoint)) + 0.5 * np.sum(
            weights_b * np.log(weights_b / midpoint)
        )
        rows.append({
            "sequence_hash": sequence_hash,
            "weight_cosine": float(np.dot(weights_a, weights_b) /
                                    (np.linalg.norm(weights_a) * np.linalg.norm(weights_b))),
            "weight_correlation": float(np.corrcoef(weights_a, weights_b)[0, 1]),
            "weight_max_abs_error": float(np.max(np.abs(weights_a - weights_b))),
            "weight_js_distance": float(np.sqrt(max(js, 0.0))),
            "vector_cosine": float(np.dot(vector_a, vector_b) /
                                    (np.linalg.norm(vector_a) * np.linalg.norm(vector_b))),
            "vector_max_abs_error": float(np.max(np.abs(vector_a - vector_b))),
            "finite": bool(np.isfinite(weights_b).all() and np.isfinite(vector_b).all()),
        })
    return {
        "mean_weight_cosine": float(np.mean([r["weight_cosine"] for r in rows])),
        "min_weight_cosine": min(r["weight_cosine"] for r in rows),
        "mean_weight_correlation": float(np.mean([r["weight_correlation"] for r in rows])),
        "max_weight_abs_error": max(r["weight_max_abs_error"] for r in rows),
        "max_weight_js_distance": max(r["weight_js_distance"] for r in rows),
        "mean_vector_cosine": float(np.mean([r["vector_cosine"] for r in rows])),
        "max_vector_abs_error": max(r["vector_max_abs_error"] for r in rows),
        "all_finite": all(r["finite"] for r in rows),
        "identical_sequence_hash_order": list(reference) == list(optimized),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--count", type=int, default=200)
    ap.add_argument("--max-tokens", type=int, action="append", default=None,
                    help="repeatable padded-token budgets; explicit values replace defaults")
    ap.add_argument("--max-batch-size", type=int, default=16)
    ap.add_argument("--profile", action="store_true",
                    help="include serial stage timing breakdown")
    ap.add_argument("--pagerank-backend", choices=["networkx", "tensor"],
                    default="networkx")
    ap.add_argument("--compare-pagerank", action="store_true",
                    help="run both NetworkX reference and tensor PageRank on the same subset")
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

    profile_timings = {} if args.profile else None

    def serial_fn():
        return {r.sequence_hash: extract_one(model, alphabet, r.amino_acid_sequence, device, cfg,
                                              pagerank_backend=args.pagerank_backend)
                if profile_timings is None else extract_one(
                    model, alphabet, r.amino_acid_sequence, device, cfg,
                    timings=profile_timings, pagerank_backend=args.pagerank_backend)
                for r in native}
    serial, serial_seconds, serial_peak = timed(serial_fn, device)
    results = {"subset_count": len(subset), "native_count": len(native),
               "lengths": subset.protein_length.tolist(),
               "serial": {"seconds": serial_seconds,
                          "proteins_per_hour": len(native) / serial_seconds * 3600,
                          "residues_per_second": float(sum(len(r.amino_acid_sequence) for r in native) / serial_seconds),
                          "peak_vram_gib": serial_peak}}
    if profile_timings is not None:
        results["serial"]["stage_seconds"] = profile_timings

    if args.compare_pagerank:
        def reference_fn():
            return {r.sequence_hash: extract_one(model, alphabet, r.amino_acid_sequence,
                                                  device, cfg, pagerank_backend="networkx")
                    for r in native}
        reference, ref_seconds, ref_peak = timed(reference_fn, device)
        def tensor_fn():
            return {r.sequence_hash: extract_one(model, alphabet, r.amino_acid_sequence,
                                                  device, cfg, pagerank_backend="tensor")
                    for r in native}
        optimized, opt_seconds, opt_peak = timed(tensor_fn, device)
        results["pagerank_backend_comparison"] = {
            "reference": {"backend": "networkx", "seconds": ref_seconds,
                          "peak_vram_gib": ref_peak},
            "optimized": {"backend": "tensor", "seconds": opt_seconds,
                           "peak_vram_gib": opt_peak,
                           "speedup": ref_seconds / opt_seconds},
            "comparison": compare(reference, optimized),
        }
        def weight_fn(backend):
            output = {}
            for row in native:
                states, attention = residue_embeddings_and_parti_attention_streaming(
                    model, alphabet, row.amino_acid_sequence, device, torch.float16)
                matrix = parti_attention_matrix(attention.unsqueeze(0).unsqueeze(1))
                weights = (pagerank_weights(matrix) if backend == "networkx" else
                           pagerank_weights_tensor(matrix))
                output[row.sequence_hash] = (weights.numpy(),
                                             (states * weights[:, None]).sum(0).numpy())
            return {key: output[key] for key in sorted(output)}
        weight_reference, _, _ = timed(lambda: weight_fn("networkx"), device)
        weight_optimized, _, _ = timed(lambda: weight_fn("tensor"), device)
        results["pagerank_weight_comparison"] = compare_weights(weight_reference, weight_optimized)

    def mean_sd_only_fn():
        output = {}
        for row in native:
            H = residue_embeddings(model, alphabet, row.amino_acid_sequence,
                                    cfg.layer, device,
                                    torch.float16 if device.startswith("cuda") else None)
            output[row.sequence_hash] = (mean_representation(H).numpy(),
                                         mean_sd_representation(H).numpy())
        return output
    _, msd_seconds, msd_peak = timed(mean_sd_only_fn, device)
    native_residues = sum(len(r.amino_acid_sequence) for r in native)
    results["mean_sd_only"] = {
        "seconds": msd_seconds,
        "proteins_per_hour": len(native) / msd_seconds * 3600,
        "residues_per_second": float(native_residues / msd_seconds),
        "peak_vram_gib": msd_peak,
        "relative_to_mean_sd_parti_serial": serial_seconds / msd_seconds,
    }
    for budget in budgets:
        batches = make_length_batches(native, budget, args.max_batch_size)
        def batch_fn():
            output = {}
            for batch in batches:
                values = extract_native_batch(model, alphabet,
                                              [r.amino_acid_sequence for r in batch], device, cfg,
                                              pagerank_backend=args.pagerank_backend)
                for row, value in zip(batch, values):
                    output[row.sequence_hash] = value
            # The execution order is length-sorted, but the comparison/output
            # contract is deterministic sequence_hash order.
            return {key: output[key] for key in sorted(output)}
        batched, seconds, peak = timed(batch_fn, device)
        results[f"batched_{budget}"] = {
            "seconds": seconds, "proteins_per_hour": len(native) / seconds * 3600,
            "residues_per_second": float(native_residues / seconds),
            "peak_vram_gib": peak, "batch_count": len(batches),
            "max_batch_size": max(len(b) for b in batches),
            "max_padded_tokens_per_batch": max(
                (max(len(r.amino_acid_sequence) for r in b) + 2) * len(b)
                for b in batches
            ),
            "speedup_vs_serial": serial_seconds / seconds,
            "pagerank_backend": args.pagerank_backend,
            "comparison": compare(serial, batched),
        }
    Path(args.output).write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
