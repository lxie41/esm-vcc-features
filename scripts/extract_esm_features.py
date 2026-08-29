"""Resumable hash-keyed production extractor; invocation only, never auto-run."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from protein_embeddings.esm import (
    ESMConfig, chunk_starts, load_model, reconstruction_weights,
    residue_embeddings_and_parti_attention_streaming,
    residue_embeddings_and_parti_attention_streaming_batch,
    reconstruct_residue_embeddings,
)
from protein_embeddings.parti import pool_parti
from protein_embeddings.pooling import mean_representation, mean_sd_representation
from protein_embeddings.storage import update_shard_checkpoint


def _autocast_dtype(device: str):
    return torch.float16 if device.startswith("cuda") else None


def extract_one(model, alphabet, seq, device, cfg, timings=None,
                pagerank_backend="networkx"):
    """Extract one sequence, retaining the validated long-protein path."""
    starts = chunk_starts(len(seq), cfg.window_size, cfg.overlap)
    chunks, local = [], []
    for start in starts:
        h, attention = residue_embeddings_and_parti_attention_streaming(
            model, alphabet, seq[start:start + cfg.window_size], device,
            _autocast_dtype(device), timings=timings,
        )
        chunks.append(h)
        z, weights = pool_parti(h, attention.unsqueeze(0).unsqueeze(1), timings=timings,
                                pagerank_backend=pagerank_backend)
        local.append((z, weights))
    H = reconstruct_residue_embeddings(
        chunks, starts, len(seq), cfg.overlap, cfg.weighting
    )
    if len(starts) == 1:
        parti = local[0][0]
    else:
        weights = torch.zeros(len(seq), dtype=torch.float32)
        for h, (start, (_, local_weights)) in zip(chunks, zip(starts, local)):
            taper = reconstruction_weights(len(h), cfg.overlap, cfg.weighting)
            weights[start:start + len(h)] += local_weights * taper
        weights /= weights.sum()
        parti = (H * weights[:, None]).sum(0)
    start = time.perf_counter()
    mean = mean_representation(H).numpy()
    if timings is not None:
        timings["mean_pooling_seconds"] = timings.get("mean_pooling_seconds", 0.0) + time.perf_counter() - start
    start = time.perf_counter()
    mean_sd = mean_sd_representation(H).numpy()
    if timings is not None:
        timings["sd_pooling_seconds"] = timings.get("sd_pooling_seconds", 0.0) + time.perf_counter() - start
    return (mean, mean_sd,
            parti.numpy(), len(starts))


def extract_native_batch(model, alphabet, sequences, device, cfg,
                         pagerank_backend="networkx"):
    """Extract native-context sequences in one padded model forward."""
    pairs = residue_embeddings_and_parti_attention_streaming_batch(
        model, alphabet, sequences, device, _autocast_dtype(device)
    )
    results = []
    for H, attention in pairs:
        parti, weights = pool_parti(H, attention.unsqueeze(0).unsqueeze(1),
                                    pagerank_backend=pagerank_backend)
        if not torch.isfinite(weights).all() or not torch.isfinite(parti).all():
            raise FloatingPointError("non-finite PaRTI output")
        results.append((mean_representation(H).numpy(),
                        mean_sd_representation(H).numpy(), parti.numpy(), 1))
    return results


def make_length_batches(rows, max_tokens: int, max_batch_size: int):
    """Greedily batch adjacent length-sorted rows under a token budget."""
    if max_tokens <= 0 or max_batch_size <= 0:
        raise ValueError("max_tokens and max_batch_size must be positive")
    ordered = sorted(rows, key=lambda r: (len(r.amino_acid_sequence), r.sequence_hash))
    batches, current = [], []
    for row in ordered:
        padded_tokens = (max(len(r.amino_acid_sequence) for r in current + [row]) + 2) * (len(current) + 1)
        if current and (len(current) >= max_batch_size or padded_tokens > max_tokens):
            batches.append(current)
            current = []
        current.append(row)
    if current:
        batches.append(current)
    return batches


def _atomic_save_npy(path: Path, array: np.ndarray):
    tmp = path.with_name(path.name + ".tmp.npy")
    np.save(tmp, np.asarray(array, dtype=np.float32))
    tmp.replace(path)


def _load_existing(out: Path, features):
    """Load only complete, finite prior rows for safe resume."""
    metadata_path = out / "metadata.parquet"
    if not metadata_path.exists():
        return {}, {}, []
    metadata = pd.read_parquet(metadata_path).reset_index(drop=True)
    vectors = {}
    for feature in features:
        path = out / f"{feature}.npy"
        if not path.exists():
            return {}, {}, []
        vectors[feature] = np.load(path, mmap_mode="r")
    complete = {}
    for index, row in metadata.iterrows():
        values = {feature: np.asarray(vectors[feature][index], dtype=np.float32)
                  for feature in features}
        if all(np.isfinite(value).all() for value in values.values()):
            complete[row.sequence_hash] = (values, row.to_dict())
    return complete, vectors, metadata


def _write_outputs(out: Path, records, features):
    """Write hash-sorted successful rows with atomically replaced arrays."""
    keys = sorted(records)
    metadata = pd.DataFrame([records[key][1] for key in keys])
    tmp_meta = out / "metadata.parquet.tmp"
    metadata.to_parquet(tmp_meta, index=False)
    tmp_meta.replace(out / "metadata.parquet")
    for feature in features:
        matrix = np.stack([records[key][0][feature] for key in keys]).astype(np.float32)
        _atomic_save_npy(out / f"{feature}.npy", matrix)


def _record(values, row, chunks):
    mean, mean_sd, parti, _ = values
    return ({"mean": np.asarray(mean, dtype=np.float32),
             "std": np.asarray(mean_sd[1280:], dtype=np.float32),
             "parti": np.asarray(parti, dtype=np.float32)},
            {"sequence_hash": row.sequence_hash,
             "length": len(row.amino_acid_sequence),
             "chunk_count": chunks, "status": "ok"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--features", nargs="+", choices=["mean", "std", "parti"],
                     default=["mean", "std", "parti"])
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--shard-id", default="unknown")
    ap.add_argument("--max-tokens", type=int, default=8192,
                    help="maximum padded-batch token budget for <=1022-aa proteins")
    ap.add_argument("--max-batch-size", type=int, default=16)
    ap.add_argument("--disable-batching", action="store_true",
                    help="serial native-context reference mode for benchmarking")
    ap.add_argument("--pagerank-backend", choices=["networkx", "tensor"],
                    default="tensor")
    args = ap.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    checkpoint = out / "checkpoint.json"
    old = json.loads(checkpoint.read_text()) if checkpoint.exists() else {
        "completed_sequence_hashes": [], "failures": []
    }
    failures = list(old.get("failures", []))
    completed, _, _ = _load_existing(out, args.features)
    done = set(completed)

    df = pd.read_parquet(args.input)
    required = {"sequence_hash", "amino_acid_sequence"}
    if not required.issubset(df.columns):
        raise ValueError(f"input missing columns: {sorted(required - set(df.columns))}")
    df = df.drop_duplicates("sequence_hash").sort_values("sequence_hash").reset_index(drop=True)
    model, alphabet, device = load_model(args.model, args.device)
    cfg = ESMConfig()
    pending = [row for row in df.itertuples(index=False) if row.sequence_hash not in done]
    native = [row for row in pending if len(row.amino_acid_sequence) <= cfg.window_size]
    long_rows = [row for row in pending if len(row.amino_acid_sequence) > cfg.window_size]
    batches = ([[row] for row in native] if args.disable_batching else
               make_length_batches(native, args.max_tokens, args.max_batch_size))

    def flush():
        if completed:
            _write_outputs(out, completed, args.features)
        update_shard_checkpoint(checkpoint, sorted(completed), failures)

    def mark_success(sequence_hash, record):
        nonlocal failures
        completed[sequence_hash] = record
        failures = [failure for failure in failures
                    if failure.get("sequence_hash") != sequence_hash]

    for batch in batches:
        try:
            values = ([extract_one(model, alphabet, row.amino_acid_sequence, device, cfg,
                                   pagerank_backend=args.pagerank_backend)
                       for row in batch] if args.disable_batching else
                      extract_native_batch(model, alphabet,
                                           [row.amino_acid_sequence for row in batch],
                                           device, cfg,
                                           pagerank_backend=args.pagerank_backend))
            for row, result in zip(batch, values):
                mark_success(row.sequence_hash, _record(result, row, 1))
        except Exception:
            for row in batch:
                try:
                    result = extract_one(model, alphabet, row.amino_acid_sequence, device, cfg,
                                         pagerank_backend=args.pagerank_backend)
                    mark_success(row.sequence_hash, _record(result, row, 1))
                except Exception as exc:
                    failures.append({"sequence_hash": row.sequence_hash, "error": repr(exc)})
        flush()

    for row in long_rows:
        try:
            result = extract_one(model, alphabet, row.amino_acid_sequence, device, cfg,
                                 pagerank_backend=args.pagerank_backend)
            chunks = len(chunk_starts(len(row.amino_acid_sequence), cfg.window_size, cfg.overlap))
            mark_success(row.sequence_hash, _record(result, row, chunks))
        except Exception as exc:
            failures.append({"sequence_hash": row.sequence_hash, "error": repr(exc)})
        flush()

    print(json.dumps({
        "input_rows": len(df), "new_rows": len(pending),
        "completed": len(completed), "failures": len(failures),
        "native_batches": len(batches), "long_rows": len(long_rows),
        "max_tokens": args.max_tokens, "max_batch_size": args.max_batch_size,
        "batching": not args.disable_batching, "pagerank_backend": args.pagerank_backend,
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
