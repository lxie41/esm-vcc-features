"""Resumable Colab/Drive orchestration for production shards.

This runner performs no embedding logic itself. It copies one manifest shard to
local scratch, invokes the shared extractor, validates the completed local
feature bank, then publishes that shard to the configured Drive output.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_output(output: Path, expected_rows: int) -> dict:
    checkpoint = output / "checkpoint.json"
    metadata_path = output / "metadata.parquet"
    if not checkpoint.exists() or not metadata_path.exists():
        raise RuntimeError(f"incomplete output: {output}")
    checkpoint_data = json.loads(checkpoint.read_text())
    metadata = pd.read_parquet(metadata_path)
    if len(metadata) != expected_rows:
        raise RuntimeError(f"metadata count {len(metadata)} != {expected_rows}")
    if metadata.sequence_hash.nunique() != expected_rows:
        raise RuntimeError("metadata sequence_hash values are not unique")
    if not metadata.sequence_hash.is_monotonic_increasing:
        raise RuntimeError("metadata sequence_hash ordering is not deterministic")
    if len(checkpoint_data.get("completed_sequence_hashes", [])) != expected_rows:
        raise RuntimeError("checkpoint completed count does not match metadata")
    if checkpoint_data.get("failures"):
        raise RuntimeError("output contains recorded failures")
    expected_dims = {"mean": 1280, "std": 1280, "parti": 1280}
    for feature, dimension in expected_dims.items():
        matrix_path = output / f"{feature}.npy"
        if not matrix_path.exists():
            raise RuntimeError(f"missing feature matrix: {matrix_path}")
        matrix = np.load(matrix_path, mmap_mode="r")
        if matrix.shape != (expected_rows, dimension):
            raise RuntimeError(f"{feature} shape {matrix.shape} is invalid")
        if not np.isfinite(matrix).all():
            raise RuntimeError(f"{feature} contains non-finite values")
    return {"rows": expected_rows, "features": list(expected_dims), "failures": 0}


def validate_input(input_path: Path, record: dict):
    frame = pd.read_parquet(input_path, columns=["sequence_hash", "amino_acid_sequence"])
    if len(frame) != record["row_count"]:
        raise RuntimeError(f"input row count {len(frame)} != {record['row_count']}")
    if frame.sequence_hash.nunique() != len(frame):
        raise RuntimeError("input sequence_hash values are not unique")
    ordered = frame.sort_values("sequence_hash")
    if ordered.sequence_hash.iloc[0] != record["first_sequence_hash"]:
        raise RuntimeError("input first hash does not match manifest")
    if ordered.sequence_hash.iloc[-1] != record["last_sequence_hash"]:
        raise RuntimeError("input last hash does not match manifest")
    if ordered.amino_acid_sequence.isna().any():
        raise RuntimeError("input contains missing sequences")


def publish(local_output: Path, drive_output: Path):
    drive_output.parent.mkdir(parents=True, exist_ok=True)
    temporary = drive_output.with_name(drive_output.name + ".uploading")
    if temporary.exists():
        shutil.rmtree(temporary)
    shutil.copytree(local_output, temporary)
    if drive_output.exists():
        shutil.rmtree(drive_output)
    temporary.replace(drive_output)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--input-root", required=True,
                    help="directory containing manifest-listed input shards")
    ap.add_argument("--drive-output-root", required=True)
    ap.add_argument("--local-root", required=True,
                    help="local scratch directory, e.g. /content/esm_work")
    ap.add_argument("--model", required=True)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--max-tokens", type=int, default=8192)
    ap.add_argument("--max-batch-size", type=int, default=16)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--extractor", default="scripts/extract_esm_features.py")
    ap.add_argument("--limit", type=int, default=None,
                    help="optional dry-run limit; never required for production")
    args = ap.parse_args()

    manifest_path = Path(args.manifest)
    manifest = json.loads(manifest_path.read_text())
    records = manifest["shards"][:args.limit]
    input_root, drive_root, local_root = map(Path, (args.input_root, args.drive_output_root, args.local_root))
    local_root.mkdir(parents=True, exist_ok=True)

    for record in records:
        shard_id = record["shard_id"]
        input_path = input_root / record["filename"]
        drive_output = drive_root / f"shard_{shard_id}"
        if drive_output.exists():
            try:
                result = validate_output(drive_output, record["row_count"])
                print(json.dumps({"shard_id": shard_id, "status": "already_complete", **result}))
                continue
            except Exception as exc:
                print(json.dumps({"shard_id": shard_id, "status": "retry_incomplete_drive_output", "reason": repr(exc)}))

        local_input = local_root / record["filename"]
        local_output = local_root / f"shard_{shard_id}"
        if local_input.exists():
            local_input.unlink()
        if local_output.exists():
            shutil.rmtree(local_output)
        shutil.copy2(input_path, local_input)
        if sha256(local_input) != record["sha256"]:
            raise RuntimeError(f"input checksum mismatch for shard {shard_id}")
        validate_input(local_input, record)

        command = [args.python, args.extractor, "--input", str(local_input),
                   "--output", str(local_output), "--model", args.model,
                   "--features", "mean", "std", "parti", "--device", args.device,
                   "--shard-id", shard_id, "--max-tokens", str(args.max_tokens),
                   "--max-batch-size", str(args.max_batch_size)]
        try:
            subprocess.run(command, check=True)
            result = validate_output(local_output, record["row_count"])
            publish(local_output, drive_output)
            print(json.dumps({"shard_id": shard_id, "status": "published", **result}))
        except Exception as exc:
            print(json.dumps({"shard_id": shard_id, "status": "failed", "reason": repr(exc)}))
            raise
        finally:
            if local_input.exists():
                local_input.unlink()
            if local_output.exists():
                shutil.rmtree(local_output)


if __name__ == "__main__":
    main()
