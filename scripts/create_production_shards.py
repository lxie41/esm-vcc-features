"""Create deterministic, immutable-input Parquet shards for production extraction."""
from __future__ import annotations
import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

REQUIRED = ["sequence_hash", "amino_acid_sequence", "protein_length"]

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data/processed/mapping_final/esm_sequences_final.parquet")
    ap.add_argument("--output-dir", default="data/production_shards")
    ap.add_argument("--shard-size", type=int, default=2000)
    args = ap.parse_args()
    source, output = Path(args.input), Path(args.output_dir)
    if args.shard_size <= 0: raise ValueError("--shard-size must be positive")
    df = pd.read_parquet(source)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing: raise ValueError(f"missing required extractor columns: {missing}")
    if len(df) != 174585: raise ValueError(f"frozen source row count changed: {len(df)}")
    if df.sequence_hash.isna().any() or (df.sequence_hash.astype(str).str.len() == 0).any(): raise ValueError("missing sequence_hash values")
    if df.amino_acid_sequence.isna().any() or (df.amino_acid_sequence.astype(str).str.len() == 0).any(): raise ValueError("missing amino_acid_sequence values")
    duplicate_count = int(df.sequence_hash.duplicated().sum())
    if duplicate_count: raise ValueError(f"duplicate sequence_hash values: {duplicate_count}")
    df = df.sort_values("sequence_hash", kind="mergesort").reset_index(drop=True); output.mkdir(parents=True, exist_ok=True)
    for old in output.glob("sequence_shard_*.parquet"): old.unlink()
    records = []
    for shard_num, start in enumerate(range(0, len(df), args.shard_size)):
        part = df.iloc[start:start + args.shard_size]; name = f"sequence_shard_{shard_num:04d}.parquet"; path = output / name
        part.to_parquet(path, index=False)
        records.append({"shard_id": f"{shard_num:04d}", "filename": name, "row_count": len(part), "first_sequence_hash": str(part.sequence_hash.iloc[0]), "last_sequence_hash": str(part.sequence_hash.iloc[-1]), "sha256": sha256_file(path)})
    manifest = {"source": str(source), "source_sha256": sha256_file(source), "created_utc": datetime.now(timezone.utc).isoformat(), "shard_size": args.shard_size, "total_rows": len(df), "unique_sequence_hashes": int(df.sequence_hash.nunique()), "duplicates": duplicate_count, "missing_hashes": int(df.sequence_hash.isna().sum()), "missing_sequences": int(df.amino_acid_sequence.isna().sum()), "columns": list(df.columns), "shards": records}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (output / "validation.json").write_text(json.dumps({"total_rows": len(df), "unique_sequence_hashes": int(df.sequence_hash.nunique()), "duplicates": duplicate_count, "missing_hashes": int(df.sequence_hash.isna().sum()), "missing_sequences": int(df.amino_acid_sequence.isna().sum()), "sequence_hash_sorted": bool(df.sequence_hash.is_monotonic_increasing), "shard_count": len(records)}, indent=2), encoding="utf-8")
    summary = {k: manifest[k] for k in manifest if k != "shards"}
    summary["shard_count"] = len(records)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__": main()
