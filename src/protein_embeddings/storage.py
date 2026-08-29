"""Small, hash-keyed feature-bank prototype; no CSV vector columns."""

from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np


def write_feature_bank(root: str | Path, rows: list[tuple[str, np.ndarray]], name: str) -> None:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    hashes = np.array([h for h, _ in rows], dtype="U64")
    matrix = np.stack([v.astype(np.float32, copy=False) for _, v in rows])
    np.save(root / f"{name}.npy", matrix)
    np.save(root / f"{name}_hashes.npy", hashes)
    (root / f"{name}.json").write_text(json.dumps({"name": name, "rows": len(rows),
        "dim": int(matrix.shape[1]), "dtype": str(matrix.dtype)}), encoding="utf-8")


def update_shard_checkpoint(path: str | Path, completed_hashes: list[str],
                            failures: list[dict]) -> None:
    """Atomically persist resumable shard state; completed hashes are the unit of work."""
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"updated_utc": datetime.now(timezone.utc).isoformat(),
               "completed_sequence_hashes": sorted(set(completed_hashes)),
               "failures": failures}
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)
