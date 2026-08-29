from types import SimpleNamespace
import json

import numpy as np
import pandas as pd

from scripts.extract_esm_features import make_length_batches
from scripts.run_production_shards import validate_output


def row(sequence_hash, length):
    return SimpleNamespace(sequence_hash=sequence_hash, amino_acid_sequence="A" * length)


def test_length_batches_respect_token_and_size_limits():
    rows = [row("c", 100), row("a", 50), row("b", 75)]
    batches = make_length_batches(rows, max_tokens=160, max_batch_size=2)
    assert [[r.sequence_hash for r in batch] for batch in batches] == [["a", "b"], ["c"]]
    assert all((max(len(r.amino_acid_sequence) for r in batch) + 2) * len(batch) <= 160
               for batch in batches)
    assert all(len(batch) <= 2 for batch in batches)


def test_length_batches_are_deterministic_for_equal_lengths():
    rows = [row("z", 100), row("x", 100), row("y", 100)]
    first = make_length_batches(rows, max_tokens=1000, max_batch_size=2)
    second = make_length_batches(rows, max_tokens=1000, max_batch_size=2)
    assert [[r.sequence_hash for r in b] for b in first] == [["x", "y"], ["z"]]
    assert [[r.sequence_hash for r in b] for b in first] == [[r.sequence_hash for r in b] for b in second]


def test_runner_validates_feature_bank_schema(tmp_path):
    output = tmp_path / "shard_0001"
    output.mkdir()
    hashes = ["a", "b"]
    pd.DataFrame({"sequence_hash": hashes, "length": [3, 4],
                  "chunk_count": [1, 1], "status": ["ok", "ok"]}).to_parquet(
                      output / "metadata.parquet", index=False)
    for feature in ("mean", "std", "parti"):
        np.save(output / f"{feature}.npy", np.zeros((2, 1280), dtype=np.float32))
    (output / "checkpoint.json").write_text(json.dumps({
        "completed_sequence_hashes": hashes, "failures": []
    }))
    assert validate_output(output, 2)["rows"] == 2
