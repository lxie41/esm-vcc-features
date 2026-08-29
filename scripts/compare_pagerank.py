"""Compare NetworkX and tensor PageRank on a deterministic real subset."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from protein_embeddings.esm import (load_model,
                                     residue_embeddings_and_parti_attention_streaming)
from protein_embeddings.parti import (pagerank_weights, pagerank_weights_tensor,
                                       parti_attention_matrix)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--count", type=int, default=25)
    args = ap.parse_args()
    df = pd.read_parquet(args.input).drop_duplicates("sequence_hash").sort_values("sequence_hash")
    df = df.iloc[:min(args.count, len(df))]
    model, alphabet, device = load_model(args.model, args.device)
    rows = []
    for row in df.itertuples(index=False):
        states, attention = residue_embeddings_and_parti_attention_streaming(
            model, alphabet, row.amino_acid_sequence, device,
            torch.float16 if device.startswith("cuda") else None)
        matrix = parti_attention_matrix(attention.unsqueeze(0).unsqueeze(1))
        reference = pagerank_weights(matrix).numpy()
        optimized = pagerank_weights_tensor(matrix).numpy()
        midpoint = (reference + optimized) / 2
        js = 0.5 * np.sum(reference * np.log(reference / midpoint)) + 0.5 * np.sum(
            optimized * np.log(optimized / midpoint))
        vector_a = (states * torch.from_numpy(reference)[:, None]).sum(0).numpy()
        vector_b = (states * torch.from_numpy(optimized)[:, None]).sum(0).numpy()
        rows.append({
            "sequence_hash": row.sequence_hash,
            "length": len(row.amino_acid_sequence),
            "weight_cosine": float(np.dot(reference, optimized) /
                                    (np.linalg.norm(reference) * np.linalg.norm(optimized))),
            "weight_correlation": float(np.corrcoef(reference, optimized)[0, 1]),
            "max_weight_abs_error": float(np.max(np.abs(reference - optimized))),
            "js_distance": float(np.sqrt(max(js, 0.0))),
            "vector_cosine": float(np.dot(vector_a, vector_b) /
                                    (np.linalg.norm(vector_a) * np.linalg.norm(vector_b))),
            "max_vector_abs_error": float(np.max(np.abs(vector_a - vector_b))),
            "finite": bool(np.isfinite(optimized).all() and np.isfinite(vector_b).all()),
        })
    result = {"count": len(rows), "rows": rows,
              "mean_weight_cosine": float(np.mean([r["weight_cosine"] for r in rows])),
              "min_weight_cosine": min(r["weight_cosine"] for r in rows),
              "mean_weight_correlation": float(np.mean([r["weight_correlation"] for r in rows])),
              "max_weight_abs_error": max(r["max_weight_abs_error"] for r in rows),
              "max_js_distance": max(r["js_distance"] for r in rows),
              "mean_vector_cosine": float(np.mean([r["vector_cosine"] for r in rows])),
              "max_vector_abs_error": max(r["max_vector_abs_error"] for r in rows),
              "all_finite": all(r["finite"] for r in rows),
              "sequence_hash_order": [r["sequence_hash"] for r in rows]}
    Path(args.output).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
