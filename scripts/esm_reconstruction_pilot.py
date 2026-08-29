"""Bounded direct-vs-overlap reconstruction experiment; never touches mapping."""
from __future__ import annotations
import json, time
from pathlib import Path
import pandas as pd
import torch
from protein_embeddings.esm import ESMConfig, load_model, residue_embeddings, reconstruct_residue_embeddings, chunk_starts
from protein_embeddings.pooling import mean_representation, mean_sd_representation

def metrics(ref: torch.Tensor, got: torch.Tensor) -> dict:
    cos = torch.nn.functional.cosine_similarity(ref, got, dim=1).clamp(-1, 1)
    err = (ref - got).float()
    l2 = err.norm(dim=1)
    return {"mean_cosine": float(cos.mean()), "median_cosine": float(cos.median()),
            "p01_cosine": float(torch.quantile(cos, .01)), "mean_rmse": float(err.square().mean().sqrt()),
            "mean_l2": float(l2.mean()), "p99_l2": float(torch.quantile(l2, .99))}

def main() -> None:
    out = Path("reports/esm_pilot"); out.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet("data/processed/mapping_final/esm_sequences_final.parquet")
    edges = [0, 300, 500, 700, 900, 1022]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        hit = df[(df.protein_length >= lo) & (df.protein_length < hi)].sort_values("sequence_hash").head(6)
        rows.extend([r for _, r in hit.iterrows()])
    model, alphabet, device = load_model("data/external/esm2/esm2_t33_650M_UR50D.pt", "cuda")
    results = []
    for row in rows:
        seq = row.amino_acid_sequence; ref = residue_embeddings(model, alphabet, seq, 33, device)
        # Forced 512-residue stress window is necessary because all references
        # fit native 1022-residue context; production remains 1022 residues.
        for window in [512]:
            for overlap in [64, 128, 256]:
                starts = chunk_starts(len(seq), window, overlap)
                chunks = [residue_embeddings(model, alphabet, seq[s:s+window], 33, device) for s in starts]
                for weight in ["uniform", "triangular", "cosine"]:
                    t = time.perf_counter(); got = reconstruct_residue_embeddings(chunks, starts, len(seq), overlap, weight)
                    m = metrics(ref, got); m.update({"sequence_hash": row.sequence_hash, "length": len(seq),
                        "window": window, "overlap": overlap, "weighting": weight,
                        "chunk_count": len(starts), "redundant_residue_ratio": sum(map(len, chunks))/len(seq),
                        "seconds_reconstruct": time.perf_counter()-t})
                    mean_delta = (mean_representation(ref)-mean_representation(got)).norm()
                    msd_delta = (mean_sd_representation(ref)-mean_sd_representation(got)).norm()
                    m.update({"mean_vector_l2": float(mean_delta), "mean_sd_vector_l2": float(msd_delta)})
                    results.append(m)
        # Negative controls for the same short reference.
        trunc = ref[:min(1022, len(ref))]
        results.append({"sequence_hash": row.sequence_hash, "length": len(seq), "control": "hard_truncate",
                        **metrics(ref, torch.cat((trunc, torch.zeros((len(ref)-len(trunc), ref.shape[1])))) if len(trunc)<len(ref) else trunc)})
    (out / "reconstruction_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({"rows": len(results), "proteins": len(rows), "output": str(out / 'reconstruction_results.json')}, indent=2))

if __name__ == "__main__": main()
