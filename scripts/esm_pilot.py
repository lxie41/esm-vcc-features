"""Run bounded ESM2 setup/token/pooling/reconstruction/storage pilot.

This script never iterates over the complete frozen sequence universe.
"""
from __future__ import annotations

import argparse, hashlib, json, platform, time
from pathlib import Path
import pandas as pd
import torch

from protein_embeddings.esm import ESMConfig, load_model, extract_sequence, chunk_starts
from protein_embeddings.pooling import mean_representation, mean_sd_representation
from protein_embeddings.storage import write_feature_bank


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="data/external/esm2")
    ap.add_argument("--out", default="reports/esm_pilot")
    ap.add_argument("--n", type=int, default=8)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    seqs = pd.read_parquet("data/processed/mapping_final/esm_sequences_final.parquet")
    # Deterministic bounded selection across lengths; hash is the identity.
    bins = [0, 300, 500, 700, 900, 1000, 1022]
    selected = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        hit = seqs[(seqs.protein_length >= lo) & (seqs.protein_length < hi)].sort_values("sequence_hash").head(1)
        if len(hit): selected.append(hit.iloc[0])
    selected += [seqs.sort_values("sequence_hash").iloc[i] for i in range(min(args.n-len(selected), len(seqs)))]
    selected = selected[:args.n]
    cfg = ESMConfig()
    checkpoint = Path(args.model_dir) / "esm2_t33_650M_UR50D.pt"
    model, alphabet, device = load_model(str(checkpoint), "cpu")
    report = {"python": platform.python_version(), "torch": torch.__version__,
              "torch_cuda": torch.version.cuda, "cuda_available": torch.cuda.is_available(),
              "device": device, "model": cfg.__dict__, "selected": []}
    converter = alphabet.get_batch_converter()
    token_checks = []
    for s in ["ACDXU", "X", "U", "ACD"]:
        _, _, t = converter([("synthetic", s)])
        token_checks.append({"sequence": s, "tokens": t[0].tolist(),
                             "finite_input": bool(torch.isfinite(t.float()).all())})
    report["token_checks"] = token_checks
    feature_rows = {"mean": [], "mean_sd": []}
    for row in selected:
        seq = row.amino_acid_sequence; t0 = time.perf_counter()
        h, meta = extract_sequence(model, alphabet, seq, cfg, device)
        h2, _ = extract_sequence(model, alphabet, seq, cfg, device)
        mean = mean_representation(h).numpy(); msd = mean_sd_representation(h).numpy()
        report["selected"].append({"sequence_hash": row.sequence_hash, "length": len(seq),
            "sha256_check": hashlib.sha256(seq.encode()).hexdigest(), "residue_shape": list(h.shape),
            "chunk_count": meta["chunk_count"], "seconds_two_pass": time.perf_counter()-t0,
            "max_abs_repeat": float((h-h2).abs().max()), "mean_abs_repeat": float((h-h2).abs().mean()),
            "cosine_repeat": float(torch.nn.functional.cosine_similarity(h.flatten(), h2.flatten(), dim=0)),
            "finite": bool(torch.isfinite(h).all()), "mean_dim": len(mean), "mean_sd_dim": len(msd)})
        feature_rows["mean"].append((row.sequence_hash, mean)); feature_rows["mean_sd"].append((row.sequence_hash, msd))
    for name, rows in feature_rows.items(): write_feature_bank(out / "feature_bank", rows, name)
    (out / "pilot_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

if __name__ == "__main__": main()
