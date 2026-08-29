"""GPU smoke and precision pilot; bounded lengths only."""
from __future__ import annotations
import json, time
from pathlib import Path
import torch
from protein_embeddings.esm import load_model, residue_embeddings
from protein_embeddings.pooling import mean_representation, mean_sd_representation

def main() -> None:
    out = Path("reports/esm_pilot"); out.mkdir(parents=True, exist_ok=True)
    seq = "ACDEFGHIKLMNPQRSTVWY" * 52
    lengths = [100, 300, 500, 800, 1022]
    results = []
    for precision in ["fp32", "fp16"]:
        model, alphabet, device = load_model("data/external/esm2/esm2_t33_650M_UR50D.pt", "cuda")
        for length in lengths:
            s = seq[:length]; torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
            start = time.perf_counter(); h = residue_embeddings(model, alphabet, s, 33, device,
                torch.float16 if precision == "fp16" else None)
            elapsed = time.perf_counter() - start
            results.append({"precision": precision, "length": length, "tokens": length + 2,
                "shape": list(h.shape), "seconds": elapsed, "residues_per_second": length / elapsed,
                "peak_vram_bytes": torch.cuda.max_memory_allocated(), "finite": bool(torch.isfinite(h).all())})
        del model
    # Compare pooled FP16 to FP32 on one representative sequence using CPU tensors from calls above.
    m32, a32, _ = load_model("data/external/esm2/esm2_t33_650M_UR50D.pt", "cuda")
    h32 = residue_embeddings(m32, a32, seq[:500], 33, "cuda")
    h16 = residue_embeddings(m32, a32, seq[:500], 33, "cuda", torch.float16)
    for name, fn in [("mean", mean_representation), ("mean_sd", mean_sd_representation)]:
        x, y = fn(h32), fn(h16); diff = (x-y).abs();
        results.append({"precision_comparison": name, "length": 500,
            "cosine": float(torch.nn.functional.cosine_similarity(x, y, dim=0)),
            "max_abs": float(diff.max()), "mean_abs": float(diff.mean()),
            "relative_l2": float(diff.norm() / x.norm())})
    (out / "cuda_precision_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({"device": torch.cuda.get_device_name(0), "torch": torch.__version__,
        "cuda": torch.version.cuda, "results": len(results)}, indent=2))

if __name__ == "__main__": main()
