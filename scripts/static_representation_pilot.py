"""Compare deterministic static summaries on a bounded GPU pilot."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import torch
from protein_embeddings.esm import load_model, residue_embeddings
from protein_embeddings.pooling import mean_representation, mean_sd_representation

def quantile_representation(h: torch.Tensor) -> torch.Tensor:
    x = h.float(); q = torch.quantile(x, torch.tensor([.25, .50, .75]), dim=0)
    return torch.cat((x.mean(0), x.std(0, correction=0), q.flatten()), 0)

def main() -> None:
    out = Path("reports/esm_pilot"); out.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet("data/processed/mapping_final/esm_sequences_final.parquet")
    edges = [0, 300, 500, 700, 900, 1022]; rows=[]
    for lo, hi in zip(edges[:-1], edges[1:]):
        rows.extend([r for _, r in df[(df.protein_length>=lo)&(df.protein_length<hi)].sort_values('sequence_hash').head(6).iterrows()])
    model, alphabet, device = load_model("data/external/esm2/esm2_t33_650M_UR50D.pt", "cuda")
    feats={"mean":[],"mean_sd":[],"mean_sd_q25_q50_q75":[]}; lengths=[]
    for row in rows:
        h=residue_embeddings(model,alphabet,row.amino_acid_sequence,33,device); lengths.append(len(h))
        feats['mean'].append(mean_representation(h)); feats['mean_sd'].append(mean_sd_representation(h)); feats['mean_sd_q25_q50_q75'].append(quantile_representation(h))
    report={"proteins":len(rows),"lengths":lengths,"families":{}}
    for name, values in feats.items():
        X=torch.stack(values); norms=X.norm(dim=1); centered=X-X.mean(0); s=torch.linalg.svdvals(centered);
        pair=torch.nn.functional.cosine_similarity(X[:,None,:],X[None,:,:],dim=2); off=pair[~torch.eye(len(X),dtype=bool)]
        report['families'][name]={"dimension":X.shape[1],"mean_norm":float(norms.mean()),"norm_length_corr":float(torch.corrcoef(torch.stack((norms,torch.tensor(lengths,dtype=torch.float))))[0,1]),"median_pair_cosine":float(off.median()),"effective_rank":float((s.square().sum()**2/s.pow(4).sum()).item()),"top5_pca_fraction":float((s[:5].square().sum()/s.square().sum()).item()),"finite":bool(torch.isfinite(X).all())}
    (out/'static_representation_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8'); print(json.dumps(report,indent=2))
if __name__ == '__main__': main()
