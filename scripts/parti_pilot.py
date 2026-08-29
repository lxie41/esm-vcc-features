"""Bounded Pool PaRTI validation on direct-context proteins."""
from __future__ import annotations
import json, time
from pathlib import Path
import pandas as pd
import torch
from protein_embeddings.esm import load_model, residue_embeddings_and_attention
from protein_embeddings.parti import pool_parti
from protein_embeddings.pooling import mean_representation

def main() -> None:
    out=Path('reports/esm_pilot'); out.mkdir(parents=True,exist_ok=True)
    df=pd.read_parquet('data/processed/mapping_final/esm_sequences_final.parquet')
    bins=[(250,350),(450,550),(650,750),(850,950),(980,1023)]; rows=[]
    for lo,hi in bins: rows.extend([r for _,r in df[(df.protein_length>=lo)&(df.protein_length<hi)].sort_values('sequence_hash').head(5).iterrows()])
    model,alphabet,device=load_model('data/external/esm2/esm2_t33_650M_UR50D.pt','cuda'); results=[]
    for row in rows:
        torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats(); t=time.perf_counter()
        h,att=residue_embeddings_and_attention(model,alphabet,row.amino_acid_sequence,33,device)
        z,w=pool_parti(h,att); elapsed=time.perf_counter()-t; pos=torch.arange(len(w),dtype=torch.float32)
        entropy=-(w*w.clamp_min(1e-12).log()).sum(); mean=mean_representation(h)
        results.append({'sequence_hash':row.sequence_hash,'length':len(w),'attention_shape':list(att.shape),'seconds':elapsed,'residues_per_second':len(w)/elapsed,'peak_vram_bytes':torch.cuda.max_memory_allocated(),'alpha_sum':float(w.sum()),'alpha_min':float(w.min()),'alpha_max':float(w.max()),'entropy':float(entropy),'effective_residues':float(entropy.exp()),'position_corr':float(torch.corrcoef(torch.stack((w,pos)))[0,1]),'parti_norm':float(z.norm()),'mean_norm':float(mean.norm()),'mean_parti_cosine':float(torch.nn.functional.cosine_similarity(z,mean,dim=0)),'finite':bool(torch.isfinite(z).all() and torch.isfinite(w).all())})
    (out/'parti_pilot_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8'); print(json.dumps({'proteins':len(results),'output':str(out/'parti_pilot_results.json')},indent=2))
if __name__=='__main__': main()
