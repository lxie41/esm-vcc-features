"""Compare bounded long-protein PaRTI candidates against direct PaRTI."""
from __future__ import annotations
import json, time
from pathlib import Path
import pandas as pd
import torch
from protein_embeddings.esm import load_model, residue_embeddings_and_parti_attention_streaming, reconstruct_residue_embeddings, reconstruction_weights, chunk_starts
from protein_embeddings.parti import pool_parti

def js(p, q):
    m=(p+q).clamp_min(1e-12)/2
    return float((0.5*(p*(p/m).log()).sum()+0.5*(q*(q/m).log()).sum()))

def main():
    out=Path('reports/esm_pilot'); df=pd.read_parquet('data/processed/mapping_final/esm_sequences_final.parquet'); bins=[(250,350),(450,550),(650,750),(850,950),(980,1023)]; rows=[]
    for lo,hi in bins: rows.extend([r for _,r in df[(df.protein_length>=lo)&(df.protein_length<hi)].sort_values('sequence_hash').head(4).iterrows()])
    model,alphabet,device=load_model('data/external/esm2/esm2_t33_650M_UR50D.pt','cuda'); results=[]
    for row in rows:
        seq=row.amino_acid_sequence; full_h,full_matrix=residue_embeddings_and_parti_attention_streaming(model,alphabet,seq,device); zfull,wfull=pool_parti(full_h,full_matrix.unsqueeze(0).unsqueeze(1))
        for overlap in [64,128,256]:
            starts=chunk_starts(len(seq),512,overlap); chunks=[]; locals=[]; t=time.perf_counter()
            for s in starts:
                h,matrix=residue_embeddings_and_parti_attention_streaming(model,alphabet,seq[s:s+512],device); z,w=pool_parti(h,matrix.unsqueeze(0).unsqueeze(1)); chunks.append(h); locals.append((z,w))
            H=reconstruct_residue_embeddings(chunks,starts,len(seq),overlap,'triangular'); global_w=torch.zeros(len(seq))
            for h,(s,(z,w)) in zip(chunks,zip(starts,locals)):
                taper=reconstruction_weights(len(h),overlap,'triangular'); global_w[s:s+len(h)] += w*taper
            global_w/=global_w.sum(); za=(H*global_w[:,None]).sum(0); chunk_weights=torch.tensor([len(h) for h in chunks],dtype=torch.float32); chunk_weights/=chunk_weights.sum(); zc=sum(c[0]*q for c,q in zip(locals,chunk_weights))
            results.extend([{'sequence_hash':row.sequence_hash,'length':len(seq),'candidate':'A_merge_local_importance','overlap':overlap,'chunks':len(starts),'seconds':time.perf_counter()-t,'cosine_to_full':float(torch.nn.functional.cosine_similarity(za,zfull,dim=0)),'l2_to_full':float((za-zfull).norm()),'js_weights':js(global_w,wfull),'weight_corr':float(torch.corrcoef(torch.stack((global_w,wfull)))[0,1]),'boundary_weight_mean':float(global_w[0]),'interior_weight_mean':float(global_w[len(seq)//4:3*len(seq)//4].mean()),'redundant_ratio':sum(len(h) for h in chunks)/len(seq)}, {'sequence_hash':row.sequence_hash,'length':len(seq),'candidate':'C_chunk_vector_control','overlap':overlap,'chunks':len(starts),'cosine_to_full':float(torch.nn.functional.cosine_similarity(zc,zfull,dim=0)),'l2_to_full':float((zc-zfull).norm()),'js_weights':None,'weight_corr':None,'redundant_ratio':sum(len(h) for h in chunks)/len(seq)}])
    (out/'parti_long_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8'); print(json.dumps({'proteins':len(rows),'rows':len(results),'output':str(out/'parti_long_results.json')},indent=2))
if __name__=='__main__': main()
