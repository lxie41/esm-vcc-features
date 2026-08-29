"""Resumable hash-keyed production extractor; invocation only, never auto-run."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from protein_embeddings.esm import (ESMConfig, load_model,
    residue_embeddings_and_parti_attention_streaming, reconstruct_residue_embeddings,
    reconstruction_weights, chunk_starts)
from protein_embeddings.pooling import mean_representation, mean_sd_representation
from protein_embeddings.parti import pool_parti
from protein_embeddings.storage import update_shard_checkpoint

def extract_one(model, alphabet, seq, device, cfg):
    starts=chunk_starts(len(seq),cfg.window_size,cfg.overlap); chunks=[]; local=[]
    for s in starts:
        h,mat=residue_embeddings_and_parti_attention_streaming(model,alphabet,seq[s:s+cfg.window_size],device,torch.float16 if device.startswith('cuda') else None)
        chunks.append(h)
        z,w=pool_parti(h,mat.unsqueeze(0).unsqueeze(1)); local.append((z,w))
    H=reconstruct_residue_embeddings(chunks,starts,len(seq),cfg.overlap,cfg.weighting)
    if len(starts)==1: parti=local[0][0]; weights=local[0][1]
    else:
        weights=torch.zeros(len(seq))
        for h,(s,(_,w)) in zip(chunks,zip(starts,local)):
            weights[s:s+len(h)] += w*reconstruction_weights(len(h),cfg.overlap,cfg.weighting)
        weights/=weights.sum(); parti=(H*weights[:,None]).sum(0)
    return mean_representation(H).numpy(), mean_sd_representation(H).numpy(), parti.numpy(), len(starts)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--output',required=True); ap.add_argument('--model',required=True); ap.add_argument('--features',nargs='+',choices=['mean','std','parti'],default=['mean','std','parti']); ap.add_argument('--device',default='cuda'); ap.add_argument('--shard-id',default='unknown'); args=ap.parse_args()
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True); ck=out/'checkpoint.json'
    old=json.loads(ck.read_text()) if ck.exists() else {'completed_sequence_hashes':[],'failures':[]}; done=set(old['completed_sequence_hashes']); failures=list(old['failures'])
    df=pd.read_parquet(args.input); df=df.drop_duplicates('sequence_hash').sort_values('sequence_hash'); model,alphabet,device=load_model(args.model,args.device); cfg=ESMConfig(); rows=[]; matrices={f:[] for f in args.features}
    for row in df.itertuples(index=False):
        if row.sequence_hash in done: continue
        try:
            mean,msd,parti,chunks=extract_one(model,alphabet,row.amino_acid_sequence,device,cfg); vals={'mean':mean,'std':msd[1280:],'parti':parti}; rows.append({'sequence_hash':row.sequence_hash,'length':len(row.amino_acid_sequence),'chunk_count':chunks,'status':'ok'}); [matrices[f].append(vals[f]) for f in args.features]; done.add(row.sequence_hash)
        except Exception as exc:
            failures.append({'sequence_hash':row.sequence_hash,'error':repr(exc)}); rows.append({'sequence_hash':row.sequence_hash,'length':len(row.amino_acid_sequence),'status':'failed','error':repr(exc)})
        update_shard_checkpoint(ck,sorted(done),failures)
    if rows:
        meta=pd.DataFrame(rows)
        if (out/'metadata.parquet').exists(): meta=pd.concat([pd.read_parquet(out/'metadata.parquet'),meta],ignore_index=True).drop_duplicates('sequence_hash',keep='last')
        tmp=out/'metadata.parquet.tmp'; meta.to_parquet(tmp,index=False); tmp.replace(out/'metadata.parquet')
        for f,values in matrices.items():
            if values:
                current=np.load(out/f'{f}.npy') if (out/f'{f}.npy').exists() else None
                matrix=np.stack(values).astype('float32') if current is None else np.concatenate([current,np.stack(values).astype('float32')])
                np.save(out/f'{f}.npy.tmp',matrix); Path(out/f'{f}.npy.tmp.npy').replace(out/f'{f}.npy')
    print(json.dumps({'input_rows':len(df),'new_rows':len(rows),'completed':len(done),'failures':len(failures),'output':str(out)},indent=2))
if __name__=='__main__': main()
