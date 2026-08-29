import torch
from protein_embeddings.parti import parti_attention_matrix, pagerank_weights, pool_parti

def test_parti_excludes_special_tokens_and_normalizes():
    att=torch.rand(33,20,5,5); h=torch.rand(3,8); z,w=pool_parti(h,att)
    assert z.shape==(8,); assert w.shape==(3,); assert torch.isclose(w.sum(),torch.tensor(1.0)); assert torch.isfinite(z).all()

def test_parti_attention_reduction_and_determinism():
    att=torch.arange(33*20*4*4,dtype=torch.float32).reshape(33,20,4,4); x=parti_attention_matrix(att)
    assert x.shape==(4,4); assert torch.equal(x,parti_attention_matrix(att)); assert torch.equal(pagerank_weights(x),pagerank_weights(x))
