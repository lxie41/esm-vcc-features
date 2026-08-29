import torch
from protein_embeddings.parti import parti_attention_matrix, pagerank_weights, pool_parti

def test_parti_excludes_special_tokens_and_normalizes():
    att=torch.rand(33,20,5,5); h=torch.rand(3,8); z,w=pool_parti(h,att)
    assert z.shape==(8,); assert w.shape==(3,); assert torch.isclose(w.sum(),torch.tensor(1.0)); assert torch.isfinite(z).all()

def test_parti_attention_reduction_and_determinism():
    att=torch.arange(33*20*4*4,dtype=torch.float32).reshape(33,20,4,4); x=parti_attention_matrix(att)
    assert x.shape==(4,4); assert torch.equal(x,parti_attention_matrix(att)); assert torch.equal(pagerank_weights(x),pagerank_weights(x))
import numpy as np
import torch

from protein_embeddings.parti import pagerank_weights, pagerank_weights_tensor


def test_tensor_pagerank_matches_networkx_reference():
    rng = np.random.default_rng(17)
    matrix = rng.random((12, 12), dtype=np.float64)
    matrix[3] = 0.0  # explicit dangling-node behavior
    reference = pagerank_weights(torch.from_numpy(matrix))
    optimized = pagerank_weights_tensor(torch.from_numpy(matrix))
    assert torch.isfinite(optimized).all()
    assert torch.allclose(reference, optimized, atol=2e-6, rtol=2e-5)
    assert torch.allclose(optimized.sum(), torch.tensor(1.0), atol=1e-6)
