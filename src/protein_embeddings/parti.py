"""Faithful Pool PaRTI implementation for ESM2 attention outputs."""

from __future__ import annotations

import numpy as np
import networkx as nx
import torch


def parti_attention_matrix(attentions: torch.Tensor) -> torch.Tensor:
    """Match the official implementation: max-head, then max-layer attention."""
    if attentions.ndim != 4:
        raise ValueError("expected layers x heads x tokens x tokens attention")
    if attentions.shape[-1] != attentions.shape[-2]:
        raise ValueError("attention matrices must be square")
    # Official extraction first computes mean/max per layer and then Pool PaRTI
    # selects the max-head branch and max-pools across all 33 layers.
    return attentions.max(dim=1).values.max(dim=0).values


def pagerank_weights(attention: torch.Tensor, alpha: float = 0.85,
                     tol: float = 1e-6, max_iter: int = 100) -> torch.Tensor:
    """Run the official NetworkX weighted PageRank and remove BOS/EOS."""
    matrix = attention.detach().cpu().numpy().astype(np.float64, copy=False)
    graph = nx.from_numpy_array(matrix, create_using=nx.DiGraph)
    scores = nx.pagerank(graph, alpha=alpha, tol=tol, weight="weight", max_iter=max_iter)
    values = np.array([scores[i] for i in range(matrix.shape[0] - 1) if i != 0], dtype=np.float64)
    values /= values.sum()
    return torch.from_numpy(values.astype(np.float32))


def pool_parti(residue_states: torch.Tensor, attentions: torch.Tensor,
               alpha: float = 0.85) -> tuple[torch.Tensor, torch.Tensor]:
    """Return 1280-d weighted pooling and normalized residue importance."""
    if residue_states.ndim != 2 or attentions.ndim != 4:
        raise ValueError("invalid residue or attention dimensions")
    weights = pagerank_weights(parti_attention_matrix(attentions), alpha=alpha)
    if len(weights) != residue_states.shape[0]:
        raise ValueError("BOS/EOS-stripped attention and residue lengths differ")
    vector = (residue_states.float() * weights[:, None]).sum(dim=0)
    return vector, weights
