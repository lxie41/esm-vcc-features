"""Faithful Pool PaRTI implementation for ESM2 attention outputs."""

from __future__ import annotations

import numpy as np
import networkx as nx
import torch
import time


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
                     tol: float = 1e-6, max_iter: int = 100, timings=None) -> torch.Tensor:
    """Run the official NetworkX weighted PageRank and remove BOS/EOS."""
    matrix = attention.detach().cpu().numpy().astype(np.float64, copy=False)
    start = time.perf_counter()
    graph = nx.from_numpy_array(matrix, create_using=nx.DiGraph)
    if timings is not None:
        timings["networkx_graph_seconds"] = timings.get("networkx_graph_seconds", 0.0) + time.perf_counter() - start
    start = time.perf_counter()
    scores = nx.pagerank(graph, alpha=alpha, tol=tol, weight="weight", max_iter=max_iter)
    if timings is not None:
        timings["networkx_pagerank_seconds"] = timings.get("networkx_pagerank_seconds", 0.0) + time.perf_counter() - start
    values = np.array([scores[i] for i in range(matrix.shape[0] - 1) if i != 0], dtype=np.float64)
    values /= values.sum()
    return torch.from_numpy(values.astype(np.float32))


def pagerank_weights_tensor(attention: torch.Tensor, alpha: float = 0.85,
                            tol: float = 1e-6, max_iter: int = 100,
                            device: str | None = None, timings=None) -> torch.Tensor:
    """Tensor PageRank matching NetworkX's weighted directed implementation.

    The update is the same power iteration used by NetworkX: row-normalize
    outgoing edge weights, send dangling mass to the uniform personalization
    vector, apply alpha=0.85 teleportation, and stop at L1 error < N*tol.
    Float64 is intentional for agreement with NetworkX's NumPy float64 path;
    only the final residue weights are returned as float32.
    """
    start = time.perf_counter()
    matrix = attention.detach().to(device=device or attention.device, dtype=torch.float64)
    n = matrix.shape[0]
    if matrix.ndim != 2 or matrix.shape[1] != n:
        raise ValueError("attention must be a square matrix")
    personalization = torch.full((n,), 1.0 / n, dtype=torch.float64, device=matrix.device)
    row_sum = matrix.sum(dim=1)
    non_dangling = row_sum != 0
    transition = torch.zeros_like(matrix)
    transition[non_dangling] = matrix[non_dangling] / row_sum[non_dangling, None]
    x = personalization.clone()
    for _ in range(max_iter):
        last = x
        dangling_mass = last[~non_dangling].sum()
        x = alpha * (last @ transition)
        x = x + alpha * dangling_mass * personalization
        x = x + (1.0 - alpha) * personalization
        if torch.abs(x - last).sum().item() < n * tol:
            break
    values = x[1:-1]
    values = values / values.sum()
    if timings is not None:
        timings["tensor_pagerank_seconds"] = timings.get("tensor_pagerank_seconds", 0.0) + time.perf_counter() - start
    result = values.detach().to("cpu", dtype=torch.float32)
    if not torch.isfinite(result).all():
        raise FloatingPointError("non-finite tensor PageRank weights")
    return result


def pool_parti(residue_states: torch.Tensor, attentions: torch.Tensor,
               alpha: float = 0.85, timings=None,
               pagerank_backend: str = "networkx") -> tuple[torch.Tensor, torch.Tensor]:
    """Return 1280-d weighted pooling and normalized residue importance."""
    if residue_states.ndim != 2 or attentions.ndim != 4:
        raise ValueError("invalid residue or attention dimensions")
    start = time.perf_counter()
    matrix = parti_attention_matrix(attentions)
    if pagerank_backend == "networkx":
        weights = pagerank_weights(matrix, alpha=alpha, timings=timings)
    elif pagerank_backend == "tensor":
        weights = pagerank_weights_tensor(matrix, alpha=alpha, timings=timings)
    else:
        raise ValueError(f"unknown PageRank backend: {pagerank_backend}")
    if timings is not None:
        timings["pagerank_total_seconds"] = timings.get("pagerank_total_seconds", 0.0) + time.perf_counter() - start
    if len(weights) != residue_states.shape[0]:
        raise ValueError("BOS/EOS-stripped attention and residue lengths differ")
    start = time.perf_counter()
    vector = (residue_states.float() * weights[:, None]).sum(dim=0)
    if timings is not None:
        timings["parti_weighted_pooling_seconds"] = timings.get("parti_weighted_pooling_seconds", 0.0) + time.perf_counter() - start
    return vector, weights
