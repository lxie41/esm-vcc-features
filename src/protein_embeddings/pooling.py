"""Deterministic residue-to-protein pooling used by the ESM pilot."""

from __future__ import annotations

import torch


def mean_representation(residue_states: torch.Tensor) -> torch.Tensor:
    return residue_states.float().mean(dim=0)


def mean_sd_representation(residue_states: torch.Tensor) -> torch.Tensor:
    x = residue_states.float()
    return torch.cat((x.mean(dim=0), x.std(dim=0, correction=0)), dim=0)


def simple_swe_representation(residue_states: torch.Tensor, reference: torch.Tensor,
                              projections: torch.Tensor) -> torch.Tensor:
    """Frozen SWE_Simple pilot: projected sorted quantile differences.

    The published SWE method uses a reference set and slicing projections;
    this function is intentionally static: reference and projections are
    supplied artifacts, never trained on the VCC loss.
    """
    x = residue_states.float() @ projections.T
    r = reference.float() @ projections.T
    n = max(x.shape[0], r.shape[0])
    q = torch.linspace(0, 1, n, device=x.device)
    xs = torch.sort(x, dim=0).values
    rs = torch.sort(r, dim=0).values
    xi = torch.linspace(0, 1, xs.shape[0], device=x.device)
    ri = torch.linspace(0, 1, rs.shape[0], device=x.device)
    def interp(values, grid, query):
        idx = torch.searchsorted(grid, query).clamp(1, len(grid) - 1)
        left, right = idx - 1, idx
        alpha = (query - grid[left]) / (grid[right] - grid[left]).clamp_min(1e-12)
        return values[left] + alpha * (values[right] - values[left])
    xq = torch.stack([interp(xs[:, j], xi, q) for j in range(xs.shape[1])], dim=1)
    rq = torch.stack([interp(rs[:, j], ri, q) for j in range(rs.shape[1])], dim=1)
    return (xq - rq).mean(dim=0)
