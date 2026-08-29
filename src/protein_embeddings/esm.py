"""Reproducible, inference-only ESM2 residue extraction.

This module deliberately operates on unique sequence hashes and never edits
the frozen biological mapping. Long proteins are handled by residue-level
reconstruction, not by averaging chunk-level protein vectors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence
from contextlib import nullcontext
from pathlib import Path

import torch


@dataclass(frozen=True)
class ESMConfig:
    model_name: str = "esm2_t33_650M_UR50D"
    layer: int = 33
    window_size: int = 1022
    overlap: int = 128
    weighting: str = "triangular"
    residue_dtype: str = "float32"


def load_model(checkpoint: str | None = None, device: str | None = None):
    """Load the official ESM2 checkpoint in evaluation/inference mode."""
    import esm

    if checkpoint:
        # Load the model state directly: ordinary residue/attention extraction
        # does not require FAIR's separate contact-regression weights.
        model_data = torch.load(checkpoint, map_location="cpu", weights_only=False)
        model, alphabet = esm.pretrained.load_model_and_alphabet_core(
            Path(checkpoint).stem, model_data, None
        )
    else:
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    target = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(target).eval()
    return model, alphabet, target


def residue_embeddings(model, alphabet, sequence: str, layer: int = 33,
                       device: str = "cpu", autocast_dtype=None) -> torch.Tensor:
    """Return exactly ``len(sequence) x hidden_dim`` residue states."""
    converter = alphabet.get_batch_converter()
    _, _, tokens = converter([("sequence", sequence)])
    tokens = tokens.to(device)
    with torch.inference_mode():
        if autocast_dtype is not None and device.startswith("cuda"):
            with torch.autocast(device_type="cuda", dtype=autocast_dtype):
                out = model(tokens, repr_layers=[layer], return_contacts=False)
        else:
            out = model(tokens, repr_layers=[layer], return_contacts=False)
    states = out["representations"][layer][0]
    # ESM alphabet places BOS at index 0 and EOS immediately after residues.
    residue = states[1:len(sequence) + 1]
    if residue.shape[0] != len(sequence):
        raise RuntimeError(f"residue indexing failed: expected {len(sequence)}, got {residue.shape[0]}")
    if not torch.isfinite(residue).all():
        raise FloatingPointError("non-finite ESM residue representation")
    return residue.detach().to("cpu", dtype=torch.float32)


def residue_embeddings_and_attention(model, alphabet, sequence: str, layer: int = 33,
                                     device: str = "cuda", autocast_dtype=None):
    """Return residue states and raw ESM attention heads including BOS/EOS."""
    converter = alphabet.get_batch_converter()
    _, _, tokens = converter([("sequence", sequence)])
    tokens = tokens.to(device)
    with torch.inference_mode():
        context = (torch.autocast(device_type="cuda", dtype=autocast_dtype)
                   if autocast_dtype is not None and device.startswith("cuda") else nullcontext())
        with context:
            out = model(tokens, repr_layers=[layer], need_head_weights=True,
                        return_contacts=False)
    states = out["representations"][layer][0, 1:len(sequence) + 1]
    attentions = out["attentions"][0]  # layers x heads x tokens x tokens
    if states.shape[0] != len(sequence) or not torch.isfinite(states).all():
        raise RuntimeError("invalid residue output from ESM")
    return states.detach().float().cpu(), attentions.detach().float().cpu()


def residue_embeddings_and_parti_attention_streaming(model, alphabet, sequence: str,
                                                      device: str = "cuda", autocast_dtype=None):
    """Run ESM2 while reducing PaRTI attention immediately, without retention.

    This is algebraically equivalent to the official max-over-heads followed
    by max-over-layers attention reduction; it avoids retaining all 33 attention
    tensors simultaneously.
    """
    converter = alphabet.get_batch_converter()
    _, _, tokens = converter([("sequence", sequence)])
    tokens = tokens.to(device)
    padding_mask = tokens.eq(model.padding_idx)
    with torch.inference_mode():
        x = model.embed_scale * model.embed_tokens(tokens)
        if model.token_dropout:
            x.masked_fill_((tokens == model.mask_idx).unsqueeze(-1), 0.0)
            mask_ratio_train = 0.15 * 0.8
            src_lengths = (~padding_mask).sum(-1)
            mask_ratio_observed = (tokens == model.mask_idx).sum(-1).to(x.dtype) / src_lengths
            x = x * (1 - mask_ratio_train) / (1 - mask_ratio_observed)[:, None, None]
        x = x * (1 - padding_mask.unsqueeze(-1).type_as(x))
        x = x.transpose(0, 1)
        pmask = None if not padding_mask.any() else padding_mask
        running = None
        context = (torch.autocast(device_type="cuda", dtype=autocast_dtype)
                   if autocast_dtype is not None and device.startswith("cuda") else nullcontext())
        with context:
            for layer in model.layers:
                x, attn = layer(x, self_attn_padding_mask=pmask, need_head_weights=True)
                layer_max = attn[:, 0].max(dim=0).values
                running = layer_max if running is None else torch.maximum(running, layer_max)
            x = model.emb_layer_norm_after(x).transpose(0, 1)
    states = x[0, 1:len(sequence) + 1].detach().float().cpu()
    if states.shape[0] != len(sequence) or not torch.isfinite(states).all():
        raise RuntimeError("invalid streaming ESM residue output")
    return states, running.detach().float().cpu()


def chunk_starts(length: int, window_size: int, overlap: int) -> list[int]:
    if length <= window_size:
        return [0]
    if not 0 < overlap < window_size:
        raise ValueError("overlap must be positive and smaller than window_size")
    stride = window_size - overlap
    starts = list(range(0, max(1, length - window_size + 1), stride))
    final = length - window_size
    if starts[-1] != final:
        starts.append(final)
    return starts


def reconstruction_weights(size: int, overlap: int, kind: str) -> torch.Tensor:
    """Deterministic positive taper; edges remain weight 1 for safe coverage."""
    if kind == "uniform":
        return torch.ones(size)
    x = torch.linspace(-1.0, 1.0, size)
    if kind == "triangular":
        return (1.0 - x.abs()).clamp_min(1.0 / max(size, 1))
    if kind in {"cosine", "hann"}:
        return (0.5 * (1.0 + torch.cos(torch.pi * x))).clamp_min(1e-3)
    raise ValueError(f"unknown weighting: {kind}")


def reconstruct_residue_embeddings(chunks: Sequence[torch.Tensor], starts: Sequence[int],
                                   length: int, overlap: int, weighting: str = "triangular") -> torch.Tensor:
    """Reconstruct full-length states by deterministic weighted residue averaging."""
    if len(chunks) != len(starts):
        raise ValueError("chunks and starts must have equal length")
    hidden = chunks[0].shape[1]
    total = torch.zeros((length, hidden), dtype=torch.float32)
    weights = torch.zeros(length, dtype=torch.float32)
    for chunk, start in zip(chunks, starts):
        w = reconstruction_weights(chunk.shape[0], overlap, weighting)
        stop = start + chunk.shape[0]
        total[start:stop] += chunk.float() * w[:, None]
        weights[start:stop] += w
    if (weights == 0).any():
        raise RuntimeError("uncovered residue during reconstruction")
    return total / weights[:, None]


def extract_sequence(model, alphabet, sequence: str, config: ESMConfig, device: str):
    """Extract direct or reconstructed residue states without silent truncation."""
    starts = chunk_starts(len(sequence), config.window_size, config.overlap)
    chunks = [residue_embeddings(model, alphabet,
                                  sequence[s:s + config.window_size], config.layer, device)
              for s in starts]
    if len(starts) == 1 and len(sequence) <= config.window_size:
        return chunks[0], {"chunk_count": 1, "starts": starts}
    return reconstruct_residue_embeddings(chunks, starts, len(sequence), config.overlap,
                                          config.weighting), {"chunk_count": len(starts), "starts": starts}
