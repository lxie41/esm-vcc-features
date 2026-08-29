# Protein Embeddings for Perturbation Modeling

Preparation repository for a focused protein-representation component of a single-cell perturbation prediction project. The intended future flow is gene → transcript(s) → protein isoform(s) → sequence → frozen ESM2 inference → traceable isoform and optional gene representations.

**PROJECT STATUS: MAPPING FROZEN; PRODUCTION SHARDS PREPARED; ESM EXTRACTION NOT STARTED**

The final GENCODE + MANE + HGNC mapping snapshot is frozen. ESM2, embeddings, external fallback, and downstream modeling have not started.

## Scope

This repository owns annotation handling, frozen gene/transcript/protein lineage, ESM2 inference, isoform-preserving feature storage, pooling, QC, provenance, and lookup interfaces. It does not own single-cell preprocessing, perturbation-model training, flow matching, or downstream architecture.

## Planned workflow

Annotations → mapping → sequences → ESM2 inference → individual isoform representations → optional derived gene representations → storage/lookup → downstream handoff.

## Repository layout

- `src/protein_embeddings/`: reusable ESM2, pooling, PaRTI, and storage code.
- `configs/`: configuration templates for local/Colab/HPC execution.
- `docs/`: scientific scope, architecture, data contracts, decisions, and operating rules.
- `data/`, `embeddings/`, and generated production artifacts: local-only and ignored by Git.
- `tests/`: unit and integration tests.

## Start here

Read `PLAN.md`, then `docs/METHODS.md`, `docs/DECISIONS.md`, and `docs/NOTES.md` before substantial changes. The files under `data/processed/mapping_final/` are frozen inputs and must not be rewritten.

## Development

Install the pinned dependencies from `requirements.txt` in a fresh environment. The Colab launcher at `colab/esm_full_extraction.ipynb` calls the same production script used locally; it does not contain a second scientific implementation.

## Colab production workflow

1. Create deterministic 2000-protein Parquet shards locally from the frozen sequence table:

   `python scripts/create_production_shards.py --input <frozen_sequences.parquet> --output-dir <shard_dir> --shard-size 2000`

2. Validate `manifest.json` and upload the shards plus the already verified ESM2 checkpoint to Google Drive.
3. Open `colab/esm_full_extraction.ipynb`, mount Drive, clone this repository, install pinned dependencies, verify CUDA and the checkpoint checksum, and select one input shard.
4. The notebook invokes `scripts/run_production_shards.py`, which calls the shared `scripts/extract_esm_features.py` with length-aware native-context batching, validates each local result, publishes only QC-valid binary Mean, SD, and PaRTI arrays plus metadata/checkpoint files to Drive, and skips completed shards.
5. Download or sync completed feature shards back to the local workspace and run final completeness/hash/QC merge checks.

The extractor defaults to `--max-tokens 8192 --max-batch-size 16`; benchmark
these controls on the target A100 with `scripts/benchmark_extraction.py` using a
separate 100–300 protein subset before changing the production budget. Use
`--disable-batching` for the serial reference comparison. The production
PageRank backend is the tensor implementation, which was numerically equivalent
to NetworkX on the fixed pilot; use `--pagerank-backend networkx` for reference
checks. Long proteins remain on the validated per-protein chunk/reconstruction
path.

Production extraction is intentionally not launched by the shard-preparation task. Model weights, frozen mapping tables, sequence shards, feature outputs, and checkpoints remain outside Git.
