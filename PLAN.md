# Implementation Roadmap

## Current phase: Phase 6 — Production extraction (next explicit task)

**Goal:** build the auditable VCC gene → transcript → unique protein-sequence mapping.

The current phase is bounded pilot work only. It must not modify
`data/processed/mapping_final/`, train downstream models, generate gene-level
representations, or iterate over all 174,585 sequences. Its exit gate is a
versioned ESM2 contract, quantitative long-protein evidence, validated static
pooling, storage/resume prototype, and an explicit readiness decision.

Artifacts include `src/protein_embeddings/esm.py`, `pooling.py`, `storage.py`,
and `scripts/esm_pilot.py`; pilot outputs are confined to
`reports/esm_pilot/`. Long proteins are never silently truncated: residue
reconstruction precedes pooling. Window/overlap/weight/dtype/batching/storage
choices remain unfrozen until controlled ablation.

## Future phase template

Each phase below remains unexecuted. For each, the listed outputs, components, and validation criteria must be finalized before implementation.

## Phase 1 — Resource selection and local setup — COMPLETE

**Goal:** choose and freeze annotation, protein sequence, ESM2, and Arc baseline resources.

**Inputs:** project requirements and candidate sources. **Outputs:** frozen resource choices and provenance plan.

## Phase 1b — VCC mapping strategy design — COMPLETE

**Goal:** specify the VCC universe, mapping hierarchy, identifier rules, isoform definition, deduplication, inclusion/exclusion, conflict handling, schemas, and resource plan.

**Validation:** `gene_names.csv` and `pert_counts.csv` inspected; 18,533 unique VCC genes and 300 in-vocabulary perturbation targets confirmed; no mapping or ESM execution.

## Phase 2 — Annotation acquisition

**Goal:** acquire approved sources reproducibly. **Inputs:** Phase 1 decisions. **Outputs:** immutable raw files and source manifest. **Dependencies:** network, licenses, storage. **Components:** download/checksum tooling. **Validation:** checksums and release/build match. **Risks:** release drift and licensing. **Open questions:** source selection.

**Status:** COMPLETE. The three approved mapping resources are frozen in
`data/raw/`; UniProt, full RefSeq releases, ESM2, and embeddings were not
downloaded.

## Phase 3 — Gene/transcript/protein mapping — COMPLETE

**Goal:** build explicit one-to-many mappings. **Inputs:** approved annotations. **Outputs:** mapping tables and classified failures. **Validation:** completed in `data/processed/mapping/`; external fallback and ESM remain deferred. Outputs include `gene_protein_mapping.parquet`, `protein_sequences.parquet`, `gene_mapping_summary.csv`, and `mapping_qc.json`.

## Phase 3b — Mapping refinement and final ESM input definition — COMPLETE

HGNC complete and withdrawn bulk TSVs are frozen in `data/raw/hgnc/current/` and
recorded in the resource manifest. Exact approved, previous, alias,
withdrawn/merged, and stable-ID resolution is deterministic. GENCODE contigs
are classified by locus context, MANE Select patch records are explicit
exceptions, and amino-acid `*` characters are audited without modification.
Refined outputs are under `data/processed/mapping_refined/`; only `primary`
rows populate the normalized ESM input tables. ESM2 was not downloaded or run.

## Phase 3c — Pre-ESM readiness audit — COMPLETE

The read-only audit is complete under `data/processed/pre_esm_audit/`. It found
that `gene_mapping_summary_refined.csv` contains 18,395 rows rather than the
required 18,533 VCC rows, leaving 138 VCC identifiers absent from the gene-level
summary. This is reported but not repaired in the audit. ESM2 setup and
embedding extraction waited for a targeted correction and rerun of the
refinement outputs. The correction is complete in `data/processed/mapping_final/`.

## Phase 3d — Final mapping freeze gate — COMPLETE

The final summary is a left join over all 18,533 VCC identifiers. Sequence and
relationship tables reconcile with zero missing primary hashes, zero orphans,
zero hash conflicts, and zero unexpected non-primary loci. All five tests pass.
ESM representation design is documented only; ESM2 remains unstarted.

## Phase 4 — Representation design documentation — COMPLETE

The design preserves three static residue-to-protein candidates (Mean,
Mean+SD, SWE/distribution-aware) and three static isoform-to-gene candidates
(Isoform Mean, Isoform Mean+Variance, MANE-Anchored Diversity). Per-isoform
features remain the durable interface for future learned aggregation. Long
proteins, storage, versioning, and intrinsic validation are documented gates;
none is implemented.

The near-term representation contract is:

```text
protein sequence
    ↓
ESM2 residue representation
    ↓
explicit long-protein strategy
    ↓
Mean / Mean+SD / SWE
    ↓
per-isoform protein feature bank
    ↓
Isoform Mean / Isoform Mean+Variance / MANE-Anchored
```

Future extensions include alignment-based local MANE splice-difference
features, protein-length bias QC, and a small-subset ESM2 multi-layer ablation.
None is implemented or selected as the default pipeline.

## Future phases

1. Long-protein strategy research and extraction-contract freeze.
2. ESM2 setup and small-scale validation.
3. Static per-isoform feature generation.
4. Intrinsic representation QC.
5. Optional deterministic gene feature bank.
6. Downstream ablation and learned aggregation.

## Phase 4 — Small-scale ESM validation

**Goal:** validate approved ESM execution and representation semantics. **Inputs:** approved sequences and runtime. **Outputs:** validation artifacts and decisions. **Dependencies:** Phases 1–3. **Components:** loader, tokenizer, pooling, resource logging. **Validation:** dimensions, determinism, no truncation. **Risks:** VRAM limits and API differences. **Open questions:** layer, dtype, long-sequence policy.

## Phase 5 — Human-scale ESM extraction

**Goal:** extract reproducible isoform vectors for the approved human universe. **Inputs:** validated sequences and Phase 4 configuration. **Outputs:** vectors, metadata, failures, manifest. **Dependencies:** Phases 1–4. **Components:** resumable batching and storage writer. **Validation:** completeness, dimensions, provenance, resource limits. **Risks:** failures, disk use, interrupted runs. **Open questions:** residue retention.

## Phase 6 — Derived representations and pooling

**Goal:** create explicitly selected gene-level products without deleting isoform vectors. **Inputs:** isoform vectors and mapping. **Outputs:** versioned derived vectors. **Dependencies:** Phase 5 and downstream requirements. **Components:** canonical/mean/variability-aware pooling. **Validation:** reproducible aggregation and lineage. **Risks:** information loss and policy confusion. **Open questions:** preferred representations.

## Phase 7 — QC and benchmarking

**Goal:** assess biological mappings, embeddings, and separately-provenanced baselines. **Inputs:** artifacts and manifests. **Outputs:** QC and benchmark reports. **Dependencies:** Phases 2–6 and legal approval. **Components:** QC suite and resource profiling. **Validation:** thresholds and review sign-off. **Risks:** incomparable baselines. **Open questions:** Arc access and comparison criteria.

## Phase 8 — Downstream handoff

**Goal:** publish a stable contract for the flow-matching team. **Inputs:** accepted artifacts and reports. **Outputs:** versioned metadata/tensors, lookup interface, handoff guide. **Dependencies:** Phases 1–7. **Components:** packaging and compatibility checks. **Validation:** consumer integration tests and missing-data semantics. **Risks:** dimension/version mismatch. **Open questions:** one vector versus multiple isoforms.

The implementation phase remains planned and must wait for review of the Phase 1b design.

## Phase 4b — ESM2 pilot validation — COMPLETE

The bounded smoke and reconstruction pilots completed under
`reports/esm_pilot/`. Model provenance, token indexing, X/U behavior,
Mean/Mean+SD dimensions, deterministic CPU repeats, and a hash-keyed storage
prototype passed. The provisional long-protein rule is 1022/128/triangular/
float32 after residue reconstruction. Full extraction is blocked until SWE is
authoritatively configured and the intended CUDA environment is verified. No
full extraction or gene-level aggregation was started.

CUDA / long-protein / representation follow-up: CUDA inference, FP16 autocast
stability, the expanded reconstruction benchmark, and a deterministic quantile
candidate pilot are complete. Static bank decision is Mean plus Mean+SD only;
SWE remains future learned pooling. Full extraction is still explicitly
unstarted because the native production profile and final SWE/feature contract
need no further ambiguity before scale-up.

## Final pre-production validation — COMPLETE

Final static bank: Mean, Mean+SD, and official Pool PaRTI. Default gene
aggregation: ordinary mean across valid isoforms. Long PaRTI: local PaRTI per
chunk, triangular merge, global normalization, weighted reconstructed H.
SWE/quantile remain deferred. The resumable binary feature schema and Colab
launcher are prepared; the 2-row shard dry run and resume check passed. Full
extraction remains prohibited in this task and is the next explicit phase.

## Production throughput preparation — IN PROGRESS

- Dynamic native-context batching is implemented with configurable token and
  batch-size limits; long proteins retain the frozen serial reconstruction path.
- A separate serial-vs-batched benchmark harness is prepared for a deterministic
  100–300 protein subset. It must be run on Colab A100 before freezing the A100
  execution budget; shard 0000 is not a benchmark input.
- The manifest-driven resumable Colab runner is implemented and validates local
  outputs before Drive publication.
- Production extraction has not started. The next task may run the benchmark
  and, only after its numerical-equivalence gate passes, launch all shards.
