# Architecture decision log

## 2026-08-29 — ESM2 setup and pilot gate

Preserve the final mapping and every valid isoform; single-isoform population
variance is zero; MANE is an anchor rather than an isoform filter; no MANE
means no automatic anchor; and silent truncation is prohibited. X and U remain
unchanged, meaning unknown residue and selenocysteine respectively; `*` is
distinct stop-related content.

Use official FAIR `esm2_t33_650M_UR50D`, final layer 33, and residue-only
pooling. Direct inference is used through the supported 1,022-residue limit.
Long proteins use overlapping windows followed by residue-level reconstruction
and then whole-protein pooling. Uniform, triangular, and cosine/Hann weighting
with overlaps 64/128/256 are measured against direct references; selection is
evidence-based. The bounded pilot cannot launch full extraction.

The official source is the `facebookresearch/esm` checkpoint. The available
PyPI package on this host is `fair-esm==2.0.0`; requested `2.0.1` is not
published for the selected Python runtime, so this compatibility discrepancy
must remain visible in provenance. Model file, resolved source, date, size,
checksum, actual config, and alphabet behavior are recorded after download.

SWE remains a pilot-only, frozen/static gate. Its reference set and slicing
projections must be versioned artifacts; no VCC/downstream training is allowed.
No gene-level nine-way aggregation is implemented in this phase. Dense vectors
use hash-keyed arrays with metadata in JSON/Parquet, and resume tracks completed
hashes, atomic shards, failures, and retries.

## 2026-08-18 — Repository scope is protein representation

Status: Accepted

Decision: Own annotation-to-protein representation and handoff, not downstream single-cell modeling.

Reason: clear responsibilities and interfaces.

Alternatives considered: a unified perturbation-model repository.

Consequences: downstream fusion and flow matching remain external.

## 2026-08-18 — Preserve individual isoforms

Status: Accepted

Decision: retain isoform records and vectors before any gene pooling; mean pooling is only a baseline.

Reason: averaging loses information and aggregation is undecided.

Consequences: storage supports one-to-many relationships.

## 2026-08-18 — Frozen ESM2 inference initially

Status: Accepted

Decision: pretrained inference only; no training or initial fine-tuning. Candidate is ESM2 650M, not 15B, for the current workstation.

Consequences: runtime, layer, precision, and checkpoint still require verification.

## 2026-08-18 — Separate baseline provenance

Status: Accepted

Decision: Arc/STATE-style and locally generated embeddings have separate provenance and namespaces.

## 2026-08-18 — Production code under src

Status: Accepted

Decision: reusable code lives under `src/`; scripts are thin and notebooks exploratory.

Consequences: no biological pipeline modules are implemented in Phase 0.

## 2026-08-22 — VCC mapping backbone and release freeze

Status: Accepted

Decision: use GENCODE human release 50 / GRCh38.p14 / Ensembl release 116 as
the complete gene-transcript-protein backbone. Use its comprehensive GTF and
matching protein translation FASTA; do not build a second parallel Ensembl
universe.

Reason: it provides broad transcript annotation, stable Ensembl identifiers,
protein IDs, biotypes, and matching translations while preserving isoforms.

## 2026-08-22 — Representative and auxiliary sources

Status: Accepted

Decision: use MANE v1.5 to flag MANE Select and MANE Plus Clinical and retain
RefSeq IDs. Use UniProtKB only as a secondary reviewed/canonical
cross-reference and sequence-QC source. Use targeted RefSeq/NCBI fallback only
for classified failures; no separate full RefSeq universe is required now.

Consequence: disagreements remain visible and source provenance is retained;
MANE/UniProt never filter GENCODE alternatives.

## 2026-08-22 — ESM protein identity and deduplication

Status: Accepted

Decision: define one ESM protein isoform as one distinct valid amino-acid
sequence. Preserve every transcript-to-sequence relationship, but deduplicate
identical sequences by SHA-256 before future ESM inference.

Consequence: ESM runs once per unique sequence while transcript and gene lineage
remains queryable. No mapping or ESM execution is part of this decision task.

## 2026-08-22 — Mapping resource snapshot

Status: Accepted

Decision: freeze GENCODE 50 / GRCh38.p14 / Ensembl 116 comprehensive ALL
annotation as `gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz`, the matching
`gencode.v50.pc_translations.fa.gz`, and MANE GRCh38 v1.5
`MANE.GRCh38.v1.5.summary.txt.gz`. Do not download UniProt or a full RefSeq
release until after initial mapping coverage is known.

Provenance: exact URLs, file sizes, SHA-256 values, and structural validation
results are recorded in `data/raw/resource_manifest.json`. The approved VCC
CSV files remain unchanged at `data/gene_names.csv` and `data/pert_counts.csv`.

## 2026-08-22 — VCC mapping execution policy

Status: Accepted

Decision: use GENCODE GTF biotypes as authoritative for ESM eligibility; do not
equate FASTA presence with eligibility. Include only valid `protein_coding`
translations. Retain special translated categories with
`requires_policy_review`, preserve noncoding/unmapped records in the summary,
and deduplicate only by SHA-256 of the amino-acid sequence.

Result: the executed GENCODE + MANE-only mapping and all QC outputs are under
`data/processed/mapping/`. MANE enriches records and never filters alternative
isoforms. UniProt, RefSeq fallback, and ESM2 remain unstarted.

## 2026-08-22 — Deterministic HGNC and ESM-input refinement

Status: Accepted

Decision: freeze the official HGNC complete set and withdrawn TSV bulk files;
resolve only exact approved/previous/alias/withdrawn-or-merged/stable-ID
evidence; and leave ambiguous or unresolved identifiers for review. Classify
GENCODE loci by contig, select one primary reference locus when unambiguous,
allow a MANE Select patch exception, and preserve all other locus contexts in
the audit mapping. Preserve sequences containing `*` exactly and mark them
review.

Consequence: `data/processed/mapping_refined/esm_sequences.parquet` and
`esm_sequence_relationships.parquet` are the frozen pre-ESM contract; the long
refined mapping and QC files preserve secondary, excluded, and review records.

## 2026-08-22 — Final mapping freeze gate

Status: Accepted

Decision: create a separate `data/processed/mapping_final/` snapshot whose gene
summary is left-joined to the immutable 18,533-value VCC master vocabulary. The
final sequence table contains only hashes referenced by selected primary
relationships. HGNC-rescued GENCODE translations are included when their exact
sequence hash is recovered from the frozen GENCODE FASTA; no biological policy
was changed.

Result: 18,533 unique summary rows, 174,585 final primary sequences, zero
missing primary hashes, zero orphan sequences, zero hash conflicts, zero
unexpected non-primary loci, and 300/300 perturbation targets with primary
sequences. Unresolved and ambiguous cases remain explicitly classified.

## 2026-08-22 — SWE promoted to primary static candidate

Status: Accepted design direction; implementation deferred

Decision: retain Mean, Mean+SD, and SWE as the three primary static
residue-to-protein candidates. Peer-reviewed work supports sliced-Wasserstein,
distribution-aware aggregation of PLM residue embeddings, with especially
relevant evidence for longer proteins. SWE is not claimed to be best for VCC
and is not an active-site weighting method.

Consequence: future design and ablation plans treat SWE alongside Mean and
Mean+SD while leaving its exact formulation and downstream validation open.

## 2026-08-22 — Representation QC and layer extensions

Status: Accepted design direction; implementation deferred

Decision: future intrinsic QC must test length dependence without imposing
automatic normalization; future MANE-Anchored Diversity may add
alignment-based local splice-difference features with explicit quality
metadata; and ESM2 multi-layer comparisons belong on a small representative
subset before any scale-up. The final layer remains the default baseline.

No alignment, length normalization, multi-layer extraction, or representation
implementation is performed by this decision.

## 2026-08-22 — Representation design boundary

Status: Accepted design direction; implementation deferred

Decision: preserve the separation between residue→protein representation and
protein-isoform→gene aggregation. Planned static candidates are protein Mean,
Mean+SD, and SWE/distribution-aware; planned deterministic gene candidates are
Isoform Mean, Isoform Mean+Variance, and MANE-Anchored Diversity. Mean is a
baseline, not a presumed winner. Per-isoform features must remain available for
future learned aggregation, including Set Transformers, isoform attention, and
cell-conditioned weighting.

Long-protein handling is a required design gate; silent truncation is prohibited.
No pooling, SWE, chunking, learned aggregation, ESM2 download, or embedding
generation is part of this decision.

## 2026-08-29 — ESM2 pilot readiness outcome

The bounded pilot selected window 1022, overlap 128, triangular deterministic
taper, float32 reconstruction, final layer 33, residue-only pooling, and
population Mean/Mean+SD statistics. Selection is provisional until native CUDA
measurements are available. The hash-keyed NumPy plus metadata prototype is
usable for the next phase; no vectors were written to CSV and no gene-level
features were created.

Readiness is **NOT READY — SWE IMPLEMENTATION UNRESOLVED**. SWE_Simple is only
scaffolded with explicit reference/projection inputs; its project-specific
reference distribution, projection count/dimension, normalization, and
streaming/storage contract are not frozen. The host is also CPU-only, so the
requested RTX 3070 Ti/CUDA gate remains unverified.

CUDA follow-up: a separate Python 3.12 environment with PyTorch 2.4.1+cu124
successfully verified the RTX 3070 Ti. FP16 is permitted only through CUDA
autocast around an FP32 model; full-model `.half()` is rejected because it
caused a fair-esm mixed-dtype attention error. The expanded 30-protein pilot
kept the 1022/128/triangular/float32 rule. A tested 6400-dimensional quantile
family did not provide clear intrinsic benefit over Mean+SD, so the guaranteed
static bank is Mean plus Mean+SD only. SWE is moved to future task-adapted
pooling because its published formulation includes learned components.

## 2026-08-29 — Final pre-production representation freeze

The primary static per-isoform bank is frozen as P0 Mean (1280), P1 Mean+SD
(mean and SD stored separately, 2560 when concatenated), and P2 Pool PaRTI
(1280). PaRTI is reproduced from the official implementation using all 33
layers, per-layer max over 20 heads, max over layers, weighted directed graph,
NetworkX PageRank alpha 0.85/tol 1e-6/max_iter 100, BOS/EOS removal,
renormalization, and weighted final-layer residue pooling.

Long PaRTI is frozen as local chunk PaRTI weights merged with the frozen
triangular overlap taper and globally renormalized before weighting the
reconstructed H. Window 1022, overlap 128, and float32 are retained. The
streaming reduction is algebraically identical to the official retained
attention reduction and avoids holding all attention layers at once. PaRTI
cannot recover beyond-context interactions.

Default gene aggregation is ordinary mean across valid protein isoforms for
G0/G1/G2. Isoform variance and MANE-anchor aggregations are optional future
products. SWE/SWE-Simple, quantile pooling, and other learned pooling are not
primary precomputed features; SWE remains future task-adaptive pooling.

## 2026-08-29 — Production throughput execution contract

Native-context sequences may use dynamic length-aware batches under configurable
`max_tokens` and `max_batch_size` limits. The default implementation uses
`max_tokens=8192` and `max_batch_size=16` as an execution starting point only;
the A100 benchmark harness must select the measured production budget. Hash-keyed
records are restored to deterministic sequence_hash order before output, and
per-sequence BOS/EOS/padding exclusion is preserved.

No batching is applied to the validated >1022-residue chunk/reconstruction path.
FP16 CUDA autocast and float32 outputs/accumulation remain frozen. The runner
`scripts/run_production_shards.py` copies one manifest shard to local scratch,
invokes the shared extractor, validates all three matrices/checkpoint/counts,
publishes only QC-valid outputs to Drive, and skips already valid outputs.

The A100 shard-0000 result (3,351.6 proteins/hour, 0 failures, ~3.4/40 GiB)
is retained as the serial baseline. It is not rerun here. The A100 batched
throughput and cross-environment numerical-equivalence report remain an explicit
Colab benchmark gate, not an assumed result.
