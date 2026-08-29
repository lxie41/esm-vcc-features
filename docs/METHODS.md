# Methods and project design

## Scope

This repository owns gene/protein mapping, protein sequence handling, frozen ESM2 embedding extraction, pooling and representation generation, embedding storage, QC, provenance, and the downstream handoff contract.

It does not own single-cell preprocessing, flow matching, perturbation-model training, cell-state prediction, loss functions, or downstream model architecture. Biological execution has not started.

## Planned workflow

```text
Annotation sources
    ↓
Gene → transcript → protein mapping
    ↓
Validated protein sequences
    ↓
Frozen ESM2 inference
    ↓
Individual isoform representations
    ↓
Optional gene-level pooling
    ↓
Storage, QC, and provenance
    ↓
Downstream lookup interface
```

The repository stays flat until implementation complexity justifies further modules. Scripts are thin entry points; production code belongs in `src/protein_embeddings/`; notebooks are not part of the current scaffold.

## Mapping strategy (VCC 2026 validation panel)

The local VCC universe is `data/gene_names.csv`, not the perturbation-target
file: it contains 18,533 unique values in the `gene_name` column. The 300
unique values in `data/pert_counts.csv` are annotated later as perturbation
targets. The original VCC string is always retained as `vcc_gene_name`.

The primary mapping backbone is GENCODE human release 50 on GRCh38.p14,
corresponding to Ensembl release 116. Use the comprehensive ALL-region GTF for
gene/transcript/protein relationships and the matching GENCODE protein-coding
translation FASTA for amino-acid sequences. Ensembl is not a second parallel
annotation universe: its release-116 metadata is considered represented by the
GENCODE/Ensembl IDs and attributes, with a separate Ensembl download/API used
only if a required canonical/APPRIS field is absent from the frozen GTF.

MANE GRCh38 v1.5 is an auxiliary representative/cross-reference layer. Store
MANE Select and MANE Plus Clinical flags plus RefSeq transcript/protein IDs;
never filter the GENCODE isoform universe to MANE. UniProtKB is a secondary
protein-level QC and cross-reference source, preferably restricted to the
human reference proteome and reviewed Swiss-Prot records. RefSeq is accessed
through MANE first; a separate full RefSeq annotation release is not required
unless a classified GENCODE/MANE failure needs fallback review.

The hierarchy is therefore:

```text
VCC string → deterministic symbol/accession resolution → GENCODE/Ensembl 50/116
          → all eligible translated transcripts → distinct amino-acid sequences
          ├─ MANE v1.5 representative and RefSeq cross-references
          ├─ UniProt reviewed/canonical cross-references and sequence QC
          └─ RefSeq/NCBI targeted fallback for unresolved cases
```

## Mapping and isoforms

The future mapping is an explicit graph:

```text
gene symbol → stable gene ID → transcript ID → protein ID → amino-acid sequence
```

A gene is not a protein. Preserve zero, one, or many transcripts and protein products, including multiple transcripts that produce identical sequences. Handle historical symbols, missing protein IDs, non-protein-coding genes, pseudogenes, ambiguous symbols, and differing canonical definitions. Mapping failures must be classified and visible.

For this project, an ESM-level `protein_isoform` means one distinct valid
amino-acid sequence produced by one or more eligible translated transcripts
(Option B). The transcript table remains one row per transcript/protein link,
so identical translations are not lost. Create `protein_sequence_id` from the
lowercase-free, validated amino-acid string and a SHA-256
`protein_sequence_hash`; run future ESM inference once per unique hash.

Include normal protein-coding transcripts with a valid translation. Retain
nonsense-mediated-decay, non-stop-decay, retained-intron, incomplete-CDS, and
other edge biotypes in the audit mapping, but set `include_for_esm=false`
unless the GENCODE record supplies a complete, valid translation and the
policy review explicitly accepts it. Pseudogenes, non-coding genes, and
transcripts without a protein ID/translation are not normal ESM candidates.
Mitochondrial protein-coding genes are included when a valid GENCODE translation
exists and are flagged separately. Immunoglobulin/TCR, readthrough, PAR, and
patch-located records are retained with explicit flags rather than silently
discarded.

MANE Select and Ensembl Canonical are selection policies, not reasons to discard
other isoforms. Individual isoform embeddings must be retained before any gene
pooling. Canonical convenience vectors are derived later and never replace the
long-form mapping.

## Identifier normalization and evidence

Normalization is deterministic and non-destructive: preserve the exact VCC
string, trim only surrounding whitespace, detect Ensembl-like IDs separately,
and split stable IDs from version suffixes only into additional fields. Do not
rewrite symbols in place. Resolve symbols through the frozen GENCODE/Ensembl
gene-name attributes, then HGNC/Ensembl aliases if needed; retain the source
(`exact_symbol`, `alias`, `stable_id`, `accession`, or `ambiguous`). Duplicate
or conflicting resolutions become `ambiguous_gene_identifier`, never an
arbitrary choice. Accession-like VCC values such as `AC118549.1` are treated as
possible stable/accession identifiers and require explicit annotation evidence.

Use transparent evidence flags rather than a numeric confidence score:
`in_gencode`, `in_ensembl`, `mane_select`, `mane_plus_clinical`,
`uniprot_reviewed`, `uniprot_match`, `refseq_match`, and
`sequence_crossvalidated`.

Sequence reconciliation is asymmetric and auditable: GENCODE is the primary
ESM sequence source; MANE matching is strong evidence when the paired
transcript/protein agrees; UniProt and RefSeq sequences are compared by hash
and never silently substituted. A mismatch is a warning when an alternative
database has a known release/isoform difference, `manual_review` when the same
claimed product disagrees unexpectedly, and `mapping_failure` only when no
valid primary sequence can be established.

## Planned mapping schemas

The long-form mapping should contain only fields needed for lineage and QC:
`vcc_gene_name`, `gene_id`, `gene_symbol`, `gene_symbol_source`,
`transcript_id`, `transcript_version`, `transcript_biotype`, `protein_id`,
`protein_version`, `protein_sequence_id`, `protein_sequence_hash`,
`protein_length`, `is_protein_coding`, `is_mitochondrial`, `is_mane_select`,
`is_mane_plus_clinical`, `is_ensembl_canonical`, `refseq_transcript_id`,
`refseq_protein_id`, `uniprot_accession`, `uniprot_reviewed`,
`uniprot_canonical`, `uniprot_sequence_match`, evidence flags,
`sequence_source`, `annotation_release`, `include_for_esm`, `exclusion_reason`,
`mapping_status`, and `review_flag`. Store transcript-to-sequence links even
when multiple transcripts share one `protein_sequence_id`.

The gene-level QC summary contains `vcc_gene_name`, `protein_coding`,
`num_transcripts`, `num_unique_protein_sequences`, `mane_available`,
`uniprot_reviewed_available`, `esm_candidate_available`,
`perturbation_target`, `mapping_status`, and `review_flag`.

## Mapping execution result

The implemented GENCODE + MANE-only run produced exactly 18,533 summary rows.
It matched 17,974 genes to at least one GENCODE gene ID, identified 17,922
protein-coding genes, 17,943 genes with translated proteins, and 17,919 genes
with at least one valid conservative ESM candidate. There were 559 exact-name
misses, 1,311 multi-gene-ID cases, 36 non-protein-coding summaries, one
policy-review summary, and one translated-sequence-invalid summary.

The mapping contains 419,453 transcript rows, 357,518 translated transcript
records, and 226,301 distinct VCC-linked amino-acid sequences. 131,217
translated transcript/protein rows share a sequence hash with another row;
transcript and protein identities remain intact. 15,416 genes have multiple
unique sequences and 2,503 have exactly one. MANE Select is available for
17,888 genes, with one MANE reconciliation lacking a joined sequence. The
perturbation-target QC is 299 of 300 with a candidate; `TMEM104` is unresolved.

The 75 versioned/accession-like VCC identifiers identified by the exact QC rule
were all unmatched. Preliminary multi-ID classification found 662
patch/alternate-locus cases and 649 true-ambiguous/other cases. These are
reported, not collapsed. Full counts by transcript type and exclusion reason
are in `data/processed/mapping/mapping_qc.json`.

## ESM2 strategy

Current candidate: `esm2_t33_650M_UR50D`, approximately 650M parameters with expected hidden dimension 1280. Initial use is frozen inference only: no training or fine-tuning. ESM2-15B is outside the current local plan.

The chosen implementation/checkpoint is recorded in the Phase 1 resource
manifest, but model execution remains outside this mapping-design task.

No model was installed, downloaded, loaded, or run during preparation.

## Required resources and compatibility

The acquired snapshot, with release and SHA-256 recorded, contains GENCODE 50
`gencode.v50.annotation.gtf.gz`, GENCODE 50
`gencode.v50.pc_translations.fa.gz`, GENCODE RefSeq metadata, MANE v1.5
`MANE.GRCh38.v1.5.summary.txt.gz`, and a pinned UniProt human reference
proteome/ID-mapping export for QC. GENCODE 50 / Ensembl 116 / GRCh38.p14 / MANE
v1.5 is a coherent current snapshot: Ensembl release-116 documentation reports
GENCODE 50 and MANE v1.5. UniProt and RefSeq exports must be pinned to their
download release and are auxiliary, not replacements for the backbone.

The future enrichment plan is:

| Resource | Exact file/request | Approximate compressed size | Role |
|---|---|---:|---|
| GENCODE | `gencode.v50.annotation.gtf.gz` | tens of MB | required backbone |
| GENCODE | `gencode.v50.pc_translations.fa.gz` | about 10 MB | required ESM sequence source |
| GENCODE | `gencode.v50.metadata.RefSeq.gz` | a few MB | useful cross-reference |
| MANE | `MANE.GRCh38.v1.5.summary.txt.gz` | about 1 MB | required representative flags |
| UniProtKB | human reference proteome `UP000005640` FASTA, pinned release | tens of MB | optional QC/cross-reference |
| UniProtKB | matching ID-mapping TSV export for Ensembl/RefSeq/UniProt accessions | tens to hundreds of MB | optional QC/cross-reference |
| RefSeq/NCBI | targeted accession records only, if failures require them | case-dependent | optional fallback |

Do not use mutable `current` URLs in a released snapshot; record the resolved
URL, date, release, license, and SHA-256 for every acquired file.

No exact VCC coverage percentage is asserted before mapping. The honest
expectation is broad GENCODE symbol coverage, with lower coverage for accession
like, withdrawn/alias, non-coding, pseudogene, and unresolved entries; the
mapping report must calculate percentages against all 18,533 VCC genes and
separately against the 300 perturbation targets.

## Storage and provenance

Future artifacts must represent genes, transcripts, protein isoforms, canonical designation, individual vectors, and metadata without forcing one vector per gene. A candidate is metadata in Parquet plus dense matrices in safetensors or NumPy, keyed by stable row IDs; the final format remains open.

Every embedding should be traceable through gene → transcript → protein ID → sequence and source release → checkpoint → layer → pooling → final vector. A future manifest should include schema version, creation/tool versions, source records and checksums, model details, sequence policy, pooling, dimensions, counts, artifact paths, and run ID. Per-row metadata should include stable IDs, biotypes, canonical flags, sequence hash/length, embedding key, status, and failure reason.

## QC requirements

Future QC must report counts of genes, protein-coding genes, transcripts, unique protein sequences, multi-isoform genes, identical protein products, missing and ambiguous mappings, duplicate IDs, invalid sequences, and protein-length distributions. Embedding QC must check failures, NaN/Inf values, dimensions, row/key uniqueness, sequence-to-row consistency, reproducibility metadata, and missing-embedding semantics.

## Downstream handoff

The downstream team should not need ESM internals. A future implementation may expose:

```python
has_protein_embedding(gene)
get_gene_isoforms(gene)
get_isoform_embedding(protein_id)
get_canonical_embedding(gene)
get_gene_embedding(gene, representation="canonical")
```

Known non-protein-coding genes, unknown genes, and known genes without an available embedding must be distinguishable. Multiple isoforms remain enumerable. Missing vectors are never silently replaced with zeros. Returned data carries dimension, representation, source, version, and provenance.

## Development rules

Use configuration rather than hard-coded constants, structured logging, public type hints, biological docstrings, deterministic behavior where practical, and explicit provenance. Never silently truncate sequences, collapse isoforms, change annotation releases/checkpoints/layers, or start large computation. Add architectural changes to `DECISIONS.md` and unresolved items to `NOTES.md`.

## Refinement and final ESM input contract

The refinement snapshot adds official HGNC bulk resolution to the unchanged
GENCODE/MANE backbone. Resolution is exact and non-fuzzy: approved symbols,
previous symbols, aliases, withdrawn/merged symbols, and Ensembl stable IDs are
recorded with method and evidence. GENCODE contigs are categorized by locus;
an unambiguous primary reference locus is selected, with a MANE Select patch
exception. Other loci remain auditable and are not automatically primary ESM
inputs.

Eligibility is `primary`, `secondary`, `excluded`, or `review`. Valid ordinary
protein-coding records on the selected locus are primary; valid non-primary or
special translated records are secondary; missing/invalid/non-protein records
are excluded; ambiguous loci and sequences containing terminal/internal/
multiple `*` characters are review. Sequences are never silently modified.
Final ESM input is normalized by SHA-256 in
`data/processed/mapping_refined/esm_sequences.parquet`, with lineage in
`esm_sequence_relationships.parquet`.

## Future protein-language-model representation design

This section is design documentation only. The frozen mapping answers which
biological protein sequences belong to each VCC gene; representation answers
how each sequence is encoded; aggregation answers how multiple isoforms become
a gene feature; prediction evaluates usefulness for the perturbation task.
These questions must remain separate.

### Two representation levels

For a protein of length `L`, ESM2 produces residue states
`H ∈ R^(L×1280)` for `esm2_t33_650M_UR50D`. A separate residue-to-protein
operation is required. The initial static candidates are:

- R0 Mean: the 1,280-dimensional residue mean; simple, cheap, and the published
  baseline, but it loses heterogeneity, local organization, and distributional
  shape.
- R1 Mean+SD: concatenate per-dimension residue mean and explicitly defined
  population standard deviation for 2,560 dimensions. This preserves first and
  second marginal moments without learned parameters, but not multimodality,
  covariance, or domain organization.
- R2 SWE/distribution-aware: a strong primary candidate treating residue states
  as a set/distribution rather than immediately reducing them to moments. A
  peer-reviewed study supports sliced-Wasserstein aggregation as a fixed-length,
  distribution-aware alternative to average pooling, particularly for longer
  proteins ([NaderiAlizadeh and Singh, 2025](https://pubmed.ncbi.nlm.nih.gov/40170888/)).
  This supports SWE as a serious candidate, not a claim that it is best for VCC.
  Its exact formulation, dimensionality, normalization, streaming behavior,
  and compute cost remain open. SWE should not be described as assigning scalar
  importance weights to active sites; that is a different residue-importance
  method.

Advanced residue attention, local/window-aware pooling, set pooling, and
domain-aware pooling remain future research. Learned residue attention belongs
to downstream modeling because its weights require a task objective.

After protein-level vectors exist, isoform-to-gene aggregation is a separate
level. Static candidates are G0 Isoform Mean, G1 Isoform Mean+Variance, and G2
MANE-Anchored Diversity. G2 means a MANE representative plus explicit
alternative-isoform deviations; it does not discard non-MANE isoforms. MANE
absence must not trigger a silent fallback. The exact anchor/fallback policy is
open.

The conceptual design is therefore 3 protein representations × 3 gene
aggregations, or up to 9 deterministic gene candidates. It does not mean that
each isoform has nine vectors. Dimensionality varies: Mean→Isoform Mean is
1,280, Mean+SD→Isoform Mean is 2,560, and adding isoform variance doubles the
corresponding protein dimension. Higher dimensionality is not evidence of
better downstream performance.

### Feature-bank and validation principles

Future work should preserve a protein feature bank keyed by exact sequence hash
and retain supporting VCC/ENSG/ENST/ENSP/MANE relationships. A separate static
gene feature bank may expose deterministic ablation candidates, but per-isoform
features must remain available for Set Transformers, learned isoform attention,
or cell-conditioned weighting without rerunning ESM. Learned weights must not be
interpreted as measured isoform abundance because the VCC inputs are generally
gene-level.

Intrinsic validation can assess isoform distinguishability, sequence-change and
domain/local-change sensitivity, representation collapse, neighborhood
preservation, compression behavior, stability to controlled edits, and
gene-level isoform-diversity preservation. These diagnostics identify collapse
or instability; they cannot determine the best VCC representation without the
same downstream split, model, training procedure, and controlled ablation.

Long proteins are a mandatory gate before full extraction. Silent first-1,022
residue truncation is prohibited. Candidate future approaches include
overlapping/non-overlapping windows, chunk aggregation, domain-aware splitting,
or a longer-context encoder; none is selected here. The likely storage pattern
is one ESM forward pass with temporary residue states, immediate computation of
multiple static features, compact per-isoform persistence, and selective raw
residue caching. Future artifacts must version checkpoint, layer, normalization,
long-protein policy, pooling definition, dtype, dimensions, code version, and
generation date; representations must never be overwritten in place.

## CUDA and static-representation pilot update (2026-08-29)

The existing FAIR checkpoint was not redownloaded. Its local SHA-256 remains
`EA9D0522B335A8778DEA6535A65301F10208DECE28CD5865482B0B1FC446168C`. The
separate Python 3.12 `.esm-cuda-venv` uses PyTorch 2.4.1+cu124 and fair-esm
2.0.0. CUDA inference succeeded on an NVIDIA GeForce RTX 3070 Ti; PyTorch
reported 8,589,410,304 bytes of device memory.

GPU smoke tests at 100/300/500/800/1022 residues returned exact residue shapes
and finite values. FP16 is selected only through CUDA autocast with an FP32
model; full-model `.half()` is prohibited because fair-esm attention produced
a mixed-dtype error. At length 500, FP16 versus FP32 pooled Mean had cosine
1.0, mean absolute error 5.04e-05, and relative L2 error 3.32e-04; Mean+SD
relative L2 error was 3.00e-04. Reconstruction and accumulation remain
float32.

The expanded reconstruction benchmark contains 30 deterministic sequences
spanning 86–979 residues in five strata, with 90 cases per weighting across
overlaps 64/128/256. Using a forced 512-residue stress window, aggregate mean
cosine was 0.995137 uniform, 0.996742 triangular, and 0.996802 cosine; mean
per-residue L2 error was 0.611453, 0.526066, and 0.519321. The selected rule
remains 1022-window/128-overlap/triangular/float32; the stress window is not a
lossless native-context proof.

The deterministic `[Mean, SD, Q25, Q50, Q75]` candidate is streaming but
6400-dimensional versus 2560 for Mean+SD. On the 30-protein intrinsic pilot,
norm-length correlations were 0.181, 0.374, and 0.286 for Mean, Mean+SD, and
the quantile candidate; effective ranks were 4.91, 4.54, and 4.47. Similar
PCA concentration and added storage do not justify a third static family in
this phase. The quantile family is therefore deferred; the final static bank is
Mean, Mean+SD, and Pool PaRTI as specified below.

## Final pre-production representation contract (2026-08-29)

Primary per-isoform features are P0 ESM2 residue Mean (1280), P1 ESM2
residue Mean+SD (2560 consumed dimension, physically mean[1280] and sd[1280]),
and P2 Pool PaRTI (1280). Pool PaRTI follows the official implementation:
all 33 ESM2 layers; per-layer maximum across 20 heads; maximum across layers;
weighted directed graph; NetworkX PageRank with alpha=0.85, tol=1e-6, and
max_iter=100; remove BOS and EOS PageRank nodes; renormalize residue weights;
and weighted-average final-layer residue states. Padding is excluded.

The implementation was verified against the official streaming reduction: the
running max over per-layer head-max matrices exactly matched the retained
attention reduction on a test sequence. PaRTI output is finite, deterministic,
1280-dimensional, and its residue weights sum to one. The 25-protein direct
pilot had mean effective weighted residues ~584, maximum individual weight
~0.027, mean position correlation ~0.021, and mean-vs-PaRTI cosine ~0.990.

For long proteins, residue Mean/SD uses the frozen 1022/128/triangular/float32
reconstruction rule. PaRTI uses candidate A: run PaRTI per overlapping chunk,
merge local residue weights with the same triangular taper, globally normalize,
and apply the resulting weights to reconstructed full-length H. On 20 direct
references and overlaps 64/128/256, candidate A mean final-vector cosine to
full PaRTI was 0.994384, 0.994438, and 0.994914; mean JS weight distance was
0.050790, 0.050621, and 0.050229. Overlap 128 is frozen as the balance of
agreement and redundant computation. Candidate C chunk-vector pooling was
retained only as a control; candidate B dense global graph construction is
not used because its quadratic memory cost is unsuitable for very long
proteins. Long-protein PaRTI is an approximation and cannot restore context
beyond the native ESM2 window.

The permanent handoff stores per-isoform sequence_hash, length, mean, sd,
parti, and QC/provenance metadata. Mean+SD is concatenated on consumption.
The default future gene operation is ordinary mean across valid protein
isoforms, producing G0/G1/G2 for Mean/Mean+SD/PaRTI. Isoform mean+variance and
MANE anchoring remain optional later aggregations; per-isoform features are
never discarded.

SWE/SWE-Simple is no longer a primary static feature. The published method and
official implementation describe learned/objective-dependent reference,
slicing, and final projection components. It remains optional future
task-adaptive pooling, as do learned residue attention, Set Transformers, and
cell-conditioned weighting. Quantile pooling remains deferred.

## Future ESM representation design (not implemented)

## ESM2 pilot and extraction contract (2026-08-29)

The frozen biological input is the 174,585-row sequence universe under
`data/processed/mapping_final/`. It is immutable; ESM failures are recorded in
a separate processing layer and never remove mapping rows. `X` means an
unknown/unspecified residue and `U` means selenocysteine; neither is a stop
symbol (`*`) and neither is substituted during extraction.

ESM2 is inference-only, final layer 33, with residue states only: BOS, EOS, and
padding are excluded by slicing token positions `1:1+L`. For `L <= 1022`, the
sequence is run directly. Longer sequences use overlapping windows, infer each
window, reconstruct every residue by deterministic weighted averaging, and
only then compute protein features. Chunk-level protein pooling is a negative
control, never the production method. Mean is the population residue mean;
Mean+SD concatenates the population standard deviation (`correction=0`).

Single-isoform future gene variance is the zero vector. MANE, when available,
is an anchor plus alternative-isoform deviations; when unavailable,
MANE-anchored output is unavailable and no medoid/canonical/longest/first
isoform substitution is permitted. A medoid is a real isoform minimizing total
distance, and its identity can differ across representation spaces; it is
future work only. SWE is pilot-only until its published static formulation is
reproduced and validated.

Mapping and representation remain separate contracts. Mean residue pooling is a
useful fixed-dimensional baseline, but preserves only a first moment and loses
variance, domain heterogeneity, rare functional regions, and distribution shape.
A second mean across isoforms further removes isoform variation.

The first future comparison should include residue mean versus mean plus
standard deviation (2,560 dimensions for ESM2-650M), with population standard
deviation explicitly defined over residue tokens only. Later ablations may
examine quantile or distribution-aware pooling. Isoform-level vectors should
remain individually available; canonical/MANE, isoform mean, isoform variance,
and learned set aggregation are separate future choices.

Proteins above the ESM2 context require an explicit future policy: overlapping
windows, non-overlapping windows, domain-aware chunks, or a longer-context model.
No silent truncation or chunking is implemented. Raw residue storage is
estimated at roughly 394 GiB FP32 / 197 GiB FP16 for the final universe, so the
likely default is streaming pooled/statistical representations with selective
raw-residue caching for method development.

### Future representation extensions

MANE-Anchored Diversity may later include local splice-difference features. For
each alternative isoform, a proper protein sequence alignment to MANE could
identify changed regions; the method could retain the MANE whole-protein vector,
global isoform diversity, and a separately pooled difference-region feature.
Simple string diff or LCS is not the assumed final method. Future alignment
metadata should include identity, coverage, gap-block count, largest gap, and
alignment quality/status. Poor-quality alignments must fall back to a general
isoform representation.

Intrinsic QC must test protein-length bias for Mean, Mean+SD, and SWE: embedding
norm versus length, major PCA components versus length, pairwise similarity versus
length difference, and representation statistics versus length. Detection leads
to magnitude assessment and comparison of raw and normalized alternatives, not
automatic normalization.

ESM2 layer choice is a future ablation. The final layer remains the default
baseline, while a manageable subset may compare final-only, intermediate layers,
simple layer averages, weighted fusion, and concatenation followed by dimension
reduction. Multi-layer extraction is not part of the first-pass full-universe
pipeline because it increases memory, storage, pooling complexity, and feature

### ESM2 pilot results and gate (2026-08-29)

The verified checkpoint is stored under `data/external/esm2/`, with provenance
in `reports/esm_pilot/model_provenance.json`. Five deterministic sequences
(276, 493, 501, 716, and 929 residues) were used for an artificial 512-residue
window test. Across 45 reconstructions, mean cosine was 0.996725 uniform,
0.997863 triangular, and 0.997840 cosine; mean per-residue L2 error was
0.447326, 0.385598, and 0.387805. Triangular/128 is the provisional rule.
Native 1,022-window CUDA validation is still required.

The smoke test produced exact `L x 1280` residue matrices and 1280/2560
Mean/Mean+SD vectors. CPU repeats were bitwise identical and all outputs were
finite. The feature-bank prototype uses hash-indexed NumPy matrices plus JSON
metadata; atomic completed-hash checkpoint support is implemented. SWE was not
generated because its project-specific static reference/projection parameters
are not frozen.

### A100 production execution mode (2026-08-29)

The production extractor supports length-aware dynamic batching for native-context
proteins (<=1022 residues). Rows are sorted by length within each execution
batch, constrained by a configurable total-token budget and maximum batch size,
and results are restored by sequence_hash before writing. BOS, EOS, and padding
are sliced independently for every sequence. Proteins over 1022 residues retain
the serial validated chunk/reconstruction path.

The frozen feature mathematics is unchanged: FP16 CUDA autocast is used for
model execution, while residue reconstruction, pooling, accumulation, and
stored outputs are float32. `scripts/benchmark_extraction.py` compares serial
and batched Mean, SD, and PaRTI outputs on a deterministic stratified subset.
The Colab A100 serial baseline for shard 0000 was 3,351.6 proteins/hour with
approximately 3.4/40 GiB observed GPU memory. A100 batch-budget measurements
must be recorded by running the benchmark harness on a separate subset; this
repository preparation does not rerun shard 0000 or launch remaining shards.

The local reference environment is RTX 3070 Ti, Python 3.12, PyTorch 2.4.1+cu124,
and fair-esm 2.0.0. The production target environment is Colab A100-SXM4-40GB,
Python 3.13, PyTorch 2.11.0+cu128, and fair-esm 2.0.0. The official checkpoint
is loaded with `torch.load(..., weights_only=False)` under PyTorch >=2.6 as a
checkpoint-package compatibility setting; it does not change model weights or
scientific definitions. Cross-environment numerical equivalence must be
reported from a fixed subset before the A100 execution mode is frozen.
