# Notes, risks, and open questions

## ESM2 pilot execution (2026-08-29)

An isolated `.venv` contains NumPy, pandas, PyArrow, and `fair-esm==2.0.0`.
PyTorch installation and checkpoint download are execution prerequisites. The
host currently has no `nvidia-smi` executable; CUDA availability must be
measured, not assumed. This is an execution-environment issue, not a mapping
issue. Pilot deliverables include token checks, long-protein ablations,
Mean/Mean+SD integrity, SWE status, intrinsic QC, profiling, and a small
hash-keyed feature-bank test. Full extraction remains prohibited.

## Open questions

### Annotation (resolved for mapping design)

- GENCODE release 50, GRCh38.p14, comprehensive ALL-region GTF is primary.
- Matching GENCODE translation FASTA is primary protein source.
- MANE GRCh38 v1.5 is a separately cached flag/cross-reference layer; all GENCODE isoforms remain.
- UniProt is not a Phase 1 input; reconciliation is deferred.

- Resource acquisition is complete. The exact GENCODE 50 ALL GTF, GENCODE
  protein translation FASTA, and MANE v1.5 summary are frozen in
  `data/raw/resource_manifest.json`.
- Snapshot SHA-256 values are GTF `10623a32cc33a6a05e00016011cfb2af1cf72f50fd12587fee4bf4c61970f169`, protein FASTA `5d6408c3a1c22d864c96a55367cf5b69f1031dea1fe2cd6d03e366ae3fc56a45`, and MANE summary `d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c`.
- UniProt remains deferred until initial mapping coverage/QC is known. No ESM2
  checkpoint or embedding resource was downloaded in this phase.

- VCC vocabulary has 18,533 unique `gene_name` values; exact symbol/accession
  resolution and coverage are intentionally deferred to the reviewed mapping
  implementation.
- The extracted VCC files are present locally, but the original `control.zip`
  was not found in the repository; preserve the supplied manifest and hashes if
  the archive is later added.

### Isoforms (resolved for mapping design)

- Store all eligible translated transcript links and unique amino-acid sequence
  records.
- Deduplicate exact sequences by SHA-256 while retaining all transcript links.
- Represent MANE Select and Ensembl Canonical as flags; no canonical-only base
  mapping.

- Mapping execution is complete: 17,974 VCC genes matched at least one GENCODE
  gene ID; 17,919 have at least one valid conservative ESM candidate; 559 were
  not found by exact GENCODE `gene_name` matching.
- There are 226,301 unique VCC-linked protein sequences from 357,518 translated
  transcript records; 131,217 transcript/protein records collapse onto an
  existing sequence hash. 15,416 genes have multiple unique sequences.
- 17,888 genes have MANE Select annotations; one MANE-linked transcript lacked
  a joined sequence. Multi-gene-ID cases classify preliminarily as 662
  patch/alternate-locus and 649 true-ambiguous/other based on GENCODE region
  names.
- Of 300 perturbation targets, 299 have an eligible candidate; `TMEM104` is
  unresolved by GENCODE + MANE alone.

### ESM (resource freeze exists; execution deferred)

- Official `fair-esm==2.0.1`, checkpoint `esm2_t33_650M_UR50D`; no checkpoint was downloaded or run.
- Layer/token policy and hardware validation remain a separate execution task.
- How should long proteins be handled?
- fp16 or fp32, and what length is practical on 8 GB VRAM?
- Should raw residue embeddings be retained?

### Storage and downstream

- Tensor and metadata formats, compression, and expected size?
- One vector per gene or multiple isoforms?
- Preferred dimension and projection ownership?
- Required missing-data semantics?

### Arc/STATE baseline (resolved in earlier resource phase; not used for mapping)

Potential baseline: human gene → selected protein isoform(s) → ESM-derived representation → gene-level lookup vector. It must remain separate from locally generated embeddings, with independent provenance, dimensions, model, pooling, annotation lineage, and license.

Arc SE-600M `protein_embeddings.pt` is cached at an immutable Hugging Face revision. Arc declares 19,790 entries and 5,120 dimensions. It remains a separate gene-level baseline; its precise construction, normalization, and missing-gene behavior are not inferred in Phase 1.

### Mapping follow-up

- 75 versioned/accession-like VCC identifiers under the exact QC pattern were
  all unmapped; no aliases were invented.
- FASTA contained 153 records with `*`; 14 affected VCC-linked records were
  excluded as invalid. No automatic character replacement was performed.
- Before ESM extraction, review `TMEM104`, the 559 exact-name misses, the 649
  ambiguous/other multi-ID cases, and special transcript categories. Do not
  download UniProt or RefSeq until that review authorizes targeted fallback.

## Design risks and review notes

Primary risks are annotation drift, differing canonical definitions, silent sequence truncation, Arc provenance/license uncertainty, residue-storage scale, GPU memory limits, and accidental Git inclusion of artifacts.

The main design tension—gene-level representations versus preserving isoforms—is resolved by treating gene vectors as optional derived products. Configuration TODOs remain unresolved rather than guessed, and the ESM runtime remains a future verification task.

This file is intentionally a living note set, not an accepted-decision record. Accepted decisions belong in `DECISIONS.md`.

## Refinement snapshot (2026-08-22)

HGNC bulk files are frozen; no per-gene API, UniProt release, full RefSeq
release, ESM2 checkpoint, or embedding extraction was used. Refined QC reports
18,533 VCC rows, 426,770 auditable transcript rows, 171,293 unique primary
ESM sequences, 10 ambiguous HGNC tokens, 77 unresolved tokens, and 14
multiple-stop sequences. `TMEM104` is rescued through HGNC previous-symbol
evidence to `SLC38A12`. Re-review remaining unresolved/accession-like and
ambiguous cases before any ESM run.

## Final mapping freeze (2026-08-22)

The final frozen outputs are under `data/processed/mapping_final/`. The prior
18,395-row summary bug was corrected in generation logic by left-joining the
master VCC vocabulary; all 18,533 identifiers now have one explicit state. The
3,292 missing primary hashes were valid HGNC-rescued GENCODE protein
translations omitted by the prior sequence-table scope and were recovered by
exact protein-ID/hash verification. Final ESM sequence count is 174,585.

Final state counts are: 18,343 resolved protein-coding, 53 resolved
non-protein-coding, 41 resolved pseudogene, 9 resolved without a valid
translation, 10 ambiguous-but-explained, and 77 unresolved-but-explained.
Pytest 9.1.1 passes all 5 tests. ESM2 remains unstarted.

## Future ESM representation open questions

- Exact SWE formulation, output dimension, reference/projection design,
  normalization, reproducibility, and streaming feasibility.
- Exact MANE-anchored diversity formula and fallback when MANE is unavailable.
- Long-protein strategy for the 15,000+ sequences above the standard context.
- Vector storage backend, dtype, normalization, and dimensionality matching.
- Intrinsic metrics and a controlled downstream ablation protocol.

The intended sequence of future work is: research the long-protein policy,
freeze the ESM extraction contract, generate static per-isoform candidates,
perform intrinsic representation QC, then generate optional static gene-level
features. Learned residue pooling and learned isoform aggregation remain
downstream model research.

Additional open design items:

- Proper protein alignment method and quality thresholds for local
  MANE/alternative-isoform difference features; no simple string diff/LCS
  assumption.
- Length-bias diagnostics for Mean, Mean+SD, and SWE, followed by an empirical
  decision about any normalization.
- Small-subset ESM2 layer ablation: final, intermediate, averaged, weighted, or
  reduced concatenated layers before considering scale-up.
- SWE reference/projection design, computational budget, and whether a frozen
  or task-trained variant is appropriate.

Peer-reviewed SWE evidence: [Aggregating residue-level protein language model
embeddings with optimal transport](https://pubmed.ncbi.nlm.nih.gov/40170888/).
The evidence supports SWE as a strong candidate, especially for longer
proteins, but does not establish it as the best VCC representation.

## Pre-ESM readiness audit (2026-08-22)

The read-only audit is stored in `data/processed/pre_esm_audit/`. Referential
integrity for the final sequence and relationship tables passed, and all 300
perturbation targets have at least one primary sequence. However,
`gene_mapping_summary_refined.csv` has 18,395 rows instead of the required
18,533; 138 VCC identifiers are absent from that gene-level summary. This is a
real blocker because the summary is part of the mapping contract. The audit did
not modify mappings, rescue identifiers, or sequences.

The audit also found 5,162 refined mapping records marked `primary` whose
sequence hashes are absent from `esm_sequences.parquet` (3,292 distinct hashes).
The final normalized table itself has no orphan relationships, but the refined
eligibility table and final sequence table are therefore not fully reconciled.

All final primary sequences are primary-reference-locus sequences; no
alternate, haplotype, scaffold, or patch sequences leaked into the final table.
The primary universe has 15,312 sequences above the official FAIR extraction
default of 1,022 residues, so a long-sequence policy remains an engineering
decision after the summary-integrity correction. Technical research and storage
estimates are in `pre_esm_audit/technical_notes.md`.

## ESM2 pilot execution (2026-08-29)

The official checkpoint download is 2,604,537,549 bytes with SHA-256
`EA9D0522B335A8778DEA6535A65301F10208DECE28CD5865482B0B1FC446168C`.
The 3,687-byte contact-regression auxiliary file is also present and hashed;
contacts were not requested. PyPI supplied fair-esm 2.0.0, not the requested
2.0.1, for this runtime. The host has PyTorch 2.4.1+cpu and no CUDA device.

The actual blocker is SWE formulation freeze. A secondary limitation is absent
GPU tooling, so native 1,022-window CUDA profiling remains pending. The
reconstruction experiment used an artificial 512-residue context and should
not be treated as native-context validation.

CUDA follow-up: PyTorch 2.4.1+cu124 sees an NVIDIA GeForce RTX 3070 Ti and
8,589,410,304 bytes of memory. `nvidia-smi.exe` is not available on PATH or at
the standard Windows locations, so CLI driver/version text remains unrecorded.
GPU inference itself succeeds. The expanded benchmark used 30 deterministic
sequences and a forced 512-residue stress window because direct references all
fit the native 1022-residue limit.

SWE is not a guaranteed static preprocessing feature under the published
formulation: reference/slicer/final-projection components have learned or
objective-dependent roles. No arbitrary surrogate parameters were invented.
The tested static quantile extension was deferred; Mean and Mean+SD remain the
only guaranteed static families.

Superseding final contract: Pool PaRTI is now the third primary per-isoform
feature, while SWE/quantiles remain deferred. The final static bank is Mean,
Mean+SD, and PaRTI; default gene aggregation is ordinary isoform mean.

Final PaRTI gate: direct 25-protein QC and 20-protein long-protein candidate A
benchmark completed. The 8 GB GPU is near memory capacity for retained
attention (about 8.17 GB allocated at the largest short pilot), so production
must use the validated streaming attention reduction. `nvidia-smi` CLI text is
still unavailable, but PyTorch CUDA inference is operational. The remaining
pre-production validation is pipeline-level dry-run/Colab launcher checking;
no biological or representation research is reopened.
