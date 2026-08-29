# Phase 1 — Frozen resource snapshot

Status: **complete** (2026-08-22)

This phase freezes external inputs and model artifacts only. The resource
selection is retained here for provenance. The three approved mapping resources
are now acquired under `data/raw/`; actual checksums and validation results are
authoritative in `data/raw/resource_manifest.json`. No gene mapping, sequence
normalization, model loading, or embedding extraction was performed.

## Accepted decisions

1. **Genome and annotation:** human GRCh38.p14, GENCODE release 50. Use the
   comprehensive ALL-region annotation GTF as the authoritative source for
   gene → transcript → protein relationships. The GTF `gene_id`,
   `transcript_id`, and `protein_id` attributes are retained with version
   suffixes; the matching GENCODE protein-coding translation FASTA is the
   protein sequence source.
2. **MANE:** cache MANE GRCh38 v1.5 summary as a cross-reference and canonical
   selection flag. MANE Select is an optional downstream selection policy; it
   does not remove non-MANE GENCODE isoforms. MANE Plus Clinical is retained as
   a separate flag if present. The GENCODE universe remains primary.
3. **ESM2:** use the official FAIR reference implementation (`fair-esm` 2.0.1,
   MIT source license) and checkpoint `esm2_t33_650M_UR50D` (33 layers, 650M
   parameters, 1280 hidden dimensions, UR50/D 2021_04). The initial extraction
   contract is final-layer (layer 33), mean over residue tokens only, with BOS
   and EOS excluded. Long-sequence handling remains a Phase 4 validation
   decision; Phase 1 never truncates anything.
4. **Arc baseline:** cache the published Arc SE-600M human protein-embedding
   artifact `protein_embeddings.pt` at immutable Hugging Face revision
   `985310635c9d6d5a60c0eb106f836ceab29f0c99`. Treat it as a separate,
   gene-level baseline: 19,790 entries and 5,120 dimensions as declared by
   Arc's configuration. Do not merge it with local ESM2 outputs. Arc code is
   CC BY-NC-SA 4.0; model weights and outputs are under the Arc Research
   Institute State Model Non-Commercial License and acceptable-use policy.

## Resource manifest

All files below are downloaded by `scripts/phase1_snapshot.ps1`. The local
paths are relative to the repository root and are intentionally ignored by
Git. SHA-256 values are generated after download in
`data/raw/resource_manifest.json`.

| Resource | Exact version / purpose | URL | License / terms | Local path |
|---|---|---|---|---|
| GENCODE comprehensive GTF | Human release 50, GRCh38.p14, ALL regions; relationships | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz | GENCODE/Ensembl data-use terms; cite GENCODE | `data/raw/gencode/v50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz` |
| GENCODE translations | Release 50 matching protein-coding translation FASTA | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_translations.fa.gz | GENCODE/Ensembl data-use terms; cite GENCODE | `data/raw/gencode/v50/gencode.v50.pc_translations.fa.gz` |
| MANE summary | GRCh38 v1.5 cross-reference and selection annotations | https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz | NCBI/RefSeq data-use terms; cite MANE | `data/raw/mane/v1.5/MANE.GRCh38.v1.5.summary.txt.gz` |
| ESM2 checkpoint | `esm2_t33_650M_UR50D.pt`; official FAIR weights | https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt | FAIR ESM source MIT; retain upstream terms/citation | `data/external/phase1/models/esm2_t33_650M_UR50D.pt` |
| Arc baseline vectors | SE-600M human protein embeddings, immutable revision above | https://huggingface.co/arcinstitute/SE-600M/resolve/985310635c9d6d5a60c0eb106f836ceab29f0c99/protein_embeddings.pt | Arc model non-commercial license + AUP | `data/external/phase1/arc/protein_embeddings.pt` |

The ESM and Arc rows are future Phase 1 selections only; they were not
downloaded for the current mapping-resource acquisition task. The ESM/Arc
snapshot script was not run here.

## Reproduce / verify

From the repository root in PowerShell:

```powershell
pwsh -File .\scripts\phase1_snapshot.ps1
Get-FileHash .\data\external\phase1\* -Algorithm SHA256
```

The script is idempotent: existing files are not re-downloaded, and the
checksum file is regenerated. The downloaded artifacts are not committed to
Git. Before Phase 2, review the checksum file and licenses and preserve this
manifest alongside any released artifact bundle.

## Explicit non-actions

Phase 1 does **not** parse the GTF, join MANE records, map genes, deduplicate
sequences, load ESM/Arc tensors, or generate embeddings. Those actions begin
only in later phases.
