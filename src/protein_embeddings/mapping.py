"""Auditable VCC-to-GENCODE transcript and protein-sequence mapping."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

ID_RE = re.compile(r"^(ENS[GPT][0-9]+)(?:\.[0-9]+)?$")
ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
SPECIAL_TYPES = ("nonsense_mediated_decay", "non_stop_decay", "protein_coding_LoF", "polymorphic_pseudogene")


def stable_id(value: str | None) -> str | None:
    """Return an Ensembl stable ID while preserving the source elsewhere."""
    if not value:
        return None
    match = ID_RE.match(value)
    return match.group(1) if match else value


def parse_attributes(value: str) -> dict[str, str]:
    return {key: val for key, val in ATTR_RE.findall(value)}


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_accession_like(value: str) -> bool:
    return bool(re.match(r"^[A-Z]{2}\d+\.\d+$", value) or re.match(r"^H3\.[A-Z]$", value))


@dataclass
class TranscriptRecord:
    seqname: str
    gene_id: str
    gene_id_stable: str
    gene_name: str
    gene_type: str
    transcript_id: str
    transcript_id_stable: str
    transcript_name: str
    transcript_type: str
    protein_id: str | None
    protein_id_stable: str | None


def parse_gtf(path: Path) -> tuple[dict[str, dict[str, str]], list[TranscriptRecord], Counter[str]]:
    """Stream gene and transcript records from the compressed GENCODE GTF."""
    genes: dict[str, dict[str, str]] = {}
    transcripts: list[TranscriptRecord] = []
    transcript_types: Counter[str] = Counter()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            if fields[2] == "gene":
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_id_stable": stable_id(gene_id) or gene_id,
                    "gene_name": attrs.get("gene_name", ""),
                    "gene_type": gene_type,
                }
            elif fields[2] == "transcript" and attrs.get("transcript_id"):
                transcript_type = attrs.get("transcript_type", attrs.get("transcript_biotype", ""))
                transcript_types[transcript_type] += 1
                protein_id = attrs.get("protein_id") or None
                transcripts.append(TranscriptRecord(
                    fields[0], gene_id, stable_id(gene_id) or gene_id, attrs.get("gene_name", ""), gene_type,
                    attrs["transcript_id"], stable_id(attrs["transcript_id"]) or attrs["transcript_id"],
                    attrs.get("transcript_name", ""), transcript_type, protein_id, stable_id(protein_id),
                ))
    return genes, transcripts, transcript_types


def parse_fasta(path: Path) -> tuple[dict[str, dict[str, Any]], Counter[str], int]:
    """Stream translation records keyed by exact and stable protein IDs."""
    records: dict[str, dict[str, Any]] = {}
    unusual: Counter[str] = Counter()
    total = 0
    current_header: str | None = None
    sequence_parts: list[str] = []

    def finish(header: str | None, parts: list[str]) -> None:
        nonlocal total
        if not header:
            return
        fields = header[1:].strip().split("|")
        if len(fields) < 8:
            unusual["malformed_header"] += 1
            return
        protein_id, transcript_id, gene_id, _, _, transcript_name, gene_name, declared_length = fields[:8]
        sequence = "".join(parts).strip().upper()
        bad = sorted(set(sequence) - VALID_AA)
        for char in bad:
            unusual[f"amino_acid:{char}"] += 1
        record = {
            "protein_id": protein_id, "protein_id_stable": stable_id(protein_id) or protein_id,
            "transcript_id": transcript_id, "transcript_id_stable": stable_id(transcript_id) or transcript_id,
            "gene_id": gene_id, "gene_id_stable": stable_id(gene_id) or gene_id,
            "gene_name": gene_name, "transcript_name": transcript_name,
            "amino_acid_sequence": sequence, "protein_length": len(sequence),
            "declared_length": declared_length, "sequence_valid": bool(sequence) and not bad,
            "sequence_issue": ";".join(bad) if bad else ("empty_sequence" if not sequence else ""),
            "protein_sequence_hash": sequence_hash(sequence) if sequence else None,
        }
        records[protein_id] = record
        records.setdefault(record["protein_id_stable"], record)
        total += 1

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                finish(current_header, sequence_parts)
                current_header, sequence_parts = line, []
            elif current_header:
                sequence_parts.append(line.strip())
        finish(current_header, sequence_parts)
    return records, unusual, total


def parse_mane(path: Path) -> dict[str, dict[str, str]]:
    """Aggregate MANE rows by stable Ensembl transcript ID."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = stable_id(row.get("Ensembl_nuc", ""))
            if key:
                grouped[key].append(row)
    result: dict[str, dict[str, str]] = {}
    for key, rows in grouped.items():
        statuses = sorted({row.get("MANE_status", "") for row in rows if row.get("MANE_status")})
        result[key] = {
            "mane_status": ";".join(statuses),
            "is_mane_select": str(any(row.get("MANE_status") == "MANE Select" for row in rows)).lower(),
            "is_mane_plus_clinical": str(any("Plus Clinical" in row.get("MANE_status", "") for row in rows)).lower(),
            "refseq_transcript_id": ";".join(sorted({row.get("RefSeq_nuc", "") for row in rows if row.get("RefSeq_nuc")})),
            "refseq_protein_id": ";".join(sorted({row.get("RefSeq_prot", "") for row in rows if row.get("RefSeq_prot")})),
        }
    return result


def eligibility(record: TranscriptRecord, sequence_record: dict[str, Any] | None) -> tuple[bool, str, bool]:
    if not record.protein_id:
        return False, "no_translation", False
    if not sequence_record:
        return False, "sequence_not_in_fasta", False
    if not sequence_record["sequence_valid"]:
        return False, "invalid_sequence:" + sequence_record["sequence_issue"], False
    if record.transcript_type == "protein_coding":
        return True, "", True
    if record.transcript_type in SPECIAL_TYPES or record.transcript_type.startswith(("IG_", "TR_")):
        return False, "requires_policy_review:" + record.transcript_type, True
    return False, "requires_policy_review:" + (record.transcript_type or "missing_transcript_type"), True


def build_mapping(vcc_path: Path, perturbation_path: Path, gtf_path: Path, fasta_path: Path,
                  mane_path: Path, output_dir: Path) -> dict[str, Any]:
    """Build mapping, unique sequence, summary, and QC outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    vcc = pd.read_csv(vcc_path, dtype=str)["gene_name"].tolist()
    perturbations = set(pd.read_csv(perturbation_path, dtype=str)["target_gene"].tolist())
    if len(vcc) != 18533 or len(set(vcc)) != len(vcc):
        raise ValueError(f"Expected 18,533 unique VCC genes, got {len(vcc)} / {len(set(vcc))}")
    if not perturbations.issubset(set(vcc)):
        raise ValueError("Perturbation target is absent from VCC vocabulary")

    genes, transcripts, transcript_type_counts = parse_gtf(gtf_path)
    fasta, unusual_aa, fasta_record_count = parse_fasta(fasta_path)
    mane = parse_mane(mane_path)
    by_name: dict[str, set[str]] = defaultdict(set)
    for gene_id, gene in genes.items():
        if gene["gene_name"]:
            by_name[gene["gene_name"]].add(gene_id)
    transcripts_by_gene: dict[str, list[TranscriptRecord]] = defaultdict(list)
    for transcript in transcripts:
        transcripts_by_gene[transcript.gene_id].append(transcript)

    mapping_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    sequence_records: dict[str, dict[str, Any]] = {}
    status_counts: Counter[str] = Counter()
    vcc_transcript_type_counts: Counter[str] = Counter()
    multi_gene_context_counts: Counter[str] = Counter()
    accession_stats = Counter()

    for vcc_gene in vcc:
        gene_ids = sorted(by_name.get(vcc_gene, set()))
        is_target = vcc_gene in perturbations
        accession = is_accession_like(vcc_gene)
        if accession:
            accession_stats["total"] += 1
        all_transcripts = [t for gene_id in gene_ids for t in transcripts_by_gene.get(gene_id, [])]
        for transcript in all_transcripts:
            vcc_transcript_type_counts[transcript.transcript_type] += 1
        linked_sequences: set[str] = set()
        valid_candidates = 0
        translated_count = 0
        mane_select = False
        protein_coding_gene = any(genes[gene_id]["gene_type"] == "protein_coding" for gene_id in gene_ids)
        special_or_review = False
        invalid_or_missing = False

        for transcript in all_transcripts:
            sequence_record = None
            if transcript.protein_id:
                sequence_record = fasta.get(transcript.protein_id) or fasta.get(transcript.protein_id_stable or "")
                if sequence_record:
                    translated_count += 1
            include, reason, sequence_valid = eligibility(transcript, sequence_record)
            mane_meta = mane.get(transcript.transcript_id_stable, {})
            mane_select = mane_select or mane_meta.get("is_mane_select") == "true"
            special_or_review |= reason.startswith("requires_policy_review")
            invalid_or_missing |= bool(reason and reason not in {"no_translation", "sequence_not_in_fasta"})
            if sequence_record and sequence_record["protein_sequence_hash"]:
                h = sequence_record["protein_sequence_hash"]
                existing = sequence_records.get(h)
                if existing and existing["amino_acid_sequence"] != sequence_record["amino_acid_sequence"]:
                    raise ValueError(f"Sequence hash conflict for {h}")
                sequence_records.setdefault(h, {
                    "protein_sequence_hash": h,
                    "amino_acid_sequence": sequence_record["amino_acid_sequence"],
                    "protein_length": sequence_record["protein_length"],
                    "sequence_valid": sequence_record["sequence_valid"],
                    "include_for_esm": False,
                    "sequence_issue": sequence_record["sequence_issue"],
                })
                if include:
                    sequence_records[h]["include_for_esm"] = True
                    linked_sequences.add(h)
                    valid_candidates += 1
            mapping_rows.append({
                "vcc_gene_name": vcc_gene, "gene_id": transcript.gene_id, "gene_id_stable": transcript.gene_id_stable,
                "genomic_region": transcript.seqname,
                "gene_name_gencode": transcript.gene_name, "gene_type": transcript.gene_type,
                "transcript_id": transcript.transcript_id, "transcript_id_stable": transcript.transcript_id_stable,
                "transcript_name": transcript.transcript_name, "transcript_type": transcript.transcript_type,
                "protein_id": transcript.protein_id, "protein_id_stable": transcript.protein_id_stable,
                "protein_sequence_hash": sequence_record["protein_sequence_hash"] if sequence_record else None,
                "protein_length": sequence_record["protein_length"] if sequence_record else None,
                "is_mane_select": mane_meta.get("is_mane_select", "false") == "true",
                "is_mane_plus_clinical": mane_meta.get("is_mane_plus_clinical", "false") == "true",
                "mane_status": mane_meta.get("mane_status", ""),
                "refseq_transcript_id": mane_meta.get("refseq_transcript_id", ""),
                "refseq_protein_id": mane_meta.get("refseq_protein_id", ""),
                "is_perturbation_target": is_target, "sequence_valid": sequence_valid,
                "include_for_esm": include, "exclusion_reason": reason,
                "mapping_status": "mapped_with_esm_candidate" if include else reason or "mapped_no_eligible_protein",
                "review_flag": bool(reason and reason != "no_translation"),
            })

        if accession:
            accession_stats["mapped" if gene_ids else "unmapped"] += 1
            accession_stats["esm_candidate"] += int(valid_candidates > 0)
            accession_stats["protein_coding"] += int(protein_coding_gene)
        if len(gene_ids) > 1:
            regions = {t.seqname for t in all_transcripts}
            alternate = any(re.match(r"^(KI|GL|JH|HSCHR|NW_|NT_)", region) for region in regions)
            classification = "patch_or_alternate_locus" if alternate else "true_ambiguous_symbol_or_other"
            multi_gene_context_counts[classification] += 1
            status, review_reason = "ambiguous_multiple_gene_ids", "multiple_gencode_gene_ids"
        elif not gene_ids:
            status, review_reason = "gene_symbol_not_found", "no_exact_gencode_gene_name_match"
        elif valid_candidates:
            status, review_reason = "mapped_with_esm_candidate", ""
        elif not protein_coding_gene:
            status, review_reason = "non_protein_coding", "no_protein_coding_gene_or_transcript"
        elif special_or_review:
            status, review_reason = "requires_policy_review", "translated_special_transcript_biotype"
        elif invalid_or_missing and translated_count:
            status, review_reason = "translated_sequence_invalid", "translated_record_invalid_or_unusable"
        else:
            status, review_reason = "mapped_no_eligible_protein", "protein_coding_gene_without_eligible_translation"
        status_counts[status] += 1
        summary_rows.append({
            "vcc_gene_name": vcc_gene, "is_perturbation_target": is_target, "mapping_status": status,
            "num_gencode_gene_ids": len(gene_ids), "num_transcripts": len(all_transcripts),
            "num_translated_transcripts": translated_count, "num_unique_protein_sequences": len(linked_sequences),
            "num_esm_candidate_sequences": valid_candidates, "protein_coding_gene": protein_coding_gene,
            "mane_select_available": mane_select, "multi_gene_id_flag": len(gene_ids) > 1,
            "multiple_isoforms": len(linked_sequences) > 1, "esm_candidate_available": valid_candidates > 0,
            "review_flag": status in {"ambiguous_multiple_gene_ids", "requires_policy_review", "translated_sequence_invalid"},
            "review_reason": review_reason,
            "multi_gene_id_classification": (classification if len(gene_ids) > 1 else ""),
        })

    mapping_df = pd.DataFrame(mapping_rows)
    sequences_df = pd.DataFrame(sorted(sequence_records.values(), key=lambda row: row["protein_sequence_hash"]))
    summary_df = pd.DataFrame(summary_rows)
    if len(summary_df) != 18533 or summary_df["vcc_gene_name"].duplicated().any():
        raise ValueError("Gene summary failed the exact-one-row-per-VCC-gene invariant")
    if sequences_df["protein_sequence_hash"].duplicated().any():
        raise ValueError("Protein sequence hash is not unique")
    mapping_df.to_parquet(output_dir / "gene_protein_mapping.parquet", index=False)
    sequences_df.to_parquet(output_dir / "protein_sequences.parquet", index=False)
    summary_df.to_csv(output_dir / "gene_mapping_summary.csv", index=False)

    qc = {
        "vcc_total_genes": len(vcc), "vcc_perturbation_targets": len(perturbations),
        "mapped_to_gencode_gene": int((summary_df["num_gencode_gene_ids"] > 0).sum()),
        "unmapped": int((summary_df["mapping_status"] == "gene_symbol_not_found").sum()),
        "protein_coding_genes": int(summary_df["protein_coding_gene"].sum()),
        "genes_with_translated_proteins": int((summary_df["num_translated_transcripts"] > 0).sum()),
        "genes_with_valid_esm_candidate": int(summary_df["esm_candidate_available"].sum()),
        "genes_with_one_unique_sequence": int((summary_df["num_unique_protein_sequences"] == 1).sum()),
        "genes_with_multiple_unique_sequences": int((summary_df["num_unique_protein_sequences"] > 1).sum()),
        "total_transcripts": len(mapping_df), "translated_transcript_records": int(mapping_df["protein_id"].notna().sum()),
        "unique_protein_sequences": len(sequences_df),
        "transcript_protein_records_collapsed_by_sequence": int(mapping_df["protein_sequence_hash"].notna().sum() - len(sequences_df)),
        "mane_select_genes": int(summary_df["mane_select_available"].sum()),
        "mane_reconciliation_failures": int(sum(row["is_mane_select"] and not row["protein_sequence_hash"] for row in mapping_rows)),
        "perturbation_targets_with_esm_candidate": int(summary_df.loc[summary_df["is_perturbation_target"], "esm_candidate_available"].sum()),
        "unresolved_perturbation_targets": summary_df.loc[(summary_df["is_perturbation_target"]) & (~summary_df["esm_candidate_available"]), "vcc_gene_name"].tolist(),
        "genes_with_multiple_gencode_ids": int(summary_df["multi_gene_id_flag"].sum()), "accession_like": dict(accession_stats),
        "unusual_fasta_characters": dict(unusual_aa), "counts_by_transcript_type": dict(transcript_type_counts),
        "counts_by_vcc_transcript_type": dict(vcc_transcript_type_counts),
        "multi_gene_id_context_classification": dict(multi_gene_context_counts),
        "counts_by_mapping_status": dict(status_counts),
        "counts_by_exclusion_reason": dict(Counter(row["exclusion_reason"] for row in mapping_rows if row["exclusion_reason"])),
        "inputs": {"vcc_gene_names": str(vcc_path), "vcc_pert_counts": str(perturbation_path), "gencode_gtf": str(gtf_path),
                   "gencode_fasta": str(fasta_path), "mane_summary": str(mane_path), "sequence_hash_algorithm": "SHA-256",
                   "gencode_gtf_sha256": file_hash(gtf_path), "gencode_fasta_sha256": file_hash(fasta_path),
                   "mane_summary_sha256": file_hash(mane_path),
                   "mapping_generation_note": "GENCODE + MANE only; no external fallback", "fasta_records_parsed": fasta_record_count},
    }
    (output_dir / "mapping_qc.json").write_text(json.dumps(qc, indent=2, sort_keys=True), encoding="utf-8")
    return qc
