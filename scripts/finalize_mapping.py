"""Create the final, left-joined mapping freeze without overwriting prior snapshots."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data/processed/mapping_refined"
INITIAL = ROOT / "data/processed/mapping"
OUT = ROOT / "data/processed/mapping_final"
GTF_FASTA = ROOT / "data/raw/gencode/v50/gencode.v50.pc_translations.fa.gz"
OUT.mkdir(parents=True, exist_ok=True)


def sha(s: str) -> str:
    return hashlib.sha256(s.encode("ascii")).hexdigest()


def stream_needed_fasta(path: Path, ids: set[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    header = None
    parts: list[str] = []
    def finish() -> None:
        if not header:
            return
        first = header[1:].split("|", 1)[0]
        if first in ids:
            found[first] = "".join(parts).strip().upper()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                finish(); header = line.rstrip("\n"); parts = []
            elif header:
                parts.append(line.strip())
        finish()
    return found


def main() -> None:
    m = pd.read_parquet(REF / "gene_protein_mapping_refined.parquet")
    initial_seq = pd.read_parquet(INITIAL / "protein_sequences.parquet")
    sequences: dict[str, str] = initial_seq.set_index("protein_sequence_hash")["amino_acid_sequence"].to_dict()
    existing = set(sequences)
    missing_rows = m[(m.eligibility == "primary") & (~m.protein_sequence_hash.isin(existing))]
    needed_ids = set(missing_rows.protein_id.dropna()) | set(missing_rows.protein_id_stable.dropna())
    fetched = stream_needed_fasta(GTF_FASTA, needed_ids)
    added = {}
    for _, row in missing_rows.iterrows():
        sequence = fetched.get(row.protein_id) or fetched.get(row.protein_id_stable)
        if sequence is None:
            continue
        h = sha(sequence)
        if h != row.protein_sequence_hash:
            raise RuntimeError(f"sequence hash mismatch for {row.protein_id}: {h} != {row.protein_sequence_hash}")
        sequences[h] = sequence
        added[h] = sequence
    if set(missing_rows.protein_sequence_hash) - set(sequences):
        raise RuntimeError("not all primary hashes could be recovered from the frozen GENCODE FASTA")

    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    selected_hashes = set(m.loc[m.eligibility == "primary", "protein_sequence_hash"].dropna())
    final_seq = pd.DataFrame({"sequence_hash": sorted(selected_hashes)})
    final_seq["amino_acid_sequence"] = final_seq.sequence_hash.map(sequences)
    final_seq["protein_length"] = final_seq.amino_acid_sequence.str.len()
    final_seq["sequence_source"] = "GENCODE 50 protein-coding translations"
    final_seq["annotation_release"] = "GENCODE 50 / Ensembl 116 / GRCh38.p14"
    final_seq["sequence_valid"] = final_seq.amino_acid_sequence.map(lambda x: not (set(x) - allowed) and "*" not in x)
    final_seq["sequence_hash_matches_string"] = final_seq.apply(lambda x: sha(x.amino_acid_sequence) == x.sequence_hash, axis=1)
    final_seq.to_parquet(OUT / "esm_sequences_final.parquet", index=False)

    final_hashes = set(final_seq.sequence_hash)
    mf = m.copy()
    mf["mapping_valid"] = mf.gene_id_stable.notna() & mf.transcript_id_stable.notna() & mf.protein_id_stable.notna() & mf.protein_sequence_hash.notna()
    mf["locus_valid"] = mf.is_primary_locus | mf.is_patch_exception
    mf["transcript_biotype_eligible"] = mf.transcript_type.eq("protein_coding")
    mf["sequence_hash_matches_final_string"] = mf.protein_sequence_hash.map(lambda h: h in final_hashes)
    mf["selected_for_primary_esm"] = mf.eligibility.eq("primary") & mf.sequence_hash_matches_final_string
    mf["final_selection_reason"] = mf.apply(lambda x: "selected_primary_reference_or_MANe_patch" if x.selected_for_primary_esm else ("sequence_table_generation_repaired" if x.eligibility == "primary" and x.sequence_hash_matches_final_string else x.final_esm_reason), axis=1)
    mf.to_parquet(OUT / "gene_protein_mapping_final.parquet", index=False)

    rel = mf[mf.selected_for_primary_esm].copy()
    rel = rel.rename(columns={"protein_sequence_hash": "sequence_hash", "gene_id": "gene_id", "transcript_id": "transcript_id", "protein_id": "protein_id"})
    rel[["sequence_hash", "vcc_gene_name", "resolved_symbol", "gene_id", "transcript_id", "protein_id", "mane_status", "locus_class", "transcript_type", "selected_for_primary_esm"]].drop_duplicates().to_parquet(OUT / "esm_sequence_relationships_final.parquet", index=False)

    complete = pd.read_csv(ROOT / "data/raw/hgnc/current/hgnc_complete_set.txt", sep="\t", dtype=str, keep_default_na=False)
    hgnc = complete.set_index("hgnc_id").to_dict("index")
    resolution = pd.read_csv(REF / "gene_resolution.csv", dtype=str).set_index("vcc_gene_name")
    vcc = pd.read_csv(ROOT / "data/gene_names.csv", dtype=str).iloc[:, 0].drop_duplicates().tolist()
    target_set = set(pd.read_csv(ROOT / "data/pert_counts.csv", dtype=str).target_gene)
    rows = []
    for token in vcc:
        d = mf[mf.vcc_gene_name == token]
        r = resolution.loc[token].to_dict()
        h = hgnc.get(r.get("hgnc_id", ""), {})
        gene_types = set(d.gene_type.dropna())
        if d.empty:
            method = r.get("method", "")
            if method == "ambiguous_hgnc_token": state = "ambiguous_but_explained"
            elif method == "unresolved": state = "unresolved_but_explained"
            elif "pseudogene" in h.get("locus_type", "").lower() or str(h.get("symbol", "")).endswith("P"): state = "resolved_pseudogene"
            elif "non-coding" in h.get("locus_group", "").lower() or "RNA" in h.get("locus_type", ""): state = "resolved_non_protein_coding"
            else: state = "resolved_no_valid_translation"
        elif any(x == "protein_coding" for x in gene_types): state = "resolved_protein_coding"
        elif any("pseudogene" in x for x in gene_types): state = "resolved_pseudogene"
        else: state = "resolved_non_protein_coding"
        if (d.selected_for_primary_esm).any(): reason = "primary ESM sequence available"
        elif state == "ambiguous_but_explained": reason = "multiple exact HGNC candidates; no deterministic choice"
        elif state == "unresolved_but_explained": reason = "no exact approved/alias/history/stable-ID evidence"
        elif state == "resolved_pseudogene": reason = "resolved pseudogene; no primary protein representation"
        elif state == "resolved_non_protein_coding": reason = "resolved non-protein-coding locus"
        else: reason = "resolved but no valid primary translation"
        rows.append({"vcc_gene_name": token, "resolved_gene_symbol": r.get("resolved_symbol", ""), "HGNC_ID": r.get("hgnc_id", ""), "resolution_method": r.get("method", ""), "is_perturbation_target": token in target_set, "gene_resolution_status": state, "protein_coding_status": "protein_coding" if "protein_coding" in gene_types else ("pseudogene" if any("pseudogene" in x for x in gene_types) else "non_protein_coding_or_unknown"), "num_primary_gene_loci": d.loc[d.is_primary_locus, "gene_id_stable"].nunique(), "num_primary_transcripts": d.loc[d.selected_for_primary_esm, "transcript_id_stable"].nunique(), "num_primary_unique_protein_sequences": d.loc[d.selected_for_primary_esm, "protein_sequence_hash"].nunique(), "MANE_available": bool(d.is_mane_select.astype(str).str.lower().eq("true").any()), "primary_esm_available": bool(d.selected_for_primary_esm.any()), "ambiguity_status": "ambiguous_unresolved" if state == "ambiguous_but_explained" else "none", "review_required": state in {"ambiguous_but_explained", "unresolved_but_explained"}, "final_reason": reason})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "gene_mapping_summary_final.csv", index=False)
    resolution.reset_index().to_csv(OUT / "gene_resolution_final.csv", index=False)

    qc = {"before_summary_rows": 18395, "after_summary_rows": int(len(summary)), "after_summary_unique_vcc": int(summary.vcc_gene_name.nunique()), "vcc_total": len(vcc), "missing_summary_rows_repaired": int(len(set(vcc) - set(pd.read_csv(REF / 'gene_mapping_summary_refined.csv', dtype=str).vcc_gene_name))), "missing_primary_hashes_before": int(missing_rows.protein_sequence_hash.nunique()), "missing_primary_hashes_after": int(len(set(mf.loc[mf.selected_for_primary_esm, 'protein_sequence_hash']) - set(final_seq.sequence_hash))), "recovered_sequence_hashes": int(len(added)), "primary_relationships_missing_sequence": int((~rel.sequence_hash.isin(set(final_seq.sequence_hash))).sum()), "orphan_final_sequences": int(len(set(final_seq.sequence_hash) - set(rel.sequence_hash))), "hash_conflicts": int(final_seq.groupby('sequence_hash').amino_acid_sequence.nunique().gt(1).sum()), "unexpected_primary_loci": mf.loc[mf.selected_for_primary_esm & ~mf.locus_class.eq('primary_reference') & ~mf.is_patch_exception, 'locus_class'].value_counts().to_dict(), "primary_stop_sequences": int(final_seq.amino_acid_sequence.str.contains('*', regex=False).sum()), "sequence_hash_policy": "SHA-256 of exact uppercase amino-acid string passed to future ESM; no replacement/truncation", "policy_changed": False, "esm2_started": False}
    (OUT / "mapping_final_qc.json").write_text(json.dumps(qc, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
