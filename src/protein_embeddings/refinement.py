"""Deterministic refinement of the frozen VCC-to-GENCODE mapping."""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from .mapping import parse_attributes, parse_mane, sequence_hash, stable_id


PRIMARY_RE = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$", re.I)


def locus_class(seqname: str) -> str:
    s = str(seqname or "")
    u = s.upper()
    if PRIMARY_RE.match(s):
        return "primary_reference"
    if "PATCH" in u or u.startswith("CHR_") and "PATCH" in u:
        return "patch"
    if u.startswith(("KI", "GL", "JH", "HSCHR", "ALT_REF")):
        return "alternate_locus"
    if u.startswith(("HLA-", "HG", "NA_")) or "HAP" in u:
        return "haplotype"
    if u.startswith(("NW_", "NT_", "NG_", "NZ_")):
        return "scaffold"
    return "other"


def _split(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [x.strip() for x in str(value).split("|") if x.strip()]


def load_hgnc(complete: Path, withdrawn: Path) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, str]]]]:
    h = pd.read_csv(complete, sep="\t", dtype=str, keep_default_na=False)
    rows = h.to_dict("records")
    by_token: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        for col, method in (("symbol", "approved_symbol"), ("prev_symbol", "previous_symbol"), ("alias_symbol", "alias_symbol")):
            for token in _split(row.get(col, "")):
                by_token[token].append({"hgnc_id": row.get("hgnc_id", ""), "symbol": row.get("symbol", ""), "method": method, "status": row.get("status", "")})
    w = pd.read_csv(withdrawn, sep="\t", dtype=str, keep_default_na=False).to_dict("records")
    for row in w:
        token = row.get("WITHDRAWN_SYMBOL", "").strip()
        if token:
            by_token[token].append({"hgnc_id": row.get("HGNC_ID", ""), "symbol": "", "method": "withdrawn_symbol", "status": "Entry Withdrawn", "merged_into": row.get("MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)", "")})
    stable: dict[str, dict[str, Any]] = {}
    for row in rows:
        ens = stable_id(row.get("ensembl_gene_id", ""))
        if ens:
            stable[ens] = row
    return stable, by_token


def resolve_symbol(token: str, stable: dict[str, dict[str, Any]], by_token: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    if token in stable:
        row = stable[token]
        return {"resolved_symbol": row.get("symbol", ""), "hgnc_id": row.get("hgnc_id", ""), "method": "stable_ensembl_id", "evidence": token, "status": row.get("status", "")}
    candidates = by_token.get(token, [])
    approved = [x for x in candidates if x["method"] == "approved_symbol"]
    if approved:
        candidates = approved
    if len({x["hgnc_id"] for x in candidates}) == 1 and candidates:
        x = candidates[0]
        if x["method"] == "withdrawn_symbol":
            merged = x.get("merged_into", "")
            parts = merged.split("|") if merged else []
            return {"resolved_symbol": parts[1] if len(parts) > 1 else "", "hgnc_id": parts[0] if parts else x["hgnc_id"], "method": "withdrawn_or_merged_symbol", "evidence": token, "status": x["status"]}
        return {"resolved_symbol": x["symbol"], "hgnc_id": x["hgnc_id"], "method": x["method"], "evidence": token, "status": x["status"]}
    if len({x["hgnc_id"] for x in candidates}) > 1:
        return {"resolved_symbol": "", "hgnc_id": "", "method": "ambiguous_hgnc_token", "evidence": token, "status": ""}
    return {"resolved_symbol": token, "hgnc_id": "", "method": "unresolved", "evidence": "", "status": ""}


def star_class(seq: str) -> str:
    n = seq.count("*")
    if not n:
        return "none"
    if n > 1:
        return "multiple_stops"
    return "terminal_stop" if seq.endswith("*") else "internal_stop"


def targeted_records(gtf: Path, fasta: Path, mane: Path, symbols: set[str]) -> list[dict[str, Any]]:
    """Stream only rescued symbols, avoiding a second multi-gigabyte in-memory copy."""
    wanted: list[dict[str, str]] = []
    with gzip.open(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            a = parse_attributes(fields[8])
            if a.get("gene_name") in symbols and a.get("protein_id"):
                wanted.append({"seqname": fields[0], "gene_id": a.get("gene_id", ""), "gene_name": a.get("gene_name", ""), "gene_type": a.get("gene_type", a.get("gene_biotype", "")), "transcript_id": a.get("transcript_id", ""), "transcript_name": a.get("transcript_name", ""), "transcript_type": a.get("transcript_type", a.get("transcript_biotype", "")), "protein_id": a.get("protein_id", "")})
    protein_ids = {x["protein_id"] for x in wanted}
    proteins: dict[str, dict[str, Any]] = {}
    with gzip.open(fasta, "rt", encoding="utf-8") as handle:
        header = None
        parts: list[str] = []
        def finish() -> None:
            if not header:
                return
            fields = header[1:].strip().split("|")
            if len(fields) >= 8 and fields[0] in protein_ids:
                seq = "".join(parts).strip().upper()
                proteins[fields[0]] = {"sequence": seq, "length": len(seq), "valid": not (set(seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO")), "issue": "invalid_character" if set(seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO") else ""}
        for line in handle:
            if line.startswith(">"):
                finish(); header = line.rstrip("\n"); parts = []
            elif header:
                parts.append(line.strip())
        finish()
    mane_rows = parse_mane(mane)
    out = []
    for x in wanted:
        sr = proteins.get(x["protein_id"])
        if not sr:
            continue
        tid = x["transcript_id"]
        mh = mane_rows.get(tid) or mane_rows.get(stable_id(tid)) or {}
        out.append({**x, "gene_id_stable": stable_id(x["gene_id"]), "transcript_id_stable": stable_id(tid), "protein_id_stable": stable_id(x["protein_id"]), "protein_sequence_hash": sequence_hash(sr["sequence"]), "protein_length": sr["length"], "is_mane_select": mh.get("MANE_status", "") == "MANE Select", "is_mane_plus_clinical": mh.get("MANE_status", "") == "MANE Plus Clinical", "mane_status": mh.get("MANE_status", ""), "refseq_transcript_id": mh.get("RefSeq_nuc", ""), "refseq_protein_id": mh.get("RefSeq_prot", ""), "is_perturbation_target": False, "sequence_valid": sr["valid"], "include_for_esm": sr["valid"], "exclusion_reason": sr.get("issue", ""), "mapping_status": "hgnc_rescued", "review_flag": True})
    return out


def refine(initial_dir: Path, out_dir: Path, vcc_path: Path, complete: Path, withdrawn: Path, gtf: Path, fasta: Path, mane: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_parquet(initial_dir / "gene_protein_mapping.parquet")
    seqs = pd.read_parquet(initial_dir / "protein_sequences.parquet")
    seq_by_hash = seqs.set_index("protein_sequence_hash")["amino_acid_sequence"].to_dict()
    stable, by_token = load_hgnc(complete, withdrawn)
    vcc = pd.read_csv(vcc_path, dtype=str).iloc[:, 0].astype(str).tolist()
    initial_symbols = set(mapping["gene_name_gencode"].dropna())
    resolutions = []
    for token in vcc:
        r = resolve_symbol(token, stable, by_token)
        if r["method"] == "unresolved" and token in initial_symbols:
            r.update(resolved_symbol=token, method="exact_gencode_symbol", evidence="GENCODE gene_name exact", status="")
        resolutions.append({"vcc_gene_name": token, **r})
    resolution = pd.DataFrame(resolutions)
    resolution["alias_rescued"] = resolution.method.isin(["previous_symbol", "alias_symbol", "withdrawn_or_merged_symbol", "stable_ensembl_id"])
    resolution["resolution_status"] = resolution.apply(lambda x: "resolved" if x.resolved_symbol else "unresolved", axis=1)

    # Add GENCODE records for HGNC-rescued symbols not represented in the initial exact-name table.
    rescued_symbols = set(resolution.loc[resolution.alias_rescued, "resolved_symbol"]) - initial_symbols
    if rescued_symbols:
        vcc_by_symbol = {x["resolved_symbol"]: x["vcc_gene_name"] for x in resolutions if x["resolved_symbol"] in rescued_symbols}
        extra = [{"vcc_gene_name": vcc_by_symbol[x["gene_name"]], "gene_id": x["gene_id"], "gene_id_stable": x["gene_id_stable"], "genomic_region": x["seqname"], "gene_name_gencode": x["gene_name"], "gene_type": x["gene_type"], "transcript_id": x["transcript_id"], "transcript_id_stable": x["transcript_id_stable"], "transcript_name": x["transcript_name"], "transcript_type": x["transcript_type"], "protein_id": x["protein_id"], "protein_id_stable": x["protein_id_stable"], "protein_sequence_hash": x["protein_sequence_hash"], "protein_length": x["protein_length"], "is_mane_select": x["is_mane_select"], "is_mane_plus_clinical": x["is_mane_plus_clinical"], "mane_status": x["mane_status"], "refseq_transcript_id": x["refseq_transcript_id"], "refseq_protein_id": x["refseq_protein_id"], "is_perturbation_target": x["is_perturbation_target"], "sequence_valid": x["sequence_valid"], "include_for_esm": x["include_for_esm"], "exclusion_reason": x["exclusion_reason"], "mapping_status": x["mapping_status"], "review_flag": x["review_flag"]} for x in targeted_records(gtf, fasta, mane, rescued_symbols)]
        if extra:
            mapping = pd.concat([mapping, pd.DataFrame(extra)], ignore_index=True)

    res_by_vcc = resolution.set_index("vcc_gene_name").to_dict("index")
    mapping["resolved_symbol"] = mapping.vcc_gene_name.map(lambda x: res_by_vcc.get(x, {}).get("resolved_symbol", ""))
    mapping["locus_class"] = mapping.genomic_region.map(locus_class)
    mapping["star_class"] = mapping.protein_sequence_hash.map(lambda h: star_class(str(seq_by_hash.get(h, ""))))
    mapping["sequence_has_stop"] = mapping.star_class.ne("none")
    mapping["is_primary_locus"] = mapping.locus_class.eq("primary_reference")
    mapping["is_patch_exception"] = mapping.is_mane_select & mapping.locus_class.eq("patch")
    mapping["eligible_locus"] = mapping.is_primary_locus | mapping.is_patch_exception
    mapping["biotype_policy"] = mapping.transcript_type.map(lambda x: "ordinary_protein_coding" if x == "protein_coding" else ("special_translated" if x in {"nonsense_mediated_decay", "non_stop_decay", "protein_coding_LoF", "polymorphic_pseudogene"} else "other_transcript_type"))

    group = mapping.groupby("vcc_gene_name", sort=False)
    selected = {}
    for token, d in group:
        p = sorted(set(d.loc[d.is_primary_locus, "gene_id_stable"]))
        selected[token] = p[0] if len(p) == 1 else (sorted(set(d.loc[d.is_patch_exception, "gene_id_stable"]))[0] if not p and len(set(d.loc[d.is_patch_exception, "gene_id_stable"])) == 1 else "")
    mapping["selected_primary_gene_id"] = mapping.vcc_gene_name.map(selected)
    mapping["selected_for_primary_esm"] = mapping.gene_id_stable.eq(mapping.selected_primary_gene_id) & mapping.selected_primary_gene_id.ne("")
    def elig(row: pd.Series) -> str:
        if row.star_class != "none": return "review"
        if not row.sequence_valid or not row.protein_sequence_hash: return "excluded"
        if row.biotype_policy == "ordinary_protein_coding": return "primary" if row.selected_for_primary_esm else "secondary"
        if row.biotype_policy == "special_translated": return "secondary"
        return "excluded"
    mapping["eligibility"] = mapping.apply(elig, axis=1)
    mapping["final_esm_reason"] = mapping.apply(lambda r: "selected_primary_locus" if r.eligibility == "primary" and r.locus_class == "primary_reference" else ("mane_select_patch_exception" if r.eligibility == "primary" else r.star_class if r.star_class != "none" else r.biotype_policy), axis=1)
    mapping.to_parquet(out_dir / "gene_protein_mapping_refined.parquet", index=False)

    primary = mapping[mapping.eligibility == "primary"].copy()
    primary = primary[primary.protein_sequence_hash.isin(seq_by_hash)]
    esm_seq = pd.DataFrame({"protein_sequence_hash": sorted(primary.protein_sequence_hash.unique())})
    esm_seq["amino_acid_sequence"] = esm_seq.protein_sequence_hash.map(seq_by_hash)
    esm_seq["protein_length"] = esm_seq.amino_acid_sequence.str.len()
    esm_seq["sequence_has_stop"] = esm_seq.amino_acid_sequence.str.contains("*", regex=False)
    esm_seq.to_parquet(out_dir / "esm_sequences.parquet", index=False)
    primary[["vcc_gene_name", "resolved_symbol", "gene_id_stable", "transcript_id_stable", "protein_id_stable", "protein_sequence_hash", "is_mane_select", "is_mane_plus_clinical", "locus_class", "eligibility"]].drop_duplicates().to_parquet(out_dir / "esm_sequence_relationships.parquet", index=False)

    summary_rows = []
    for token, d in mapping.groupby("vcc_gene_name", sort=False):
        rr = res_by_vcc.get(token, {})
        summary_rows.append({"vcc_gene_name": token, "resolved_symbol": rr.get("resolved_symbol", ""), "resolution_method": rr.get("method", ""), "hgnc_id": rr.get("hgnc_id", ""), "resolution_status": rr.get("resolution_status", "unresolved"), "selected_primary_gene_id": selected.get(token, ""), "num_gene_ids": d.gene_id_stable.nunique(), "num_primary_gene_ids": d.loc[d.is_primary_locus, "gene_id_stable"].nunique(), "num_patch_gene_ids": d.loc[d.locus_class == "patch", "gene_id_stable"].nunique(), "num_alt_gene_ids": d.loc[d.locus_class == "alternate_locus", "gene_id_stable"].nunique(), "num_haplotype_gene_ids": d.loc[d.locus_class == "haplotype", "gene_id_stable"].nunique(), "num_scaffold_gene_ids": d.loc[d.locus_class == "scaffold", "gene_id_stable"].nunique(), "primary_esm_available": bool((d.eligibility == "primary").any()), "n_primary_sequences": d.loc[d.eligibility == "primary", "protein_sequence_hash"].nunique(), "n_secondary_records": int((d.eligibility == "secondary").sum()), "n_review_records": int((d.eligibility == "review").sum()), "n_excluded_records": int((d.eligibility == "excluded").sum()), "review_required": bool((d.eligibility == "review").any() or not rr.get("resolved_symbol")), "review_reason": "ambiguous_or_unresolved" if not rr.get("resolved_symbol") else ("multiple_primary_loci" if d.loc[d.is_primary_locus, "gene_id_stable"].nunique() > 1 else "")})
    pd.DataFrame(summary_rows).to_csv(out_dir / "gene_mapping_summary_refined.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    resolution.to_csv(out_dir / "gene_resolution.csv", index=False)
    qc = {"n_vcc": len(vcc), "n_refined_mapping_rows": int(len(mapping)), "n_final_esm_sequences": int(len(esm_seq)), "eligibility_counts": mapping.eligibility.value_counts().to_dict(), "locus_counts": mapping.locus_class.value_counts().to_dict(), "star_counts": mapping.star_class.value_counts().to_dict(), "resolution_methods": resolution.method.value_counts().to_dict(), "alias_rescued_vcc": sorted(resolution.loc[resolution.alias_rescued, "vcc_gene_name"]), "esm2_started": False}
    (out_dir / "refinement_qc.json").write_text(json.dumps(qc, indent=2, sort_keys=True), encoding="utf-8")
