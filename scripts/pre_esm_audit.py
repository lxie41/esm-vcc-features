"""Read-only pre-ESM readiness audit for the frozen refined mapping."""
from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MAP = ROOT / "data/processed/mapping"
REF = ROOT / "data/processed/mapping_refined"
OUT = ROOT / "data/processed/pre_esm_audit"
OUT.mkdir(parents=True, exist_ok=True)


def write(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False)


def qstats(s: pd.Series) -> dict:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if not len(s):
        return {k: None for k in ["count", "min", "median", "mean", "p90", "p95", "p99", "p99_9", "max"]}
    return {"count": int(len(s)), "min": int(s.min()), "median": float(s.median()), "mean": float(s.mean()), "p90": float(s.quantile(.90)), "p95": float(s.quantile(.95)), "p99": float(s.quantile(.99)), "p99_9": float(s.quantile(.999)), "max": int(s.max())}


def main() -> None:
    m = pd.read_parquet(REF / "gene_protein_mapping_refined.parquet")
    seq = pd.read_parquet(REF / "esm_sequences.parquet")
    rel = pd.read_parquet(REF / "esm_sequence_relationships.parquet")
    summary = pd.read_csv(REF / "gene_mapping_summary_refined.csv", dtype=str)
    resolution = pd.read_csv(REF / "gene_resolution.csv", dtype=str)
    initial_seq = pd.read_parquet(MAP / "protein_sequences.parquet")
    initial_map = pd.read_parquet(MAP / "gene_protein_mapping.parquet")
    targets = pd.read_csv(ROOT / "data/pert_counts.csv", dtype=str)["target_gene"].drop_duplicates().tolist()
    seq_lookup = seq.set_index("protein_sequence_hash")["amino_acid_sequence"].to_dict()
    all_seq_lookup = initial_seq.set_index("protein_sequence_hash")["amino_acid_sequence"].to_dict()
    m["sequence"] = m["protein_sequence_hash"].map(all_seq_lookup)
    m["sequence_length_actual"] = m["sequence"].str.len()
    primary = m[m["eligibility"] == "primary"].copy()
    primary_hashes = set(seq["protein_sequence_hash"])
    primary_mapping_missing_final = m[(m.eligibility == "primary") & (~m.protein_sequence_hash.isin(primary_hashes))].copy()
    write(primary_mapping_missing_final[["vcc_gene_name", "resolved_symbol", "gene_id_stable", "transcript_id_stable", "protein_id_stable", "protein_sequence_hash", "locus_class", "eligibility"]], "primary_mapping_not_in_esm.csv")
    summary_missing = resolution[~resolution.vcc_gene_name.isin(set(summary.vcc_gene_name))].copy()
    summary_missing["has_refined_mapping_row"] = summary_missing.vcc_gene_name.isin(set(m.vcc_gene_name))
    summary_missing["is_perturbation_target"] = summary_missing.vcc_gene_name.isin(targets)
    summary_missing["reason"] = summary_missing.apply(lambda x: "no refined mapping rows" if not x.has_refined_mapping_row else "summary construction omitted mapped VCC", axis=1)
    write(summary_missing, "summary_coverage_gaps.csv")

    # Core integrity.
    rel_hashes = set(rel["protein_sequence_hash"])
    core = {
        "summary_rows": int(len(summary)), "summary_unique_vcc": int(summary.vcc_gene_name.nunique()),
        "mapping_rows": int(len(m)), "sequence_rows": int(len(seq)), "relationship_rows": int(len(rel)),
        "sequence_hash_unique": bool(seq.protein_sequence_hash.is_unique),
        "sequence_hash_to_one_string": bool(seq.groupby("protein_sequence_hash").amino_acid_sequence.nunique().max() <= 1),
        "all_sequences_have_relationship": bool(primary_hashes <= rel_hashes),
        "all_relationships_have_sequence": bool(rel_hashes <= primary_hashes),
        "all_primary_mapping_hashes_have_final_sequence": bool(len(primary_mapping_missing_final) == 0),
        "primary_mapping_records_missing_final_sequence": int(len(primary_mapping_missing_final)),
        "all_primary_rows_resolved_vcc": bool(primary.resolved_symbol.notna().all() and primary.resolved_symbol.ne("").all()),
        "sequence_lengths_agree": bool((seq.protein_length == seq.amino_acid_sequence.str.len()).all()),
        "orphan_primary_hashes": sorted(primary_hashes - rel_hashes),
        "orphan_relationship_hashes": sorted(rel_hashes - primary_hashes),
    }

    # Target audit, always one row per target.
    sm = summary.set_index("vcc_gene_name")
    rr = resolution.set_index("vcc_gene_name")
    rows = []
    for t in targets:
        d = m[m.vcc_gene_name == t]
        r = rr.loc[t] if t in rr.index else pd.Series(dtype=str)
        s = sm.loc[t] if t in sm.index else pd.Series(dtype=str)
        rows.append({"vcc_gene_name": t, "resolved_gene_symbol": r.get("resolved_symbol", ""), "hgnc_resolution_method": r.get("method", ""), "gencode_gene_ids": "|".join(sorted(d.gene_id_stable.dropna().unique())), "number_primary_esm_sequences": int(d.loc[d.eligibility == "primary", "protein_sequence_hash"].nunique()), "number_secondary_review_sequences": int(d.loc[d.eligibility.isin(["secondary", "review"]), "protein_sequence_hash"].nunique()), "number_protein_sequences": int(d.protein_sequence_hash.dropna().nunique()), "mane_select_available": bool(d.is_mane_select.astype(str).str.lower().eq("true").any()), "mapping_status": s.get("resolution_status", r.get("resolution_status", "unresolved")), "review_flag": bool(d.eligibility.isin(["review"]).any() or not r.get("resolved_symbol", "")), "perturbation_exception": ""})
    target_df = pd.DataFrame(rows)
    target_df.loc[target_df.vcc_gene_name == "TMEM104", "perturbation_exception"] = "HGNC previous_symbol: SLC38A12"
    write(target_df, "perturbation_target_audit.csv")

    # Unresolved and ambiguous HGNC evidence.
    complete = pd.read_csv(ROOT / "data/raw/hgnc/current/hgnc_complete_set.txt", sep="\t", dtype=str, keep_default_na=False)
    token_cols = ["symbol", "prev_symbol", "alias_symbol"]
    unresolved = resolution[(resolution.resolved_symbol.fillna("").eq("")) | (~resolution.vcc_gene_name.isin(set(m.vcc_gene_name)))].copy()
    unresolved["is_perturbation_target"] = unresolved.vcc_gene_name.isin(targets)
    def cls(x: str) -> str:
        if re.match(r"^[A-Z]{2}\d+\.\d+$", x) or re.match(r"^H3\.[A-Z]$", x): return "accession_like_or_versioned"
        if x.startswith("LOC") and x[3:].isdigit(): return "LOC_like"
        if x.startswith(("LINC", "MIR", "SNOR", "RNU")): return "noncoding_or_RNA_like"
        if x.endswith("P") or "PSEUDO" in x.upper(): return "pseudogene_like"
        return "standard_looking_or_unknown_symbol"
    unresolved["identifier_class"] = unresolved.vcc_gene_name.map(cls)
    unresolved["likely_protein_coding_relevance"] = unresolved.identifier_class.map(lambda x: "low" if x in {"noncoding_or_RNA_like", "pseudogene_like"} else "unknown")
    unresolved["recommended_future_evidence_source"] = unresolved.identifier_class.map(lambda x: "targeted NCBI/RefSeq or accession record review" if x == "accession_like_or_versioned" else "manual HGNC/GENCODE review")
    unresolved["reason_unresolved"] = unresolved.apply(lambda x: "no refined mapping rows; no ESM-backed GENCODE record" if x.vcc_gene_name in set(summary_missing.vcc_gene_name) else ("no exact HGNC or stable-ID evidence" if x.method == "unresolved" else x.method), axis=1)
    write(unresolved, "unresolved_gene_audit.csv")
    ambiguous = resolution[resolution.method.eq("ambiguous_hgnc_token")].copy()
    amb_rows = []
    for _, x in ambiguous.iterrows():
        hits = []
        for _, g in complete.iterrows():
            if any(x.vcc_gene_name in str(g.get(c, "")).split("|") for c in token_cols):
                hits.append(f"{g.get('hgnc_id','')}:{g.get('symbol','')}:{g.get('ensembl_gene_id','')}:{g.get('locus_type','')}:{g.get('status','')}")
        amb_rows.append({"vcc_gene_name": x.vcc_gene_name, "all_hgnc_matches": "|".join(hits), "candidate_count": len(hits), "selected": "", "why_ambiguous": "token has multiple exact HGNC candidates", "perturbation_target": x.vcc_gene_name in targets})
    write(pd.DataFrame(amb_rows), "ambiguous_hgnc_audit.csv")

    # Locus and leakage audit.
    locus_counts = primary.groupby("locus_class").protein_sequence_hash.nunique().rename("unique_primary_sequences").reset_index()
    locus_counts["vcc_genes"] = locus_counts.locus_class.map(primary.groupby("locus_class").vcc_gene_name.nunique())
    write(locus_counts, "locus_audit.csv")
    nonprimary = primary[~primary.locus_class.eq("primary_reference")].copy()
    write(nonprimary[["vcc_gene_name", "resolved_symbol", "gene_id_stable", "transcript_id_stable", "protein_id_stable", "protein_sequence_hash", "locus_class", "is_mane_select", "final_esm_reason"]], "non_primary_locus_exceptions.csv")

    # Biotypes and stop characters.
    b = m[m.protein_sequence_hash.notna()].copy()
    ba = b.groupby("transcript_type").agg(transcript_records=("transcript_id_stable", "size"), unique_protein_sequences=("protein_sequence_hash", "nunique"), genes_affected=("vcc_gene_name", "nunique")).reset_index()
    elig = pd.crosstab(b.transcript_type, b.eligibility).reset_index()
    biotype = ba.merge(elig, on="transcript_type", how="left").fillna(0)
    write(biotype, "biotype_audit.csv")
    stop = m[m.sequence.fillna("").str.contains("*", regex=False)].copy()
    stop["stop_positions"] = stop.sequence.map(lambda x: ",".join(str(i + 1) for i, c in enumerate(x) if c == "*"))
    stop_summary = stop.groupby("star_class").agg(sequence_records=("protein_sequence_hash", "size"), unique_sequences=("protein_sequence_hash", "nunique"), genes=("vcc_gene_name", "nunique")).reset_index()
    write(stop_summary, "stop_character_audit.csv")
    write(stop[stop.star_class.eq("multiple_stops")][["vcc_gene_name", "transcript_id_stable", "protein_id_stable", "sequence_length_actual", "stop_positions", "transcript_type", "locus_class", "mane_status", "eligibility"]].drop_duplicates(), "multiple_stop_details.csv")

    # Alphabet and lengths.
    chars = sorted(set("".join(seq.amino_acid_sequence.astype(str))))
    char_rows = []
    for c in chars:
        hit = m[m.sequence.fillna("").str.contains(re.escape(c), regex=True)]
        char_rows.append({"character": c, "primary_sequences": int(seq[seq.amino_acid_sequence.str.contains(re.escape(c), regex=True)].protein_sequence_hash.nunique()), "all_refined_records": int(hit.protein_sequence_hash.nunique()), "primary_records": int((hit.eligibility == "primary").sum()), "secondary_records": int((hit.eligibility == "secondary").sum()), "review_records": int((hit.eligibility == "review").sum()), "excluded_records": int((hit.eligibility == "excluded").sum()), "is_standard_20aa": c in set("ACDEFGHIKLMNPQRSTVWY")})
    write(pd.DataFrame(char_rows), "sequence_character_audit.csv")
    length = seq.copy(); length["length_bin"] = pd.cut(length.protein_length, [-1, 512, 1022, 1500, 2000, 3000, 5000, 10**9], labels=["<=512", "513-1022", "1023-1500", "1501-2000", "2001-3000", "3001-5000", ">5000"])
    length_audit = length.groupby("length_bin", observed=False).agg(sequence_count=("protein_sequence_hash", "nunique"), gene_count=("protein_sequence_hash", lambda x: primary[primary.protein_sequence_hash.isin(set(x))].vcc_gene_name.nunique())).reset_index()
    length_audit["percentage"] = length_audit.sequence_count / len(seq) * 100
    write(length_audit, "protein_length_audit.csv")

    # Gene isoform distributions and outliers.
    iso = primary.groupby(["vcc_gene_name", "resolved_symbol"], dropna=False).agg(unique_primary_sequences=("protein_sequence_hash", "nunique"), transcripts=("transcript_id_stable", "nunique"), transcript_types=("transcript_type", lambda x: "|".join(sorted(set(x)))), locus_classes=("locus_class", lambda x: "|".join(sorted(set(x))))).reset_index()
    iso_counts = iso.set_index("vcc_gene_name").unique_primary_sequences.reindex(targets if False else summary.vcc_gene_name).fillna(0)
    dist = {"genes_with_0": int((iso_counts == 0).sum()), "genes_with_1": int((iso_counts == 1).sum()), "genes_with_gt1": int((iso_counts > 1).sum()), **qstats(iso_counts)}
    write(iso.sort_values(["unique_primary_sequences", "transcripts"], ascending=False).head(25), "isoform_count_outliers.csv")
    write(iso, "gene_isoform_counts.csv")
    long = seq[seq.protein_length > 1022].copy(); long["genes"] = long.protein_sequence_hash.map(lambda h: "|".join(sorted(primary.loc[primary.protein_sequence_hash == h, "vcc_gene_name"].unique())))
    write(long[["protein_sequence_hash", "protein_length", "genes"]], "long_protein_genes.csv")

    # Cross-gene sequence relationships.
    cross = primary.groupby("protein_sequence_hash").agg(gene_count=("vcc_gene_name", "nunique"), genes=("vcc_gene_name", lambda x: "|".join(sorted(set(x)))), sequence_length=("protein_length", "first")).reset_index()
    write(cross[cross.gene_count > 1], "cross_gene_sequence_hashes.csv")

    total_residues = int(seq.protein_length.sum())
    storage = {"sequences": len(seq), "total_residues": total_residues, "mean_fp32_bytes": len(seq)*1280*4, "mean_fp16_bytes": len(seq)*1280*2, "mean_std_fp32_bytes": len(seq)*2560*4, "mean_std_fp16_bytes": len(seq)*2560*2, "raw_residue_fp32_bytes": total_residues*1280*4, "raw_residue_fp16_bytes": total_residues*1280*2}
    for k, v in list(storage.items()):
        if k.endswith("bytes"): storage[k.replace("bytes", "GiB")] = v / 1024**3

    target_counts = target_df.set_index("vcc_gene_name")
    core.update({"resolved_vcc": int(resolution.resolved_symbol.fillna("").ne("").sum()), "unresolved_vcc": int(resolution.resolved_symbol.fillna("").eq("").sum()), "ambiguous_hgnc": int((resolution.method == "ambiguous_hgnc_token").sum()), "primary_sequences": len(seq), "primary_reference_sequences": int(primary[primary.locus_class == "primary_reference"].protein_sequence_hash.nunique()), "mane_patch_sequences": int(primary[primary.locus_class == "patch"].protein_sequence_hash.nunique()), "alternate_sequences": int(primary[primary.locus_class == "alternate_locus"].protein_sequence_hash.nunique()), "haplotype_sequences": int(primary[primary.locus_class == "haplotype"].protein_sequence_hash.nunique()), "scaffold_sequences": int(primary[primary.locus_class == "scaffold"].protein_sequence_hash.nunique()), "ordinary_patch_sequences": int(primary[primary.locus_class == "patch"].loc[~primary[primary.locus_class == "patch"].is_mane_select].protein_sequence_hash.nunique()), "primary_stop_sequences": int(seq[seq.amino_acid_sequence.str.contains("*", regex=False)].protein_sequence_hash.nunique()), "target_total": len(targets), "target_with_primary": int((target_counts.number_primary_esm_sequences.astype(int) > 0).sum()), "target_review": int(target_counts.review_flag.astype(str).str.lower().eq("true").sum()), "target_unresolved": int(target_counts.resolved_gene_symbol.fillna("").eq("").sum()), "length_stats": qstats(seq.protein_length), "length_over_1022": int((seq.protein_length > 1022).sum()), "length_over_2000": int((seq.protein_length > 2000).sum()), "length_over_5000": int((seq.protein_length > 5000).sum()), "long_target_genes": int(target_df[target_df.vcc_gene_name.isin(long.genes.str.split("|", expand=True).stack().unique())].vcc_gene_name.nunique()) if len(long) else 0, "isoform_distribution": dist, "storage_estimates": storage, "initial_unique_translated_sequences": int(initial_seq.protein_sequence_hash.nunique()), "initial_mapping_rows": int(len(initial_map)), "refined_mapping_rows": int(len(m)), "sequence_hashes_removed_from_initial": int(len(set(initial_seq.protein_sequence_hash) - set(seq.protein_sequence_hash))), "sequence_hash_change_without_normalization": True, "mapping_policy_modified": False, "esm2_downloaded": False, "esm2_run": False, "embedding_generated": False, "label_independent": True, "hidden_data_leakage_found": False})
    (OUT / "pre_esm_audit_summary.json").write_text(json.dumps(core, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
