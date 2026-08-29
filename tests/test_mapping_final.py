from pathlib import Path
import json

import pandas as pd


def test_final_mapping_freeze_invariants():
    root = Path(__file__).resolve().parents[1]
    final = root / "data/processed/mapping_final"
    summary = pd.read_csv(final / "gene_mapping_summary_final.csv")
    mapping = pd.read_parquet(final / "gene_protein_mapping_final.parquet")
    sequences = pd.read_parquet(final / "esm_sequences_final.parquet")
    relationships = pd.read_parquet(final / "esm_sequence_relationships_final.parquet")
    hashes = set(sequences.sequence_hash)
    assert len(summary) == 18_533
    assert summary.vcc_gene_name.nunique() == 18_533
    assert summary.is_perturbation_target.sum() == 300
    assert summary.loc[summary.is_perturbation_target, "primary_esm_available"].all()
    assert mapping.loc[mapping.selected_for_primary_esm, "protein_sequence_hash"].isin(hashes).all()
    assert set(sequences.sequence_hash).issubset(set(relationships.sequence_hash))
    assert set(relationships.sequence_hash).issubset(hashes)
    assert sequences.sequence_hash.is_unique
    assert sequences.groupby("sequence_hash").amino_acid_sequence.nunique().max() == 1
    assert not sequences.amino_acid_sequence.str.contains("*", regex=False).any()
    assert set(mapping.loc[mapping.selected_for_primary_esm, "locus_class"]) == {"primary_reference"}


def test_final_qc_records_zero_reconciliation_failures():
    root = Path(__file__).resolve().parents[1]
    qc = json.loads((root / "data/processed/mapping_final/mapping_final_qc.json").read_text())
    assert qc["after_summary_rows"] == 18_533
    assert qc["missing_primary_hashes_after"] == 0
    assert qc["primary_relationships_missing_sequence"] == 0
    assert qc["orphan_final_sequences"] == 0
    assert qc["hash_conflicts"] == 0
    assert qc["primary_stop_sequences"] == 0
