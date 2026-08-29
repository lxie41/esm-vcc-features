from pathlib import Path

import pandas as pd

from protein_embeddings.refinement import locus_class, star_class


def test_locus_classes_are_deterministic():
    assert locus_class("chr1") == "primary_reference"
    assert locus_class("KI270728.1") == "alternate_locus"
    assert locus_class("NW_012345.1") == "scaffold"
    assert locus_class("CHR_HSCHR6_MHC_MCF_PATCH") == "patch"


def test_stop_audit_does_not_normalize():
    assert star_class("MPEP*") == "terminal_stop"
    assert star_class("M*PEP") == "internal_stop"
    assert star_class("M*PEP*") == "multiple_stops"


def test_final_contract_exists_and_initial_outputs_are_preserved():
    root = Path(__file__).resolve().parents[1]
    refined = root / "data/processed/mapping_refined"
    initial = root / "data/processed/mapping"
    assert (refined / "gene_resolution.csv").exists()
    assert (refined / "esm_sequences.parquet").exists()
    assert (refined / "esm_sequence_relationships.parquet").exists()
    assert (initial / "gene_protein_mapping.parquet").exists()
    seq = pd.read_parquet(refined / "esm_sequences.parquet")
    assert seq["protein_sequence_hash"].is_unique
    assert not seq["sequence_has_stop"].any()
