"""Thin entry point for the frozen VCC mapping stage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from protein_embeddings.mapping import build_mapping


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    root = parser.parse_args().root.resolve()
    qc = build_mapping(
        root / "data/gene_names.csv", root / "data/pert_counts.csv",
        root / "data/raw/gencode/v50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz",
        root / "data/raw/gencode/v50/gencode.v50.pc_translations.fa.gz",
        root / "data/raw/mane/v1.5/MANE.GRCh38.v1.5.summary.txt.gz",
        root / "data/processed/mapping",
    )
    print(json.dumps(qc, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
