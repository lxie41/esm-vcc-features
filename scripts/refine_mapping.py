from pathlib import Path
from protein_embeddings.refinement import refine

root = Path(__file__).resolve().parents[1]
refine(root / "data/processed/mapping", root / "data/processed/mapping_refined", root / "data/gene_names.csv", root / "data/raw/hgnc/current/hgnc_complete_set.txt", root / "data/raw/hgnc/current/withdrawn.txt", root / "data/raw/gencode/v50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz", root / "data/raw/gencode/v50/gencode.v50.pc_translations.fa.gz", root / "data/raw/mane/v1.5/MANE.GRCh38.v1.5.summary.txt.gz")
