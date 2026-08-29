[CmdletBinding()]
param(
    [string]$Root
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($Root)) { $Root = Join-Path $PSScriptRoot ".." }
$Root = [IO.Path]::GetFullPath($Root)
$Raw = Join-Path $Root "data/external/phase1"
$Model = Join-Path $Root "data/external/phase1/models"
$Arc = Join-Path $Root "data/external/phase1/arc"
New-Item -ItemType Directory -Force -Path $Raw,$Model,$Arc | Out-Null

function Get-Resource([string]$Url, [string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) {
        Write-Host "Downloading $Url"
        Invoke-WebRequest -Uri $Url -OutFile $Path
    } else {
        Write-Host "Already present: $Path"
    }
}

# Primary annotation and its matching protein translations.
Get-Resource "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz" (Join-Path $Raw "gencode.v50.annotation.gtf.gz")
Get-Resource "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_translations.fa.gz" (Join-Path $Raw "gencode.v50.pc_translations.fa.gz")
Get-Resource "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.metadata.RefSeq.gz" (Join-Path $Raw "gencode.v50.metadata.RefSeq.gz")

# MANE is a selection/cross-reference layer, not the primary universe.
Get-Resource "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz" (Join-Path $Raw "MANE.GRCh38.v1.5.summary.txt.gz")

# Official FAIR ESM-2 checkpoint.
Get-Resource "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt" (Join-Path $Model "esm2_t33_650M_UR50D.pt")

# Arc/STATE baseline: only the published human protein-embedding artifact,
# pinned to an immutable Hugging Face revision. Do not load or transform it.
$ArcRevision = "9853106c9d6d5a60c0eb106f836ceab29f0c99"
Get-Resource "https://huggingface.co/arcinstitute/SE-600M/resolve/$ArcRevision/protein_embeddings.pt?download=true" (Join-Path $Arc "protein_embeddings.pt")
Get-Resource "https://huggingface.co/arcinstitute/SE-600M/resolve/$ArcRevision/README.md?download=true" (Join-Path $Arc "README.md")
Get-Resource "https://huggingface.co/arcinstitute/SE-600M/resolve/$ArcRevision/LICENSE.md?download=true" (Join-Path $Arc "LICENSE.md")
Get-Resource "https://huggingface.co/arcinstitute/SE-600M/resolve/$ArcRevision/MODEL_LICENSE.md?download=true" (Join-Path $Arc "MODEL_LICENSE.md")
Get-Resource "https://huggingface.co/arcinstitute/SE-600M/resolve/$ArcRevision/MODEL_ACCEPTABLE_USE_POLICY.md?download=true" (Join-Path $Arc "MODEL_ACCEPTABLE_USE_POLICY.md")
Get-Resource "https://huggingface.co/arcinstitute/SE-600M/resolve/$ArcRevision/config.yaml?download=true" (Join-Path $Arc "config.yaml")

$Files = Get-ChildItem -LiteralPath (Join-Path $Root "data/external/phase1") -File -Recurse | Sort-Object FullName
$Lines = foreach ($File in $Files) {
    $Hash = (Get-FileHash -LiteralPath $File.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    $Relative = [IO.Path]::GetRelativePath($Root, $File.FullName).Replace("\\", "/")
    "$Hash  $Relative"
}
Set-Content -LiteralPath (Join-Path $Root "data/external/phase1/SHA256SUMS") -Value $Lines -Encoding ascii
Write-Host "Wrote SHA256SUMS for $($Files.Count) files."
