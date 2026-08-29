[CmdletBinding()]
param([string]$Root)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($Root)) { $Root = Join-Path $PSScriptRoot ".." }
$Root = [IO.Path]::GetFullPath($Root)

$Raw = Join-Path $Root "data/raw"
$Gencode = Join-Path $Raw "gencode/v50"
$Mane = Join-Path $Raw "mane/v1.5"
New-Item -ItemType Directory -Force -Path $Gencode,$Mane | Out-Null

function Download-Verified([string]$Url, [string]$Destination) {
    $part = "$Destination.part"
    if (-not (Test-Path -LiteralPath $Destination)) {
        if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
        Write-Host "Downloading $Url"
        Invoke-WebRequest -Uri $Url -OutFile $part
        if (-not (Test-Path -LiteralPath $part)) { throw "Download did not create $part" }
        Move-Item -LiteralPath $part -Destination $Destination
    } else {
        Write-Host "Already present: $Destination"
    }
}

function Validate-Gzip([string]$Path, [string]$Kind) {
    $inputStream = [IO.File]::OpenRead($Path)
    try {
        $gzip = [IO.Compression.GzipStream]::new($inputStream, [IO.Compression.CompressionMode]::Decompress)
        try {
            $reader = [IO.StreamReader]::new($gzip)
            try {
                $firstHeaders = [Collections.Generic.List[string]]::new()
                $firstRecords = [Collections.Generic.List[string]]::new()
                $header = $null
                $hasSequence = $false
                $lineCount = 0
                while ($null -ne ($line = $reader.ReadLine())) {
                    $lineCount++
                    if ($line.Trim()) {
                        if ($line.StartsWith('>')) {
                            if ($firstHeaders.Count -lt 5) { $firstHeaders.Add($line) }
                        } elseif ($line -notmatch '^#' -and $firstRecords.Count -lt 100) {
                            $firstRecords.Add($line)
                        }
                        if (-not $header -and $Kind -eq 'MANE_SUMMARY') { $header = $line }
                        if ($Kind -eq 'GENCODE_FASTA' -and $line -notmatch '^>') { $hasSequence = $true }
                    }
                }
            } finally { $reader.Dispose() }
        } finally { $gzip.Dispose() }
    } finally { $inputStream.Dispose() }
    if ($lineCount -eq 0) { throw "$Kind decompressed to empty text: $Path" }
    switch ($Kind) {
        "GENCODE_GTF" {
            $hasGene = $false
            foreach ($record in $firstRecords) { if (($record -split "`t").Count -ge 3 -and ($record -split "`t")[2] -eq 'gene') { $hasGene = $true; break } }
            if (-not $hasGene) { throw "No GTF gene records found in first records" }
            $allRecords = $firstRecords -join "`n"
            foreach ($field in @('gene_id','transcript_id','gene_name','gene_type','transcript_type')) {
                if ($allRecords.IndexOf($field, [StringComparison]::Ordinal) -lt 0) { throw "Missing expected GTF attribute $field" }
            }
            return @{ format = "GTF"; decompressed_line_count = $lineCount; structural_check = "gzip_eof_reached_gene_records_and_expected_attributes_present" }
        }
        "GENCODE_FASTA" {
            if ($firstHeaders.Count -eq 0 -or -not $hasSequence) { throw "No FASTA headers and sequence records found" }
            return @{ format = "FASTA"; decompressed_line_count = $lineCount; structural_check = "gzip_eof_reached_protein_fasta_headers_and_sequence_present"; first_headers = @($firstHeaders) }
        }
        "MANE_SUMMARY" {
            if ($header -notmatch 'MANE|Ensembl|RefSeq|Gene') { throw "MANE header does not contain expected identifiers: $header" }
            return @{ format = "TABULAR"; decompressed_line_count = $lineCount; structural_check = "gzip_eof_reached_header_present"; header = $header }
        }
    }
}

$resources = @(
    [ordered]@{ resource_name = "GENCODE human comprehensive ALL annotation"; release = "GENCODE 50 / Ensembl 116"; assembly = "GRCh38.p14"; filename = "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"; source_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"; path = (Join-Path $Gencode "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"); purpose = "gene, transcript, metadata, biotypes, and stable relationships"; kind = "GENCODE_GTF" },
    [ordered]@{ resource_name = "GENCODE protein-coding translations"; release = "GENCODE 50 / Ensembl 116"; assembly = "GRCh38.p14"; filename = "gencode.v50.pc_translations.fa.gz"; source_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_translations.fa.gz"; path = (Join-Path $Gencode "gencode.v50.pc_translations.fa.gz"); purpose = "authoritative future amino-acid sequence input"; kind = "GENCODE_FASTA" },
    [ordered]@{ resource_name = "MANE human summary"; release = "MANE v1.5"; assembly = "GRCh38"; filename = "MANE.GRCh38.v1.5.summary.txt.gz"; source_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz"; path = (Join-Path $Mane "MANE.GRCh38.v1.5.summary.txt.gz"); purpose = "MANE Select/Plus Clinical and Ensembl-RefSeq cross-reference metadata"; kind = "MANE_SUMMARY" }
)

$manifest = [ordered]@{
    snapshot_name = "VCC mapping reference resources"
    snapshot_date = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
    vcc_inputs = @(
        [ordered]@{ resource_name = "VCC gene vocabulary"; path = "data/gene_names.csv"; sha256 = (Get-FileHash (Join-Path $Root "data/gene_names.csv") -Algorithm SHA256).Hash.ToLowerInvariant(); file_size_bytes = (Get-Item (Join-Path $Root "data/gene_names.csv")).Length; purpose = "official 18,533-gene VCC universe; preserved in place" },
        [ordered]@{ resource_name = "VCC perturbation targets"; path = "data/pert_counts.csv"; sha256 = (Get-FileHash (Join-Path $Root "data/pert_counts.csv") -Algorithm SHA256).Hash.ToLowerInvariant(); file_size_bytes = (Get-Item (Join-Path $Root "data/pert_counts.csv")).Length; purpose = "official 300 perturbation targets; not the full universe" }
    )
    resources = @()
}

foreach ($resource in $resources) {
    Download-Verified $resource.source_url $resource.path
    $check = Validate-Gzip $resource.path $resource.kind
    $item = [ordered]@{}
    foreach ($key in @('resource_name','release','assembly','filename','source_url','purpose')) { $item[$key] = $resource[$key] }
    $relativePath = $resource.path.Substring($Root.Length)
    $item.path = $relativePath.TrimStart([char[]]'\\/').Replace('\', '/')
    $item.sha256 = (Get-FileHash -LiteralPath $resource.path -Algorithm SHA256).Hash.ToLowerInvariant()
    $item.file_size_bytes = (Get-Item -LiteralPath $resource.path).Length
    $item.download_date = (Get-Item -LiteralPath $resource.path).LastWriteTimeUtc.ToString("yyyy-MM-ddTHH:mm:ssZ")
    $item.validation = $check
    $manifest.resources += $item
}

$manifestPath = Join-Path $Raw "resource_manifest.json"
$manifest | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $manifestPath -Encoding utf8
Write-Host "Wrote $manifestPath"
