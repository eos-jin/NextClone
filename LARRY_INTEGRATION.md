# LARRY Integration for NextClone

This document describes the LARRY (Lineage And RNA Recovery) barcode extraction module integrated into NextClone.

## Overview

LARRY is a Cas9-based lineage tracing method from the Klein Lab (Weinreb et al., Science 2020). This module extracts LARRY barcodes from paired R1/R2 FASTQ files and assigns cells to clones.

## How LARRY Differs from Other Modes

**LARRY does NOT use discovery mode or whitelist mode.**

Unlike DNAseq/scRNAseq where Flexiplex discovers unknown barcodes from the data, LARRY barcodes are always extracted by finding reads with a fixed prefix sequence (`GTTGCTAGGAGAGACCATATG`). There is no barcode "discovery" step — we simply extract everything that has the LARRY prefix, then filter and cluster the extracted barcodes.

The LARRY workflow is:
1. **Extract**: Find all reads with the LARRY prefix in R2, extract cell_bc + UMI from R1
2. **Filter reads**: Remove low-confidence (cell, umi, barcode) tuples
3. **Cluster barcodes**: Group similar LARRY barcodes by Hamming distance
4. **Count UMIs**: Count unique UMIs per (cell, barcode) combination
5. **Filter UMIs**: Remove barcodes with too few UMIs per cell
6. **Output**: Clone assignments per cell

## LARRY Data Structure

**Input FASTQ files:**
- **R1**: Cell barcode (16bp) + UMI (8bp) from 10X scRNA-seq
- **R2**: LARRY barcode read containing:
  - Prefix: `GTTGCTAGGAGAGACCATATG` (21bp)
  - Barcode: 40bp (with validation pattern at positions 4:6, 10:12, 16:18, 22:24, 28:30, 34:36)

## Usage

### Basic Command

```bash
nextflow run main.nf \
    --mode LARRY \
    --larry_r1_files "data/larry/*_R1_*.fastq.gz" \
    --larry_r2_files "data/larry/*_R2_*.fastq.gz" \
    --publish_dir results/larry
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--larry_r1_files` | `data/larry/*_R1_*.fastq.gz` | Glob pattern for R1 FASTQ files |
| `--larry_r2_files` | `data/larry/*_R2_*.fastq.gz` | Glob pattern for R2 FASTQ files |
| `--larry_prefix` | `GTTGCTAGGAGAGACCATATG` | LARRY barcode prefix sequence |
| `--larry_cell_bc_len` | `16` | Cell barcode length in R1 |
| `--larry_umi_len` | `8` | UMI length in R1 |
| `--larry_bc_len` | `40` | LARRY barcode length |
| `--larry_min_reads` | `10` | Minimum reads per (cell, umi, barcode) tuple |
| `--larry_min_umis` | `3` | Minimum UMIs per (cell, barcode) combination |
| `--larry_max_hamming` | `3` | Maximum Hamming distance for barcode clustering |

### Filtering Parameters Explained

These parameters come from the original LARRY paper (Weinreb et al., Science 2020):

- **`larry_min_reads` (default: 10)**: Each unique (cell, umi, barcode) combination must have at least this many reads. This removes low-confidence extractions and sequencing noise.

- **`larry_min_umis` (default: 3)**: Each (cell, barcode) combination must be supported by at least this many unique UMIs. This removes PCR duplicates and ensures the barcode is genuinely present in the cell.

- **`larry_max_hamming` (default: 3)**: LARRY barcodes within this Hamming distance are clustered together (mapped to the more abundant barcode). This corrects for sequencing errors in the 40bp barcode region.

## Workflow Steps

### Step 1: Extract (`larry_extract_barcodes.py`)
- Reads paired R1/R2 FASTQ files
- Extracts cell_bc (16bp) and UMI (8bp) from R1
- Searches R2 for the LARRY prefix
- Extracts 40bp barcode after prefix
- Validates barcode structure (conserved positions)
- Outputs FASTQ with `>@cell_bc,UMI` header format

### Step 2: Filter & Cluster (`larry_filter_and_cluster.py`)
- Counts (cell_bc, umi, larry_bc) tuples
- Filters by minimum read count per tuple
- Clusters LARRY barcodes by Hamming distance
- Counts unique UMIs per (cell, barcode)
- Filters by minimum UMI count per (cell, barcode)
- Outputs clone assignments: `sample,cell_bc,barcodes`

## Output Files

```
results/larry/
├── *_larry_barcodes.fastq.gz       # Extracted barcodes (intermediate)
├── *_larry_stats.txt               # Extraction statistics
├── sample1_larry_clones.csv        # Final clone assignments
└── sample1_larry_filter_stats.txt  # Filtering and clustering statistics
```

The clones CSV format:
```
sample,cell_bc,barcodes
sample1,ACGTACGTACGTACGT-1,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
sample1,CGTACGTACGTACGTA-1,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

Cells with multiple barcodes have them joined with `-`.

## Barcode Validation

LARRY barcodes are validated using conserved positions (from the original LARRY notebook):
- Position 4:6 = `TG`
- Position 10:12 = `CA`
- Position 16:18 = `AC`
- Position 22:24 = `GA`
- Position 28:30 = `GT`
- Position 34:36 = `AG`

This filters out non-specific amplification and sequencing errors. Disable with `--no-validate` in the extraction step if needed.

## Comparison with Original LARRY Pipeline

The original LARRY notebook (`LARRY_for_10X.ipynb`) performs the same steps:
1. Barcode extraction with prefix search and validation
2. Read count filtering (N_READS = 10)
3. Hamming distance clustering (N_HAMMING = 3)
4. UMI counting per cell per barcode
5. UMI filtering (N_UMIS = 3)
6. Output clone assignments

NextClone provides:
- Parallelized processing via Nextflow
- Configurable parameters (not hardcoded)
- Integration with the broader NextClone ecosystem
- Reproducible workflow management

## Troubleshooting

### No barcodes extracted

Check:
1. R1/R2 files are not swapped
2. LARRY prefix matches your protocol variant
3. Data is actually from LARRY protocol
4. Try `--no-validate` to skip structure check

### Low extraction rate

Possible causes:
- Poor sequencing quality (check R2 quality scores)
- Incorrect cell barcode length (adjust `--larry_cell_bc_len`)
- Prefix mutations (try `--no-validate`)

### Too few cells with barcodes

Try relaxing the filtering thresholds:
- `--larry_min_reads 5` (instead of 10)
- `--larry_min_umis 2` (instead of 3)

### Too many unique clones (possible over-clustering)

Try increasing the Hamming distance:
- `--larry_max_hamming 4` or `5`

## Files Added

```
NextClone/
├── bin/
│   ├── larry_extract_barcodes.py       # Extraction script
│   └── larry_filter_and_cluster.py     # Filtering and clustering script
├── modules/
│   ├── extract_larry_barcodes.nf       # Nextflow process: extraction
│   └── larry_filter_and_cluster.nf     # Nextflow process: filtering
├── main.nf                             # Updated with LARRY workflow
└── nextflow.config                     # Updated with LARRY parameters
```

## References

- Weinreb et al. (2020) Science: "Lineage tracing on transcriptional landscapes links state to fate during differentiation"
- LARRY GitHub: https://github.com/AllonKleinLab/LARRY
