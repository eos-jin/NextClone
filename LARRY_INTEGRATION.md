# LARRY Integration for NextClone

This document describes the LARRY (Lineage And RNA Recovery) barcode extraction module integrated into NextClone.

## Overview

LARRY is a Cas9-based lineage tracing method from the Klein Lab (Weinreb et al., Science 2020). This module extracts LARRY barcodes from paired R1/R2 FASTQ files and processes them through NextClone's standard barcode analysis pipeline.

## LARRY Data Structure

**Input FASTQ files:**
- **R1**: Cell barcode (16bp) + UMI (8bp) from 10X scRNA-seq
- **R2**: LARRY barcode read containing:
  - Prefix: `GTTGCTAGGAGAGACCATATG`
  - Barcode: 40bp (with validation pattern at positions 4:6, 10:12, 16:18, 22:24, 28:30, 34:36)

## Usage

### Basic Command

```bash
nextflow run main.nf \
    --mode LARRY \
    --larry_r1_files "data/larry/*_R1_*.fastq.gz" \
    --larry_r2_files "data/larry/*_R2_*.fastq.gz" \
    --discovery_mode true \
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
| `--discovery_mode` | `false` | Use barcode discovery (true) or whitelist (false) |
| `--clone_barcodes_reference` | - | Required if discovery_mode=false |

### Workflow Steps

1. **larry_extract_barcodes**: Extract LARRY barcodes from R1/R2 pairs
   - Validates barcode structure (conserved positions)
   - Outputs FASTQ with cell_bc,UMI in header
   
2. **larry_count_barcodes**: Count unique barcodes
   - Outputs TSV: barcode\tcount
   
3. **larry_split_reads_to_chunks**: Split for parallel processing
   - Reuses DNAseq chunking logic
   
4. **Mapping** (discovery or whitelist mode):
   - Discovery: Two-pass barcode discovery + mapping
   - Whitelist: Direct mapping to known barcodes
   
5. **dnaseq_collapse_barcodes**: Final barcode collapsing and counting

## Output Files

```
results/larry/
├── *_larry_barcodes.fastq.gz    # Extracted barcodes
├── *_larry_stats.txt            # Extraction statistics
├── *_barcodes_counts.txt        # Barcode frequency counts
├── clone_barcode_counts.csv     # Final collapsed counts
└── nextclone_qc_report.html     # QC dashboard (if enabled)
```

## Validation

LARRY barcodes are validated using conserved positions:
- Position 4:6 = `TG`
- Position 10:12 = `CA`
- Position 16:18 = `AC`
- Position 22:24 = `GA`
- Position 28:30 = `GT`
- Position 34:36 = `AG`

This filters out sequencing errors and non-specific amplification.

## Testing

Test with synthetic data:

```bash
# Create test data
python3 tests/create_larry_test_data.py

# Run pipeline
nextflow run main.nf \
    --mode LARRY \
    --larry_r1_files "tests/larry_data/*_R1_*.fastq.gz" \
    --larry_r2_files "tests/larry_data/*_R2_*.fastq.gz" \
    --discovery_mode true
```

## Integration with NextClone

The LARRY module integrates seamlessly with NextClone's existing infrastructure:
- Reuses `dnaseq_split_reads.py` for chunking
- Reuses `dnaseq_map_barcodes` and `dnaseq_map_with_discovered_barcodes` for mapping
- Reuses `dnaseq_collapse_barcodes` for final counting
- Compatible with QC report generation

## Files Added

```
NextClone/
├── bin/
│   ├── larry_extract_barcodes.py    # Extraction script
│   └── larry_count_barcodes.py      # Counting script
├── modules/
│   └── extract_larry_barcodes.nf    # Nextflow processes
├── main.nf                          # Updated with LARRY workflow
└── nextflow.config                  # Updated with LARRY parameters
```

## Troubleshooting

### No barcodes extracted

Check:
1. R1/R2 files are not swapped
2. LARRY prefix matches your protocol variant
3. Data is actually from LARRY protocol
4. Try `--larry_prefix` with different sequence if using modified LARRY

### Low extraction rate

Possible causes:
- Poor sequencing quality (check R2 quality scores)
- Incorrect cell barcode length (adjust `--larry_cell_bc_len`)
- Prefix mutations (try `--no-validate` to skip structure check)

### Comparison with original LARRY pipeline

The original LARRY notebook (`LARRY_for_10X.ipynb`) performs:
1. Barcode extraction (similar to our `larry_extract_barcodes.py`)
2. Hamming distance clustering (NextClone uses Flexiplex edit distance)
3. UMI filtering (NextClone handles this in collapse step)

NextClone provides:
- Parallelized processing via Nextflow
- Integration with standard QC reporting
- Compatibility with other protocols (DNAseq, scRNAseq)
- Reproducible workflow management

## References

- Weinreb et al. (2020) Science: "Lineage tracing on a transcriptional landscape of hematopoietic development"
- LARRY GitHub: https://github.com/AllonKleinLab/LARRY
