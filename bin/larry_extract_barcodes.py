#!/usr/bin/env python3
"""
LARRY barcode extraction for NextClone

Extracts LARRY (Lineage And RNA Recovery) barcodes from paired R1/R2 FASTQ files.
LARRY is a Cas9-based lineage tracing method (Weinreb et al., Science 2020).

Input:
  R1 FASTQ: Cell barcode (16bp) + UMI (8bp) from 10X scRNA-seq
  R2 FASTQ: LARRY barcode read containing prefix GTTGCTAGGAGAGACCATATG + 40bp barcode

Output:
  FASTQ file with extracted LARRY barcodes.
  Header format: >@cell_bc,UMI
  Sequence: 40bp LARRY barcode (after prefix)

This output feeds directly into NextClone's DNAseq workflow for barcode
counting, collapsing, and clone identification.

For cell-level clone assignment, the cell_bc and UMI in the header are
preserved through the pipeline for downstream UMI deduplication.
"""

import gzip
import sys
import argparse
from pathlib import Path


# LARRY barcode validation pattern
# The LARRY barcode has conserved bases at specific positions that serve as
# quality checkpoints. This filters out non-specific amplification and sequencing errors.
LARRY_CHECK_POSITIONS = {
    (4, 6): 'TG',
    (10, 12): 'CA',
    (16, 18): 'AC',
    (22, 24): 'GA',
    (28, 30): 'GT',
    (34, 36): 'AG',
}


def is_valid_larry_barcode(bc):
    """Validate LARRY barcode structure using conserved position checks."""
    if len(bc) < 40:
        return False
    return all(bc[start:end] == expected for (start, end), expected in LARRY_CHECK_POSITIONS.items())


def open_fastq(path, mode='rt'):
    """Open a FASTQ file, handling both gzipped and plain text."""
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    else:
        return open(path, mode.replace('t', '').replace('b', '') + 't' if 'b' not in mode else mode)


def extract_larry_barcodes(r1_path, r2_path, output_path, larry_prefix='GTTGCTAGGAGAGACCATATG',
                           cell_bc_len=16, umi_len=8, larry_bc_len=40, validate=True):
    """
    Extract LARRY barcodes from paired R1/R2 FASTQ files.
    
    For each read pair:
    1. Extract cell_bc (first 16bp of R1) and UMI (next 8bp of R1)
    2. Search R2 for the LARRY prefix
    3. Extract 40bp barcode after prefix
    4. Validate barcode structure (optional)
    5. Write output FASTQ: >@cell_bc,UMI / barcode_sequence
    
    Returns dict with extraction statistics.
    """
    stats = {
        'total_reads': 0,
        'reads_with_prefix': 0,
        'valid_barcodes': 0,
        'written_reads': 0
    }
    
    r1_handle = gzip.open(r1_path, 'rt') if str(r1_path).endswith('.gz') else open(r1_path, 'r')
    r2_handle = gzip.open(r2_path, 'rt') if str(r2_path).endswith('.gz') else open(r2_path, 'r')
    out_handle = gzip.open(output_path, 'wt') if str(output_path).endswith('.gz') else open(output_path, 'w')
    
    try:
        while True:
            # Read R1 (4 lines per FASTQ record)
            r1_header = r1_handle.readline().strip()
            r1_seq = r1_handle.readline().strip()
            r1_plus = r1_handle.readline().strip()
            r1_qual = r1_handle.readline().strip()
            
            # Read R2
            r2_header = r2_handle.readline().strip()
            r2_seq = r2_handle.readline().strip()
            r2_plus = r2_handle.readline().strip()
            r2_qual = r2_handle.readline().strip()
            
            # End of file
            if not r1_header or not r2_header:
                break
            
            stats['total_reads'] += 1
            
            # Skip malformed records
            if not r1_header.startswith('@') or not r2_header.startswith('@'):
                continue
            
            # Extract cell barcode and UMI from R1
            if len(r1_seq) < cell_bc_len + umi_len:
                continue
            
            cell_bc = r1_seq[:cell_bc_len]
            umi = r1_seq[cell_bc_len:cell_bc_len + umi_len]
            
            # Search for LARRY prefix in R2
            prefix_pos = r2_seq.find(larry_prefix)
            if prefix_pos == -1:
                continue
            
            stats['reads_with_prefix'] += 1
            
            # Extract LARRY barcode (40bp after prefix)
            larry_start = prefix_pos + len(larry_prefix)
            larry_bc = r2_seq[larry_start:larry_start + larry_bc_len]
            
            if len(larry_bc) < larry_bc_len:
                continue
            
            # Validate LARRY barcode structure
            if validate and not is_valid_larry_barcode(larry_bc):
                continue
            
            stats['valid_barcodes'] += 1
            
            # Write output FASTQ
            # Format: >@cell_bc,UMI (compatible with NextClone parsing)
            out_handle.write(f">@{cell_bc},{umi}\n")
            out_handle.write(f"{larry_bc}\n")
            out_handle.write("+\n")
            # Use R2 quality scores for the barcode region
            larry_qual = r2_qual[larry_start:larry_start + larry_bc_len]
            if len(larry_qual) < larry_bc_len:
                larry_qual = 'I' * larry_bc_len
            out_handle.write(f"{larry_qual}\n")
            
            stats['written_reads'] += 1
    
    finally:
        r1_handle.close()
        r2_handle.close()
        out_handle.close()
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Extract LARRY barcodes from paired R1/R2 FASTQ files for NextClone',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--r1', required=True, help='R1 FASTQ file (cell_bc + UMI)')
    parser.add_argument('--r2', required=True, help='R2 FASTQ file (LARRY barcode read)')
    parser.add_argument('--output', required=True, help='Output FASTQ file')
    parser.add_argument('--prefix', default='GTTGCTAGGAGAGACCATATG',
                       help='LARRY barcode prefix sequence')
    parser.add_argument('--cell-bc-len', type=int, default=16,
                       help='Cell barcode length in R1')
    parser.add_argument('--umi-len', type=int, default=8,
                       help='UMI length in R1')
    parser.add_argument('--larry-bc-len', type=int, default=40,
                       help='LARRY barcode length to extract')
    parser.add_argument('--no-validate', action='store_true',
                       help='Skip LARRY barcode structure validation')
    
    args = parser.parse_args()
    
    # Validate input files
    for path, name in [(args.r1, 'R1'), (args.r2, 'R2')]:
        if not Path(path).exists():
            print(f"ERROR: {name} file not found: {path}", file=sys.stderr)
            sys.exit(1)
    
    # Extract barcodes
    print(f"Extracting LARRY barcodes from:", file=sys.stderr)
    print(f"  R1: {args.r1}", file=sys.stderr)
    print(f"  R2: {args.r2}", file=sys.stderr)
    print(f"  Prefix: {args.prefix}", file=sys.stderr)
    print(f"  Barcode length: {args.larry_bc_len}bp", file=sys.stderr)
    print(f"  Validation: {'disabled' if args.no_validate else 'enabled'}", file=sys.stderr)
    
    stats = extract_larry_barcodes(
        args.r1, args.r2, args.output,
        larry_prefix=args.prefix,
        cell_bc_len=args.cell_bc_len,
        umi_len=args.umi_len,
        larry_bc_len=args.larry_bc_len,
        validate=not args.no_validate
    )
    
    # Print statistics
    print(f"\n=== LARRY Extraction Statistics ===", file=sys.stderr)
    print(f"Total read pairs processed: {stats['total_reads']:,}", file=sys.stderr)
    print(f"Reads with LARRY prefix:    {stats['reads_with_prefix']:,} "
          f"({100*stats['reads_with_prefix']/max(1,stats['total_reads']):.1f}%)", file=sys.stderr)
    print(f"Valid barcodes extracted:   {stats['valid_barcodes']:,} "
          f"({100*stats['valid_barcodes']/max(1,stats['reads_with_prefix']):.1f}% of prefix matches)", file=sys.stderr)
    print(f"Written to output:          {stats['written_reads']:,}", file=sys.stderr)
    
    if stats['written_reads'] == 0:
        print("\nERROR: No valid LARRY barcodes found!", file=sys.stderr)
        print("Troubleshooting:", file=sys.stderr)
        print("  1. Verify R1/R2 files are not swapped", file=sys.stderr)
        print("  2. Check that LARRY prefix matches your protocol variant", file=sys.stderr)
        print("  3. Try --no-validate to skip structure validation", file=sys.stderr)
        print("  4. Confirm data is from LARRY protocol (not indrops)", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nOutput written to: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
