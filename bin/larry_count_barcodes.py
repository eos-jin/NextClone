#!/usr/bin/env python3
"""
Count unique LARRY barcodes and their frequencies.

This script reads the extracted LARRY barcodes from a FASTQ file
and outputs a tab-separated file with barcode sequences and their counts.

Output format:
    barcode\tcount

This format is compatible with NextClone's DNAseq workflow for downstream
barcode collapsing, filtering, and clone assignment.
"""

import gzip
import sys
import argparse
from pathlib import Path
from collections import Counter


def count_barcodes(input_fastq, output_file):
    """
    Count unique barcodes from a FASTQ file.
    
    Args:
        input_fastq: Path to input FASTQ file (gzipped or plain)
        output_file: Path to output TSV file
    
    Returns:
        dict with counting statistics
    """
    stats = {
        'total_reads': 0,
        'unique_barcodes': 0,
        'max_count': 0,
        'min_count': float('inf')
    }
    
    # Count barcodes
    barcode_counts = Counter()
    
    # Open input file
    if str(input_fastq).endswith('.gz'):
        handle = gzip.open(input_fastq, 'rt')
    else:
        handle = open(input_fastq, 'r')
    
    try:
        while True:
            # Read FASTQ record (4 lines)
            header = handle.readline().strip()
            if not header:
                break
            
            seq = handle.readline().strip()
            plus = handle.readline().strip()
            qual = handle.readline().strip()
            
            if not seq:
                break
            
            stats['total_reads'] += 1
            barcode_counts[seq] += 1
    
    finally:
        handle.close()
    
    # Calculate statistics
    stats['unique_barcodes'] = len(barcode_counts)
    if barcode_counts:
        stats['max_count'] = max(barcode_counts.values())
        stats['min_count'] = min(barcode_counts.values())
    else:
        stats['min_count'] = 0
    
    # Write output (sorted by count, descending)
    with open(output_file, 'w') as out:
        for barcode, count in sorted(barcode_counts.items(), key=lambda x: -x[1]):
            out.write(f'{barcode}\t{count}\n')
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Count unique LARRY barcodes from extracted FASTQ file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--input', required=True, help='Input FASTQ file (extracted LARRY barcodes)')
    parser.add_argument('--output', required=True, help='Output TSV file (barcode\\tcount)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input).exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Count barcodes
    print(f"Counting barcodes from: {args.input}", file=sys.stderr)
    stats = count_barcodes(args.input, args.output)
    
    # Print statistics
    print(f"\n=== Barcode Counting Statistics ===", file=sys.stderr)
    print(f"Total reads processed:  {stats['total_reads']:,}", file=sys.stderr)
    print(f"Unique barcodes found:  {stats['unique_barcodes']:,}", file=sys.stderr)
    print(f"Max count per barcode:  {stats['max_count']:,}", file=sys.stderr)
    print(f"Min count per barcode:  {stats['min_count']:,}", file=sys.stderr)
    
    if stats['unique_barcodes'] == 0:
        print("\nWARNING: No barcodes found!", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nOutput written to: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
