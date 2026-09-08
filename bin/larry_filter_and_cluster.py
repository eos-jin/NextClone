#!/usr/bin/env python3
"""
LARRY barcode filtering and clustering for NextClone

Implements the complete LARRY filtering workflow from Weinreb et al. (Science 2020):
1. Extract (cell_bc, umi, larry_bc) tuples with read counts
2. Filter by minimum read count per tuple
3. Cluster LARRY barcodes by Hamming distance
4. Count UMIs per (cell, barcode) combination
5. Filter by minimum UMI count per (cell, barcode)
6. Output clone assignments per cell

Reference: Weinreb et al. (2020) Science
Lineage tracing on transcriptional landscapes links state to fate during differentiation
"""

import sys
import argparse
import gzip
from collections import defaultdict


def hamming_distance(bc1, bc2):
    """Calculate Hamming distance between two barcodes."""
    if len(bc1) != len(bc2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(bc1, bc2))


def cluster_barcodes(all_barcodes, max_hamming):
    """
    Cluster barcodes by Hamming distance.
    
    Sort barcodes by frequency (descending), then for each barcode:
    - If within max_hamming of a more abundant barcode, map it to that barcode
    - Otherwise, create a new cluster
    
    Returns:
        bc_map: dict mapping each barcode to its representative
        representatives: list of representative barcodes
    """
    # Sort by frequency (already sorted in input)
    good_bcs = []
    bc_map = {}
    
    for i, bc1 in enumerate(all_barcodes):
        if i > 0 and i % 500 == 0:
            print(f'Clustered {i} out of {len(all_barcodes)} barcodes', file=sys.stderr)
        
        mapped = False
        for bc2 in good_bcs:
            if hamming_distance(bc1, bc2) <= max_hamming:
                mapped = True
                bc_map[bc1] = bc2
                break
        
        if not mapped:
            good_bcs.append(bc1)
            bc_map[bc1] = bc1
    
    return bc_map, good_bcs


def larry_filter_and_cluster(input_file, output_file, sample_id,
                              min_reads=10, min_umis=3, max_hamming=3):
    """
    Complete LARRY filtering and clustering workflow.
    
    Args:
        input_file: FASTQ file with extracted LARRY barcodes
        output_file: Output CSV file for clone assignments
        sample_id: Sample identifier
        min_reads: Minimum reads per (cell, umi, barcode) tuple
        min_umis: Minimum UMIs per (cell, barcode) combination
        max_hamming: Maximum Hamming distance for clustering
    """
    print(f"=== LARRY Filtering and Clustering ===", file=sys.stderr)
    print(f"Input: {input_file}", file=sys.stderr)
    print(f"Parameters:", file=sys.stderr)
    print(f"  min_reads: {min_reads}", file=sys.stderr)
    print(f"  min_umis: {min_umis}", file=sys.stderr)
    print(f"  max_hamming: {max_hamming}", file=sys.stderr)
    
    # Step 1: Count (cell_bc, umi, larry_bc) tuples
    print(f"\nStep 1: Counting tuples...", file=sys.stderr)
    counts = defaultdict(int)
    
    opener = gzip.open if input_file.endswith('.gz') else open
    with opener(input_file, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # quality
            
            # Parse header: >@cell_bc,umi
            if not header.startswith('>@'):
                continue
            
            parts = header[2:].split(',')
            if len(parts) != 2:
                continue
            
            cell_bc = parts[0]
            umi = parts[1]
            larry_bc = seq
            
            counts[(cell_bc, umi, larry_bc)] += 1
    
    print(f"  Total unique tuples: {len(counts):,}", file=sys.stderr)
    
    # Step 2: Filter by minimum read count
    print(f"\nStep 2: Filtering by min_reads={min_reads}...", file=sys.stderr)
    counts_filtered = {k: v for k, v in counts.items() if v >= min_reads}
    print(f"  Retaining {len(counts_filtered):,} out of {len(counts):,} tuples", file=sys.stderr)
    
    # Step 3: Cluster LARRY barcodes by Hamming distance
    print(f"\nStep 3: Clustering barcodes (max_hamming={max_hamming})...", file=sys.stderr)
    all_larry_bcs = sorted(set(k[2] for k in counts_filtered.keys()))
    print(f"  Unique LARRY barcodes: {len(all_larry_bcs):,}", file=sys.stderr)
    
    bc_map, representatives = cluster_barcodes(all_larry_bcs, max_hamming)
    print(f"  Clustered to {len(representatives):,} representative barcodes", file=sys.stderr)
    
    # Step 4: Count UMIs per (cell, barcode) using clustered barcodes
    print(f"\nStep 4: Counting UMIs per cell...", file=sys.stderr)
    cell_data = defaultdict(lambda: defaultdict(int))
    
    for (cell_bc, umi, larry_bc), read_count in counts_filtered.items():
        # Map to representative barcode
        rep_bc = bc_map[larry_bc]
        # Count unique UMIs (each umi counts once, regardless of read count)
        cell_data[cell_bc][rep_bc] += 1
    
    print(f"  Cells with barcodes: {len(cell_data):,}", file=sys.stderr)
    
    # Step 5: Filter by minimum UMI count
    print(f"\nStep 5: Filtering by min_umis={min_umis}...", file=sys.stderr)
    final_bcs = {}
    for cell_bc, larry_bc_counts in cell_data.items():
        # Keep only barcodes with >= min_umis UMIs
        filtered = [bc for bc, umi_count in larry_bc_counts.items() if umi_count >= min_umis]
        final_bcs[cell_bc] = '-'.join(sorted(filtered))
    
    cells_with_barcodes = sum(1 for v in final_bcs.values() if v)
    unique_clones = len(set(v for v in final_bcs.values() if v))
    print(f"  Cells with final barcodes: {cells_with_barcodes:,}", file=sys.stderr)
    print(f"  Unique clones: {unique_clones:,}", file=sys.stderr)
    
    # Step 6: Write output
    print(f"\nStep 6: Writing output...", file=sys.stderr)
    with open(output_file, 'w') as f:
        f.write("sample,cell_bc,barcodes\n")
        for cell_bc, barcodes in final_bcs.items():
            f.write(f"{sample_id},{cell_bc},{barcodes}\n")
    
    print(f"Output written to: {output_file}", file=sys.stderr)
    print(f"\n=== Summary ===", file=sys.stderr)
    print(f"Total cells: {len(final_bcs):,}", file=sys.stderr)
    print(f"Cells with barcodes: {cells_with_barcodes:,}", file=sys.stderr)
    print(f"Unique clones: {unique_clones:,}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='LARRY barcode filtering and clustering for NextClone',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('input', help='Input FASTQ file with extracted LARRY barcodes')
    parser.add_argument('output', help='Output CSV file for clone assignments')
    parser.add_argument('--sample', default='sample1',
                       help='Sample identifier')
    parser.add_argument('--min-reads', type=int, default=10,
                       help='Minimum reads per (cell, umi, barcode) tuple')
    parser.add_argument('--min-umis', type=int, default=3,
                       help='Minimum UMIs per (cell, barcode) combination')
    parser.add_argument('--max-hamming', type=int, default=3,
                       help='Maximum Hamming distance for barcode clustering')
    
    args = parser.parse_args()
    
    larry_filter_and_cluster(
        args.input, args.output, args.sample,
        min_reads=args.min_reads,
        min_umis=args.min_umis,
        max_hamming=args.max_hamming
    )


if __name__ == '__main__':
    main()
