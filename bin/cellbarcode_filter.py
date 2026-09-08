#!/usr/bin/env python3
"""
CellBarcode-style filtering for NextClone discovery mode

Implements filtering strategies from the CellBarcode paper:
1. Automatic threshold filtering (k-means clustering on log-transformed counts)
2. Cluster filtering (remove barcodes similar to more abundant ones)
3. UMI filtering (if applicable)

Reference: Sun et al. (2024) Nature Computational Science
"""

import sys
import argparse
import numpy as np
from collections import Counter
from typing import List, Tuple, Dict


def auto_threshold_filter(barcode_counts: Dict[str, int], min_count: int = 1) -> Tuple[List[str], int]:
    """
    Automatic threshold filtering using 1D k-means clustering.
    
    Algorithm from CellBarcode:
    1. Remove barcodes with count below median
    2. Transform counts by log2(x+1)
    3. Apply 1D k-means clustering with k=2
    4. Use boundary between clusters as threshold
    
    Args:
        barcode_counts: Dict mapping barcode -> count
        min_count: Minimum count to consider
    
    Returns:
        Tuple of (filtered_barcodes, threshold_used)
    """
    # Filter out very low count barcodes
    filtered = {bc: count for bc, count in barcode_counts.items() if count >= min_count}
    
    if len(filtered) < 2:
        return list(filtered.keys()), min_count
    
    # Get counts and compute median
    counts = np.array(list(filtered.values()))
    median_count = np.median(counts)
    
    # Remove barcodes below median
    above_median = {bc: count for bc, count in filtered.items() if count >= median_count}
    
    if len(above_median) < 2:
        return list(filtered.keys()), min_count
    
    # Log transform
    log_counts = np.log2(np.array(list(above_median.values())) + 1)
    
    # 1D k-means with k=2
    # Initialize with min and max
    centers = np.array([log_counts.min(), log_counts.max()])
    
    # Iterate until convergence
    for _ in range(100):
        # Assign to nearest center
        assignments = np.argmin(np.abs(log_counts[:, None] - centers[None, :]), axis=1)
        
        # Update centers
        new_centers = np.array([
            log_counts[assignments == 0].mean() if np.any(assignments == 0) else centers[0],
            log_counts[assignments == 1].mean() if np.any(assignments == 1) else centers[1]
        ])
        
        # Check convergence
        if np.allclose(centers, new_centers):
            break
        centers = new_centers
    
    # Threshold is the midpoint between centers
    threshold_log = (centers[0] + centers[1]) / 2
    threshold = int(2 ** threshold_log - 1)
    
    # Filter barcodes
    filtered_barcodes = [bc for bc, count in filtered.items() if count >= threshold]
    
    return filtered_barcodes, threshold


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError("Strings must have equal length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def cluster_filter(barcode_counts: Dict[str, int], max_distance: int = 1, 
                   distance_metric: str = 'hamming') -> List[str]:
    """
    Cluster filtering: remove barcodes with small edit distance to more abundant ones.
    
    Algorithm:
    1. Sort barcodes by abundance (descending)
    2. For each barcode, compare to all more abundant barcodes
    3. If distance < threshold to any more abundant barcode, remove it
    
    Args:
        barcode_counts: Dict mapping barcode -> count
        max_distance: Maximum edit distance to consider as cluster member
        distance_metric: 'hamming' or 'levenshtein'
    
    Returns:
        List of filtered barcodes (true barcodes)
    """
    # Sort by count descending
    sorted_barcodes = sorted(barcode_counts.items(), key=lambda x: -x[1])
    
    kept_barcodes = []
    
    for barcode, count in sorted_barcodes:
        # Check if this barcode is similar to any already-kept barcode
        is_cluster_member = False
        
        for kept_bc in kept_barcodes:
            if distance_metric == 'hamming':
                if len(barcode) == len(kept_bc):
                    dist = hamming_distance(barcode, kept_bc)
                else:
                    # Different lengths, skip hamming check
                    continue
            else:
                # For now, use hamming only
                # Could implement levenshtein if needed
                continue
            
            if dist <= max_distance:
                is_cluster_member = True
                break
        
        if not is_cluster_member:
            kept_barcodes.append(barcode)
    
    return kept_barcodes


def manual_threshold_filter(barcode_counts: Dict[str, int], threshold: int) -> List[str]:
    """
    Manual threshold filtering: keep barcodes with count >= threshold.
    
    Args:
        barcode_counts: Dict mapping barcode -> count
        threshold: Minimum count to keep
    
    Returns:
        List of filtered barcodes
    """
    return [bc for bc, count in barcode_counts.items() if count >= threshold]


def apply_filtering(barcode_counts_file: str, output_file: str, 
                    filter_method: str = 'auto', threshold: int | None = None,
                    cluster_distance: int = 1, min_count: int = 1) -> None:
    """
    Apply CellBarcode-style filtering to discovered barcodes.
    
    Args:
        barcode_counts_file: Input file with barcode\\tcount format
        output_file: Output file for filtered barcodes
        filter_method: 'auto', 'manual', 'cluster', or 'combined'
        threshold: Manual threshold (required if filter_method='manual')
        cluster_distance: Max edit distance for cluster filtering
        min_count: Minimum count to consider
    """
    # Read barcode counts
    barcode_counts = {}
    with open(barcode_counts_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                barcode, count = parts
                barcode_counts[barcode] = int(count)
    
    print(f"Read {len(barcode_counts)} barcodes from {barcode_counts_file}", file=sys.stderr)
    
    # Apply filtering
    if filter_method == 'auto':
        filtered_barcodes, threshold_used = auto_threshold_filter(barcode_counts, min_count)
        print(f"Auto threshold: {threshold_used}", file=sys.stderr)
    
    elif filter_method == 'manual':
        if threshold is None:
            raise ValueError("Manual threshold must be specified")
        filtered_barcodes = manual_threshold_filter(barcode_counts, threshold)
        print(f"Manual threshold: {threshold}", file=sys.stderr)
    
    elif filter_method == 'cluster':
        filtered_barcodes = cluster_filter(barcode_counts, cluster_distance)
        print(f"Cluster filtering with max distance: {cluster_distance}", file=sys.stderr)
    
    elif filter_method == 'combined':
        # First apply auto threshold, then cluster filtering
        threshold_barcodes, threshold_used = auto_threshold_filter(barcode_counts, min_count)
        threshold_counts = {bc: barcode_counts[bc] for bc in threshold_barcodes}
        filtered_barcodes = cluster_filter(threshold_counts, cluster_distance)
        print(f"Combined: auto threshold {threshold_used}, then cluster distance {cluster_distance}", 
              file=sys.stderr)
    
    else:
        raise ValueError(f"Unknown filter method: {filter_method}")
    
    # Write output
    with open(output_file, 'w') as f:
        for barcode in filtered_barcodes:
            f.write(f"{barcode}\n")
    
    print(f"Filtered to {len(filtered_barcodes)} barcodes", file=sys.stderr)
    print(f"Output written to {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='CellBarcode-style filtering for NextClone discovery mode',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('input', help='Input barcode counts file (barcode\\tcount)')
    parser.add_argument('output', help='Output filtered barcodes file')
    parser.add_argument('--method', choices=['auto', 'manual', 'cluster', 'combined'],
                       default='auto', help='Filtering method')
    parser.add_argument('--threshold', type=int, help='Manual threshold (for method=manual)')
    parser.add_argument('--cluster-distance', type=int, default=1,
                       help='Max edit distance for cluster filtering')
    parser.add_argument('--min-count', type=int, default=1,
                       help='Minimum count to consider')
    
    args = parser.parse_args()
    
    apply_filtering(
        args.input, args.output,
        filter_method=args.method,
        threshold=args.threshold,
        cluster_distance=args.cluster_distance,
        min_count=args.min_count
    )


if __name__ == '__main__':
    main()
