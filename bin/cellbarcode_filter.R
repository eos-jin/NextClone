#!/usr/bin/env Rscript
# CellBarcode-based filtering for NextClone discovery mode
# Uses the official CellBarcode R package for barcode filtering
#
# Reference: Sun et al. (2024) Nature Computational Science
# Package: https://github.com/wenjie1991/CellBarcode

suppressPackageStartupMessages({
  library(CellBarcode)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: cellbarcode_filter.R <input_counts_file> <output_filtered_file> [method] [threshold] [cluster_distance] [min_count]")
}

input_file <- args[1]
output_file <- args[2]
method <- ifelse(length(args) >= 3, args[3], "auto")
threshold <- ifelse(length(args) >= 4, as.numeric(args[4]), NA)
cluster_distance <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
min_count <- ifelse(length(args) >= 6, as.numeric(args[6]), 1)

# Read barcode counts
message("Reading barcode counts from ", input_file)
counts_data <- fread(input_file, header = FALSE, sep = "\t")
setnames(counts_data, c("barcode_seq", "count"))

message("Read ", nrow(counts_data), " barcodes")

# Filter by minimum count first
counts_data <- counts_data[count >= min_count]
message("After min_count filter (", min_count, "): ", nrow(counts_data), " barcodes")

if (nrow(counts_data) == 0) {
  stop("No barcodes remaining after min_count filter")
}

# Create BarcodeObj from the count data
# CellBarcode expects a matrix with barcodes as rows and samples as columns
# We have a single sample, so create a matrix with one column
barcode_matrix <- as.matrix(counts_data$count)
rownames(barcode_matrix) <- counts_data$barcode_seq
colnames(barcode_matrix) <- "sample1"

# Create BarcodeObj
bc_obj <- bc_create_BarcodeObj(barcode_matrix, sample_name = "sample1")

message("Created BarcodeObj with ", length(bc_barcodes(bc_obj)), " barcodes")

# Apply filtering based on method
if (method == "auto") {
  # Automatic threshold using CellBarcode's bc_auto_cutoff
  message("Applying automatic threshold filtering (k-means)")
  
  # bc_cure_depth with depth=-1 uses automatic cutoff
  bc_obj <- bc_cure_depth(bc_obj, depth = -1)
  
  # Get the cutoff that was used
  cutoff <- bc_auto_cutoff(bc_obj, useCleanBc = FALSE)
  message("Automatic threshold: ", round(cutoff, 2))
  
} else if (method == "manual") {
  # Manual threshold
  if (is.na(threshold)) {
    stop("Manual method requires threshold parameter")
  }
  
  message("Applying manual threshold: ", threshold)
  bc_obj <- bc_cure_depth(bc_obj, depth = threshold)
  
} else if (method == "cluster") {
  # Cluster filtering using CellBarcode's bc_cure_cluster
  message("Applying cluster filtering (max distance: ", cluster_distance, ")")
  
  # First need to cure depth (use min_count as threshold)
  bc_obj <- bc_cure_depth(bc_obj, depth = min_count)
  
  # Then cluster similar barcodes
  bc_obj <- bc_cure_cluster(
    bc_obj,
    dist_threshold = cluster_distance,
    dist_method = "hamm"  # Hamming distance
  )
  
} else if (method == "combined") {
  # Combined: auto threshold + cluster filtering
  message("Applying combined filtering (auto + cluster)")
  
  # Step 1: auto threshold
  bc_obj <- bc_cure_depth(bc_obj, depth = -1)
  cutoff <- bc_auto_cutoff(bc_obj, useCleanBc = FALSE)
  message("Step 1 - Auto threshold: ", round(cutoff, 2))
  
  # Step 2: cluster filtering
  message("Step 2 - Cluster filtering (max distance: ", cluster_distance, ")")
  bc_obj <- bc_cure_cluster(
    bc_obj,
    dist_threshold = cluster_distance,
    dist_method = "hamm"
  )
  
} else {
  stop("Unknown method: ", method, ". Use 'auto', 'manual', 'cluster', or 'combined'")
}

# Extract filtered barcodes
filtered_barcodes <- bc_barcodes(bc_obj, unlist = TRUE)
message("Filtered to ", length(filtered_barcodes), " barcodes")

# Write output (just the barcode sequences, one per line)
writeLines(filtered_barcodes, output_file)
message("Output written to ", output_file)
