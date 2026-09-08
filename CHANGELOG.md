# Changelog

All notable changes to NextClone (eos-jin fork) will be documented in this file.

---

## [Unreleased] — v3.0.0 (LARRY Support + CellBarcode Filtering)

### Added
- **LARRY mode** (`--mode LARRY`): Native support for Cas9-based lineage tracing data (Weinreb et al., Science 2020)
  - Extracts LARRY barcodes from paired R1/R2 FASTQ files
  - Validates barcode structure using conserved position patterns
  - Filters by minimum reads per (cell, umi, barcode) tuple
  - Clusters barcodes by Hamming distance for error correction
  - Counts and filters UMIs per (cell, barcode) combination
  - Outputs clone assignments per cell
- **CellBarcode filtering** for discovery mode (Sun et al., 2024, Nature Computational Science)
  - Replaces flexiplex knee-plot filtering
  - Implements multiple filtering strategies: auto (k-means), manual, cluster, combined
  - Configurable parameters for threshold, cluster distance, and minimum count
- **New files**:
  - `bin/larry_extract_barcodes.py`: LARRY barcode extraction from R1/R2 pairs
  - `bin/larry_filter_and_cluster.py`: LARRY filtering and clustering pipeline
  - `bin/cellbarcode_filter.R`: CellBarcode filtering using R package for discovery mode
  - `modules/extract_larry_barcodes.nf`: Nextflow processes for LARRY extraction
  - `modules/larry_filter_and_cluster.nf`: Nextflow process for LARRY filtering
  - `modules/cellbarcode_filter.nf`: Nextflow process for CellBarcode filtering
  - `LARRY_INTEGRATION.md`: Comprehensive LARRY documentation

### Changed
- Discovery mode now uses **CellBarcode filtering** instead of flexiplex knee-plot
- Removed `discovery_filter_method` parameter (CellBarcode is now the only method)
- Updated README.md with LARRY mode documentation and CellBarcode parameters
- Updated parameter tables to include CellBarcode and LARRY parameters

### LARRY Workflow Details
LARRY mode does NOT use discovery/whitelist mode. It:
1. Extracts all barcodes with the fixed LARRY prefix (`GTTGCTAGGAGAGACCATATG`)
2. Filters by minimum read count per tuple (default: 10)
3. Clusters similar barcodes by Hamming distance (default: 3)
4. Counts unique UMIs per (cell, barcode)
5. Filters by minimum UMI count (default: 3)
6. Outputs clone assignments

### CellBarcode Filtering Details
CellBarcode implements multiple filtering strategies for discovery mode:
- **auto** (default): 1D k-means clustering on log-transformed counts
- **manual**: User-specified threshold
- **cluster**: Remove barcodes similar to more abundant ones
- **combined**: Auto threshold + cluster filtering

---

## [v2.0.0] — Discovery Mode + Reports

### Added
- **Discovery mode** (`--discovery_mode true`): Two-pass barcode identification without a whitelist reference
  - Pass 1: Flexiplex discovers barcodes from raw data (`-f 0` strict matching)
  - Pass 2: Flexiplex maps reads against discovered barcodes
- **`filter_discovered_barcodes` parameter**: Optional knee-plot inflection filtering for discovered barcodes (default: `false`, keeping all barcodes including singletons)
- **Auto-generated HTML QC reports** at end of every run:
  - Clone overlap table across samples at multiple thresholds (≥5, 10, 15, 20, 50, 100 cells)
  - Heterogeneity metrics: Gini coefficient and Shannon index
  - Clone size density plot (KDE curve)
  - Ranked clone abundance (log scale)
  - Edit distance QC (FlankEditDist & BarcodeEditDist)
  - Cross-sample clonality comparison
- **Comparison report generator** (`reports/generate_comparison_report.py`): Compare two runs side-by-side
- **Run log generator** (`generate_run_log` process): Records all parameters, software versions, and git commit for reproducibility
- **`all_barcodes.txt` output**: All discovered barcodes with counts (no filtering, for QC/debugging)
- **Synthetic test data** and 407-line test suite for discovery mode

### Changed
- Merged `sc_merge_discovered_barcodes_nofilter` into single `sc_merge_discovered_barcodes` process (was duplicated code)
- `sc_merge_discovered_barcodes` now outputs two channels: `all_barcodes` and `filtered_barcodes`
- Gini/Shannon metrics rounded to 2 decimal places
- Added barcode header comments to `all_barcodes.txt`
- Cross-sample charts sorted alphabetically
- Default `filter_discovered_barcodes = false` (recommended for lineage tracing)

### Fixed
- **Critical**: Don't call flexiplex-filter when filtering is disabled (was causing pipeline crashes)
- **Critical**: Use `cp` instead of `cat` for filtered_barcodes.txt when filtering disabled
- **Critical**: Two discovery mode pipeline bugs (detailed in commit `165b480`)
- Escape all bash `$` variables in Nextflow template strings (was causing interpolation bugs)
- Added validation for `combined_barcodes_counts.txt` + debug output
- Validation for `clone_barcodes_reference` requirement when `discovery_mode = false`
- Charts blank on load → defer auto-select to `window.load` event
- Remove duplicate 'Sample: xxx' heading below dropdown
- Remove average metrics, fix run mode display, fix density chart

### Removed
- Removed `defaults` channel from conda environments (WEHI HPC Anaconda policy)
- Removed `tenx_whitelist` parameter (simplified discovery mode)
- Removed dead `sc_filter_discovered_barcodes` process
- Removed backup file `generate_report.py.bak`

### Infrastructure
- Enabled `useMamba = true` for faster conda environment creation
- Added `generate_report` and `generate_run_log` processes to scRNAseq workflow
- Added `reports/README.md` with usage instructions
- Updated README.md with discovery mode documentation, parameter table, and output file list

---

## [v1.0.0] — Initial Release (upstream)

Original release by Phipson Lab. Whitelist mode only. See [upstream repo](https://github.com/phipsonlab/NextClone).
