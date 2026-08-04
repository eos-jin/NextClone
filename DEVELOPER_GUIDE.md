# NextClone Developer Guide

> **Fork:** [eos-jin/NextClone](https://github.com/eos-jin/NextClone)  
> **Upstream:** [phipsonlab/NextClone](https://github.com/phipsonlab/NextClone)  
> **Last updated:** 2026-08-04

---

## What This Fork Is

The original [NextClone](https://github.com/phipsonlab/NextClone) is a Nextflow pipeline by the Phipson Lab for extracting and counting clonal barcodes from scRNA-seq and bulk DNA-seq data. It only supported **whitelist mode** — you had to provide a known barcode reference file.

This fork adds **discovery mode** (barcode identification without a reference), **auto-generated HTML QC reports**, and **WEHI HPC compatibility**. All changes are additive and backward-compatible — the original whitelist mode works exactly as before.

The fork is **20 commits ahead of upstream/main** with +3,247 lines added.

---

## Architecture

### Pipeline Structure

```
main.nf                          ← Entry point, workflow orchestration
├── modules/extract_dnaseq_barcodes.nf    ← DNA-seq processes (9 processes)
├── modules/extract_sc_clone_barcodes.nf  ← scRNA-seq processes (11 processes)
├── conda_env/                           ← Conda environments per workflow
│   ├── extract_dnaseq_env.yaml
│   └── extract_sc_env.yaml
├── reports/                             ← HTML report generators
│   ├── generate_report.py               ← Single-run dashboard (1,117 lines)
│   └── generate_comparison_report.py    ← Cross-sample comparison (881 lines)
└── tests/                               ← Test suite
    ├── test_discovery_mode.py           ← 407-line integration tests
    └── data/                            ← Synthetic test FASTQs
```

### Data Flow (scRNAseq, discovery mode)

```
Input BAM files
    │
    ▼
┌─────────────────────────────────────────┐
│ sc_get_unmapped_reads (sambamba)        │ → Pull unmapped reads from BAM
│ sc_remove_low_qual_reads (Python)       │ → Filter by Phred quality
│ sc_retain_reads_with_CB_tag (sambamba)  │ → Keep reads with CB+UB tags
│ sc_split_unmapped_reads (Python)        │ → Split into chunks for parallelism
└─────────────────────────────────────────┘
    │
    ├─→ Pass 1: sc_discover_barcodes (flexiplex -f 0, no -k)
    │       ↓
    │   sc_merge_discovered_barcodes       → Combine counts, output:
    │       ├── all_barcodes.txt           → ALL discovered barcodes (QC)
    │       └── filtered_barcodes.txt      → Barcodes for Pass 2 mapping
    │
    └─→ Pass 2: sc_map_with_discovered_barcodes (flexiplex -k filtered_barcodes.txt)
            ↓
        sc_merge_barcodes (Python)         → Assign clones to cells
            ↓
        generate_report (Python)           → HTML QC dashboard
        generate_run_log                   → Reproducibility log
```

### Key Processes

| Process | Tool | Purpose | Label |
|---------|------|---------|-------|
| `sc_get_unmapped_reads` | sambamba | Extract unmapped reads from BAM | medium_mem |
| `sc_remove_low_qual_reads` | Python/pysam | Filter low-quality reads | small_mem |
| `sc_retain_reads_with_CB_tag` | sambamba | Keep reads with cell+UMI barcodes | medium_mem |
| `sc_split_unmapped_reads` | Python | Chunk reads for parallel processing | small_mem |
| `sc_discover_barcodes` | flexiplex | Pass 1: discover barcodes from data | varies |
| `sc_merge_discovered_barcodes` | bash+flexiplex-filter | Combine + filter discovered barcodes | small |
| `sc_map_with_discovered_barcodes` | flexiplex | Pass 2: map reads to discovered barcodes | varies |
| `sc_map_unmapped_reads` | flexiplex | Whitelist mode: map reads to known barcodes | varies |
| `sc_merge_barcodes` | Python | Collapse barcode assignments per cell | small_mem |
| `generate_report` | Python | Generate HTML QC dashboard | small |

---

## Two Modes of Operation

### Whitelist Mode (original, default)
```bash
nextflow run main.nf \
    --mode scRNAseq \
    --clone_barcodes_reference /path/to/barcodes.txt \
    --scrnaseq_bam_files /path/to/bams
```

### Discovery Mode (our addition)
```bash
# Keep all discovered barcodes (recommended for lineage tracing)
nextflow run main.nf \
    --mode scRNAseq \
    --discovery_mode true \
    --scrnaseq_bam_files /path/to/bams

# With knee-plot filtering (removes low-count barcodes)
nextflow run main.nf \
    --mode scRNAseq \
    --discovery_mode true \
    --filter_discovered_barcodes true \
    --scrnaseq_bam_files /path/to/bams
```

### How Discovery Mode Works

**Pass 1 (Discovery):** Run Flexiplex without a barcode reference (`-k` flag omitted), with strict flanking sequence matching (`-f 0`). This identifies candidate barcode sequences directly from the data.

**Merge + Filter:** Combine barcode counts from all chunks. By default, keep ALL discovered barcodes (including singletons). Optionally apply knee-plot filtering via `flexiplex-filter`.

**Pass 2 (Mapping):** Re-read the original reads and map them against the discovered barcode list. This is the same as whitelist mode but with a discovered reference instead of a provided one.

---

## What We Changed (Full List)

### 1. Discovery Mode (the big feature)
- Added `discovery_mode` parameter to `main.nf` and `nextflow.config`
- Added `filter_discovered_barcodes` parameter (knee-plot inflection filtering)
- DNAseq: `dnaseq_discover_barcodes`, `dnaseq_filter_discovered_barcodes`, `dnaseq_map_with_discovered_barcodes` processes
- scRNAseq: `sc_discover_barcodes`, `sc_merge_discovered_barcodes`, `sc_map_with_discovered_barcodes` processes
- Two-channel output from `sc_merge_discovered_barcodes`: `all_barcodes` (QC) + `filtered_barcodes` (mapping input)
- Parameter validation: error if no reference and no discovery mode; warning if both provided

### 2. HTML QC Reports
- `reports/generate_report.py` — auto-generated at end of every run (1,117 lines, pure Python stdlib)
- `reports/generate_comparison_report.py` — compare two runs (881 lines)
- Charts: clone overlap across samples, Gini/Shannon heterogeneity, clone size density, ranked abundance, edit distance QC, cross-sample comparison
- `generate_run_log` process — records all parameters, versions, and git commit for reproducibility

### 3. WEHI HPC Compatibility
- Removed `defaults` channel from conda environments (WEHI Anaconda policy)
- Enabled `useMamba = true` in `nextflow.config` for faster env creation
- Fixed all `$` variable escaping in Nextflow template strings (bash interpolation bug)

### 4. Bug Fixes
- Don't call flexiplex-filter when filtering is disabled (was causing crashes)
- Use `cp` instead of `cat` for filtered_barcodes.txt when filtering disabled
- Added validation for combined_barcodes_counts.txt
- Two critical discovery mode pipeline bugs fixed
- Gini/Shannon metrics rounded to 2 decimals, added barcode header to output
- Sort cross-sample charts alphabetically
- Charts blank on load → defer auto-select to window.load event
- Remove duplicate 'Sample: xxx' heading below dropdown

### 5. Testing
- 407-line test suite (`tests/test_discovery_mode.py`)
- Synthetic test data for both whitelist and discovery modes
- Parameter validation tests

### 6. Cleanup
- Removed dead `sc_filter_discovered_barcodes` process
- Removed backup file `generate_report.py.bak`
- Merged `sc_merge_discovered_barcodes_nofilter` into single process (was duplicated code)
- Removed `tenx_whitelist` parameter (simplified discovery mode)

---

## Known Issues & TODOs

### High Priority

1. **No Singularity/Docker profile** — currently only supports conda/mamba. For production use on HPC systems without conda, container support is needed.

2. **Flexiplex version not pinned** — the conda environments don't specify a Flexiplex version. Future Flexiplex updates could change behavior silently.

3. **Python 3.8 is EOL** — the conda environments use Python 3.8, which reached end-of-life in October 2024. Should upgrade to 3.10+.

4. **No CI/CD** — there are no automated tests running on PRs. The test suite (`test_discovery_mode.py`) requires Flexiplex to be installed, which makes it hard to run in CI.

5. **Hardcoded adapter sequences** — the default `adapter_5prime` and `adapter_3prime` in `nextflow.config` are hardcoded for a specific experimental setup. Users need to know to change these.

### Medium Priority

6. **Report generators use pure stdlib only** — this was intentional (no pip installs) but means no Plotly, no seaborn. The HTML is generated with inline JS/Chart.js. Future reports could benefit from richer visualization.

7. **No input validation on BAM files** — the pipeline doesn't check if BAM files are valid, have the expected tags (CB/UB), or contain unmapped reads before starting.

8. **`n_chunks` default is 2** — may not be optimal for large datasets. Should be documented or auto-calculated based on input size.

9. **Comparison report requires manual invocation** — not integrated into the main workflow. Could be auto-generated when multiple output directories exist.

### Low Priority

10. **No performance benchmarking** — we don't have benchmarks for how the pipeline scales with dataset size. The `long_mapping` profile (96h timeout) suggests some runs take a long time, but we don't know why.

11. **`mapping_process_profile` is confusing** — users need to choose between `regular_mapping` and `long_mapping` but it's not clear what determines this. Could be auto-selected based on chunk size.

12. **No versioned releases** — the repo has a `v1.0.0` tag but it was created early. Should tag the current state as v2.0.0 (discovery mode + reports is a major version bump).

---

## Development Workflow

### Running the test suite
```bash
# Requires flexiplex to be installed
python tests/test_discovery_mode.py
```

### Running the pipeline locally
```bash
# With test data (if you add it)
nextflow run main.nf -profile test
```

### Adding a new parameter
1. Add to `nextflow.config` under `params {}`
2. Add validation in `main.nf` workflow block
3. Update the `generate_run_log` process to include it
4. Update README.md and this developer guide

### Submitting changes upstream
This fork has significant additions that the upstream may or may not want. Before submitting a PR:
1. Rebase on latest `phipsonlab/NextClone` main
2. Ensure all tests pass
3. Consider splitting discovery mode and reports into separate PRs
4. Check that upstream hasn't already implemented similar features

---

## Who Built This

- **Original:** Phipson Lab (https://github.com/phipsonlab/NextClone)
- **This fork (discovery mode, reports, bug fixes):** Eos
