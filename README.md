# NextClone [![DOI](https://zenodo.org/badge/679986518.svg)](https://doi.org/10.5281/zenodo.15686558)

NextClone is a Nextflow pipeline to facilitate rapid extraction and quantification 
of clonal barcodes from both DNA-seq and scRNAseq data.
DNA-seq data refers to dedicated DNA barcoding data which exclusively sequences 
the synthetic lineage tracing clone barcode reads using Next Generation Sequencing.

<p> <img src="Nextclone_diagram_v5.png" width="500"/> </p>

The pipeline comprises two distinct workflows, one for DNA-seq data and the other for scRNAseq data. 
Both workflows are highly modular and adaptable, with software that can easily be substituted as required, 
and with parameters that can be tailored through the nextflow.config file to suit diverse needs.
It is heavily optimised for usage in high-performance computing (HPC) platforms.

## Documentation

For instructions on how to use *NextClone*, please visit the [user guide](https://phipsonlab.github.io/NextClone/).

## Discovery Mode

NextClone now supports **discovery mode**, which enables barcode identification without requiring a pre-defined whitelist of known barcodes. This is particularly useful when:

- The exact barcode sequences are unknown
- You want to discover novel barcodes from your data
- You're working with a new clonal barcoding system

### How Discovery Mode Works

Discovery mode uses a two-pass approach powered by [Flexiplex](https://github.com/DavidsonGroup/flexiplex):

1. **Pass 1 (Discovery):** Run Flexiplex without a known barcode list (`-k` flag) to identify all potential barcodes in the data. Uses strict flanking sequence matching (`-f 0`) to reduce barcode errors.

2. **Filtering:** Use `flexiplex-filter` to identify high-quality barcodes using the knee-plot inflection point method. Optionally, discovered barcodes can be intersected with a 10x barcode whitelist.

3. **Pass 2 (Mapping):** Run Flexiplex with the filtered barcode list to perform final read assignments with standard edit distance parameters.

### Usage

Enable discovery mode by setting the `discovery_mode` parameter:

```bash
nextflow run main.nf --discovery_mode true
```

Optionally, provide a 10x barcode whitelist to filter discovered barcodes:

```bash

```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `discovery_mode` | `false` | Enable two-pass barcode discovery mode |
| `filter_discovered_barcodes` | `false` | Apply knee-plot filtering to discovered barcodes (see below) |

When `discovery_mode = false` (default), the pipeline requires `clone_barcodes_reference` to be provided with a list of known barcodes.

## HTML Reports

### Standard report (auto-generated)

NextClone automatically generates an interactive HTML dashboard at the end of every run. The report is saved to your `publish_dir` as `nextclone_report.html`.

The report includes:
- Sample overview table (reads, cells, unique clones, clonality)
- Ranked clone abundance plot (log scale)
- Clone size distribution (singleton → dominant)
- Top 20 clones per sample
- Edit distance QC (FlankEditDist & BarcodeEditDist)
- Cross-sample clonality comparison

To customise the report title:
```bash
nextflow run main.nf --report_title "My Experiment — ZR751 2026"
```

### Comparison report (manual, two runs)

To compare two runs (e.g. reference mode vs discovery mode), use the comparison script after both runs are complete:

```bash
python3 reports/generate_comparison_report.py \
    /path/to/run_a/clone_barcodes.csv \
    /path/to/run_b/clone_barcodes.csv \
    --label-a "Reference" \
    --label-b "Discovery" \
    --output comparison_report.html \
    --title "Reference vs Discovery — My Experiment"
```

The comparison report shows:
- Δ reads, cells, and clones between the two runs
- Per-sample ranked abundance overlay (both modes on one log-scale plot)
- Clone size distribution side by side
- Top clone overlap (concordance between modes)
- Clonality metrics comparison (top1%, top3%, top10%)
- Cell recovery validation across samples

> **No pip installs required.** Both scripts use Python stdlib only, with Chart.js loaded via CDN for charts.

---

### Barcode filtering in discovery mode

By default (`filter_discovered_barcodes = false`), **all barcodes discovered in Pass 1 are passed to Pass 2**, including singletons. This is the recommended setting for lineage tracing experiments where rare clones are biologically meaningful and should not be discarded.

Setting `filter_discovered_barcodes = true` enables knee-plot inflection filtering via `flexiplex-filter`, which removes low-count barcodes. This can be useful for noisy datasets but **will discard singleton and low-count clones** that may be genuine:

```bash
nextflow run main.nf --discovery_mode true --filter_discovered_barcodes true
```

<!-- ## Citation -->

<!-- If you use NextClone in your study, please kindly cite our preprint on bioRxiv. -->
