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

## Modes

### Whitelist mode (default)

Provide a list of known barcode sequences. Flexiplex maps all reads against the whitelist.

```bash
nextflow run main.nf --clone_barcodes_reference /path/to/barcodes.txt
```

### Discovery mode

NextClone supports **discovery mode**, which identifies barcodes directly from the data without a pre-defined whitelist. This is useful when:

- The exact barcode sequences are unknown
- You are working with a new or custom clonal barcoding system
- You want to validate or supplement a known barcode list

Discovery mode uses a two-pass approach powered by [Flexiplex](https://github.com/DavidsonGroup/flexiplex):

1. **Pass 1 (Discovery):** Run Flexiplex without a barcode list (`-k` flag) using strict flanking sequence matching (`-f 0`) to identify candidate barcodes.
2. **Pass 2 (Mapping):** Run Flexiplex with the discovered barcode list using standard edit distance parameters.

```bash
nextflow run main.nf --discovery_mode true
```

#### Barcode filtering in discovery mode

By default (`filter_discovered_barcodes = false`), **all barcodes discovered in Pass 1 are passed to Pass 2**, including singletons. This is recommended for lineage tracing experiments where rare clones are biologically meaningful.

Setting `filter_discovered_barcodes = true` applies `flexiplex-filter` knee-plot inflection filtering, which removes low-count barcodes. Use this only for noisy datasets — **it will discard singleton and low-count clones**:

```bash
nextflow run main.nf --discovery_mode true --filter_discovered_barcodes true
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mode` | `"scRNAseq"` | Workflow mode: `"scRNAseq"` or `"DNAseq"` |
| `clone_barcodes_reference` | — | Path to known barcode whitelist (required when `discovery_mode = false`) |
| `discovery_mode` | `false` | Enable two-pass barcode discovery mode |
| `filter_discovered_barcodes` | `false` | Apply knee-plot filtering to discovered barcodes (see above) |
| `barcode_edit_distance` | `2` | Maximum edit distance for barcode matching |
| `adapter_edit_distance` | `6` | Maximum edit distance for flanking adapter matching |
| `adapter_5prime` | — | 5′ flanking adapter sequence |
| `adapter_3prime` | — | 3′ flanking adapter sequence |
| `barcode_length` | `20` | Expected barcode length (bp) |
| `n_chunks` | `2` | Number of read chunks for parallel processing |
| `publish_dir` | `output/` | Output directory |
| `report_title` | — | Custom title for the HTML report (defaults to date-stamped title) |

## HTML Reports

### Standard report (auto-generated)

NextClone automatically generates an interactive HTML dashboard at the end of every run, saved to your `publish_dir` as `nextclone_qc_report.html`.

The report includes:
- Sample overview table (reads, cells, unique clones, clonality)
- Ranked clone abundance plot (log scale)
- Clone size distribution (singleton → dominant)
- Top 20 clones per sample
- Edit distance QC (FlankEditDist & BarcodeEditDist)
- Cross-sample clonality comparison

To set a custom title:
```bash
nextflow run main.nf --report_title "My Experiment — ZR751 2026"
```

### Comparison report (manual)

To compare two runs side by side (e.g. reference mode vs discovery mode), use the comparison script after both runs are complete:

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
- Per-sample ranked abundance overlay (both modes, log-scale)
- Clone size distribution side by side
- Top clone overlap (concordance between modes)
- Clonality metrics comparison (top1%, top3%, top10%)
- Cell recovery validation across samples

> **No pip installs required.** Both report scripts use Python stdlib only, with Chart.js loaded via CDN.

<!-- ## Citation -->

<!-- If you use NextClone in your study, please kindly cite our preprint on bioRxiv. -->
