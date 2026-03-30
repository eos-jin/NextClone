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

3. **Pass 2 (Mapping):** Run Flexiplex with the discovered barcode list to perform final read assignments with standard edit distance parameters.

### Usage

Enable discovery mode by setting the `discovery_mode` parameter:

```bash
nextflow run phipsonlab/Nextclone -r main --discovery_mode true
```

By default, discovered barcodes are filtered using a knee-plot inflection method (via `flexiplex-filter`) to remove low-confidence barcodes. If you expect a **low number of clones** in your data, disable filtering to retain all discovered barcodes:

```bash
nextflow run phipsonlab/Nextclone -r main --discovery_mode true --filter_discovered_barcodes false
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `discovery_mode` | `false` | Enable two-pass barcode discovery mode |
| `filter_discovered_barcodes` | `true` | Filter discovered barcodes using knee-plot method. Set to `false` for datasets with a low expected number of clones. |

When `discovery_mode = false` (default), the pipeline requires `clone_barcodes_reference` to be provided with a list of known barcodes.

<!-- ## Citation -->

<!-- If you use NextClone in your study, please kindly cite our preprint on bioRxiv. -->
