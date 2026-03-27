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
nextflow run main.nf --discovery_mode true --tenx_whitelist /path/to/3M-february-2018.txt
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `discovery_mode` | `false` | Enable two-pass barcode discovery mode |
| `tenx_whitelist` | `null` | Optional path to 10x barcode whitelist for filtering |

When `discovery_mode = false` (default), the pipeline requires `clone_barcodes_reference` to be provided with a list of known barcodes.

<!-- ## Citation -->

<!-- If you use NextClone in your study, please kindly cite our preprint on bioRxiv. -->
