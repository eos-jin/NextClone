# NextClone Report Generator

Self-contained Python scripts to generate interactive HTML dashboards from NextClone output. No external dependencies — pure Python stdlib + Chart.js via CDN.

## Single-run report (v2)

Generates a per-sample HTML dashboard from a single `clone_barcodes.csv`.

### Quick Start

```bash
# Basic usage (outputs report.html)
python3 generate_report.py clone_barcodes.csv

# Custom output filename and title
python3 generate_report.py clone_barcodes.csv \
  --output my_report.html \
  --title "ZR751 Clonal Analysis — 2026-04-09"
```

### From NextClone Output

After running NextClone, generate the report from your results directory:

```bash
# If NextClone output is in results_discoverymode_260331/
cd /path/to/nextclone/results_discoverymode_260331
python3 /path/to/NextClone/reports/generate_report.py clone_barcodes.csv \
  --output nextclone_qc_report.html \
  --title "Discovery Mode — ZR751"
```

### Command-Line Options

```bash
python3 generate_report.py <input_csv> [OPTIONS]

Positional:
  input_csv              Path to clone_barcodes.csv from NextClone output

Options:
  --output FILE          Output HTML file (default: report.html)
  --title TEXT           Report title (default: "NextClone Report")
  --help                 Show help message and exit

Examples:
  # Default output (report.html)
  python3 generate_report.py clone_barcodes.csv
  
  # Custom output and title
  python3 generate_report.py clone_barcodes.csv \
    --output qc_report.html \
    --title "Sample ABC — Discovery Mode"
```

**New in v2 (2026-04-09):**
- **Clone overlap table** — shows how many clones are shared across samples at different cell thresholds (≥5, 10, 15, 20, 50, 100 cells)
- **Heterogeneity metrics** — Gini coefficient and Shannon index for each sample
- **Clone size density plot** — KDE-style curve showing clone size distribution
- **Reversed top 20 clones** — largest clones now at top of chart (easier to read)

**All charts:**
- Sample overview table (reads, cells, clones, Gini, Shannon)
- Clone overlap across samples (new!)
- Heterogeneity metrics summary (new!)
- Ranked clone abundance (log scale, top 3 annotated)
- Clone size density curve (new!)
- Top 20 clones (horizontal bar, reversed, with % labels)
- Edit distance QC (FlankEditDist + BarcodeEditDist)
- Cross-sample clonality comparison

## Comparison report

Compares two runs side by side (e.g. reference mode vs discovery mode).

```bash
python3 generate_comparison_report.py run_a.csv run_b.csv \
  --label-a "Reference" \
  --label-b "Discovery" \
  --output comparison.html \
  --title "Reference vs Discovery — ZR751"
```

**Charts included:**
- Summary header with Δ metrics (reads, cells, clones)
- Per-sample delta table (click row to drill in)
- Ranked abundance overlay (both modes on one log-scale plot)
- Clone size distribution side by side
- Top clone overlap (are the same clones found in both modes?)
- Clonality metrics comparison (top1%, top3%, top10%)
- Cross-sample clone count and cell recovery comparison
