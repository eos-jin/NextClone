# NextClone Report Generator

Self-contained Python scripts to generate interactive HTML dashboards from NextClone output. No external dependencies — pure Python stdlib + Chart.js via CDN.

## Single-run report

Generates a per-sample HTML dashboard from a single `clone_barcodes.csv`.

```bash
python3 generate_report.py clone_barcodes.csv \
  --output report.html \
  --title "My Run"
```

**Charts included:**
- Sample overview table (reads, cells, clones, clonality)
- Ranked clone abundance (log scale)
- Clone size distribution (singleton → dominant)
- Top 20 clones (horizontal bar)
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
