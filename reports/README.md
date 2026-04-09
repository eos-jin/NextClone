# NextClone Report Generator

Self-contained Python scripts to generate interactive HTML dashboards from NextClone output. No external dependencies — pure Python stdlib + Chart.js via CDN.

## Single-run report (v2)

Generates a per-sample HTML dashboard from a single `clone_barcodes.csv`.

```bash
python3 generate_report.py clone_barcodes.csv \
  --output report.html \
  --title "My Run"
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
