#!/usr/bin/env python3
"""
NextClone Comparison Report Generator
Reads two clone_barcodes.csv files and generates a self-contained HTML comparison dashboard.

Usage:
    python3 generate_comparison_report.py <csv_a> <csv_b> \
        --label-a "Reference" --label-b "Discovery (No Filter)" \
        --output report_comparison.html \
        --title "NextClone: Reference vs Discovery Mode — ZR751"
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from datetime import datetime


# ---------------------------------------------------------------------------
# Data loading & stats computation
# ---------------------------------------------------------------------------

def load_data(csv_path):
    """Parse the CSV and return a dict of per-sample data structures."""
    samples = defaultdict(lambda: {
        "reads": 0,
        "cells": set(),
        "clone_cells": defaultdict(set),  # clone_barcode -> set of cell barcodes
    })

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["SourceBAMFile"]
            cell = row["CellBarcode"]
            clone = row["CloneBarcode"]

            s = samples[sample]
            s["reads"] += 1
            s["cells"].add(cell)
            s["clone_cells"][clone].add(cell)

    return dict(samples)


def compute_sample_stats(sample_data):
    """Compute derived stats for a single sample dict."""
    clone_cells = sample_data["clone_cells"]
    total_cells = len(sample_data["cells"])

    # Clone sizes: number of unique cells per clone
    clone_sizes = {clone: len(cells) for clone, cells in clone_cells.items()}
    sorted_clones = sorted(clone_sizes.items(), key=lambda x: -x[1])

    n_clones = len(sorted_clones)
    n_cells = total_cells

    # Top1, Top3, Top10 %
    def pct_top_n(n):
        if n_cells == 0:
            return 0.0
        top_cells = sum(sz for _, sz in sorted_clones[:n])
        return round(100.0 * top_cells / n_cells, 2)

    top1_pct = pct_top_n(1)
    top3_pct = pct_top_n(3)
    top10_pct = pct_top_n(10)

    # Ranked sizes for top 100 (log abundance plot)
    ranked_sizes = [sz for _, sz in sorted_clones[:100]]

    # Size buckets
    buckets = {"singleton": 0, "small": 0, "medium": 0, "large": 0, "dominant": 0}
    for _, sz in sorted_clones:
        if sz == 1:
            buckets["singleton"] += 1
        elif sz <= 5:
            buckets["small"] += 1
        elif sz <= 20:
            buckets["medium"] += 1
        elif sz <= 100:
            buckets["large"] += 1
        else:
            buckets["dominant"] += 1

    # Top 20 clones
    top20 = []
    for clone, sz in sorted_clones[:20]:
        pct = round(100.0 * sz / n_cells, 2) if n_cells > 0 else 0.0
        top20.append({
            "barcode": clone[:12] + "…" if len(clone) > 12 else clone,
            "barcode_full": clone,
            "n_cells": sz,
            "pct": pct,
        })

    return {
        "reads": sample_data["reads"],
        "cells": n_cells,
        "clones": n_clones,
        "top1_pct": top1_pct,
        "top3_pct": top3_pct,
        "top10_pct": top10_pct,
        "ranked_sizes": ranked_sizes,
        "buckets": buckets,
        "top20": top20,
        "clone_sizes": clone_sizes,  # full dict for cross-run lookup
    }


def build_comparison_data(data_a, data_b, label_a, label_b):
    """Build the full comparison dataset for the HTML template."""
    samples_a = {name: compute_sample_stats(sd) for name, sd in data_a.items()}
    samples_b = {name: compute_sample_stats(sd) for name, sd in data_b.items()}

    all_samples = sorted(set(list(samples_a.keys()) + list(samples_b.keys())))

    # Per-sample comparison rows
    sample_rows = []
    for sample in all_samples:
        sa = samples_a.get(sample)
        sb = samples_b.get(sample)

        def delta_pct(a, b):
            if a is None or b is None or a == 0:
                return None
            return round(100.0 * (b - a) / a, 1)

        row = {
            "sample": sample,
            "reads_a": sa["reads"] if sa else 0,
            "reads_b": sb["reads"] if sb else 0,
            "delta_reads": delta_pct(sa["reads"] if sa else None, sb["reads"] if sb else None),
            "cells_a": sa["cells"] if sa else 0,
            "cells_b": sb["cells"] if sb else 0,
            "delta_cells": delta_pct(sa["cells"] if sa else None, sb["cells"] if sb else None),
            "clones_a": sa["clones"] if sa else 0,
            "clones_b": sb["clones"] if sb else 0,
            "delta_clones": delta_pct(sa["clones"] if sa else None, sb["clones"] if sb else None),
        }
        sample_rows.append(row)

    # Per-sample detail data for charts
    sample_detail = {}
    for sample in all_samples:
        sa = samples_a.get(sample, {})
        sb = samples_b.get(sample, {})

        # Top clones overlap: top 10 from A, look up in B
        top10_clones_a = sa.get("top20", [])[:10]
        clone_sizes_b = sb.get("clone_sizes", {})

        overlap = []
        for cl in top10_clones_a:
            full_bc = cl["barcode_full"]
            cells_b = clone_sizes_b.get(full_bc, 0)
            overlap.append({
                "label": cl["barcode"],
                "cells_a": cl["n_cells"],
                "cells_b": cells_b,
            })

        sample_detail[sample] = {
            "ranked_a": sa.get("ranked_sizes", []),
            "ranked_b": sb.get("ranked_sizes", []),
            "clones_a": sa.get("clones", 0),
            "clones_b": sb.get("clones", 0),
            "buckets_a": sa.get("buckets", {}),
            "buckets_b": sb.get("buckets", {}),
            "overlap": overlap,
            "top1_a": sa.get("top1_pct", 0),
            "top3_a": sa.get("top3_pct", 0),
            "top10_a": sa.get("top10_pct", 0),
            "top1_b": sb.get("top1_pct", 0),
            "top3_b": sb.get("top3_pct", 0),
            "top10_b": sb.get("top10_pct", 0),
        }

    # Global summary totals
    total_reads_a = sum(s["reads"] for s in samples_a.values())
    total_reads_b = sum(s["reads"] for s in samples_b.values())
    total_cells_a = sum(s["cells"] for s in samples_a.values())
    total_cells_b = sum(s["cells"] for s in samples_b.values())
    total_clones_a = sum(s["clones"] for s in samples_a.values())
    total_clones_b = sum(s["clones"] for s in samples_b.values())

    def fmt_delta(a, b):
        if a == 0:
            return "N/A"
        d = 100.0 * (b - a) / a
        sign = "+" if d > 0 else ""
        return f"{sign}{d:.1f}%"

    summary = {
        "total_reads_a": total_reads_a,
        "total_reads_b": total_reads_b,
        "delta_reads": fmt_delta(total_reads_a, total_reads_b),
        "total_cells_a": total_cells_a,
        "total_cells_b": total_cells_b,
        "delta_cells": fmt_delta(total_cells_a, total_cells_b),
        "total_clones_a": total_clones_a,
        "total_clones_b": total_clones_b,
        "delta_clones": fmt_delta(total_clones_a, total_clones_b),
        "samples_a": len(samples_a),
        "samples_b": len(samples_b),
    }

    # Section 3 cross-sample chart data (sorted by clones_a desc)
    cross_sorted = sorted(sample_rows, key=lambda r: -r["clones_a"])

    return {
        "summary": summary,
        "sample_rows": sample_rows,
        "sample_detail": sample_detail,
        "all_samples": all_samples,
        "cross_sorted": cross_sorted,
        "label_a": label_a,
        "label_b": label_b,
    }


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>{title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet"/>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{
  font-family: 'Inter', sans-serif;
  background: #F1F5F9;
  color: #1E293B;
  font-size: 14px;
  line-height: 1.5;
}}
a {{ color: inherit; text-decoration: none; }}

/* ── Layout ── */
.container {{ max-width: 1280px; margin: 0 auto; padding: 24px 20px; }}
.card {{
  background: #fff;
  border-radius: 12px;
  box-shadow: 0 1px 4px rgba(0,0,0,.06), 0 4px 16px rgba(0,0,0,.04);
  padding: 24px;
  margin-bottom: 24px;
}}

/* ── Header ── */
.header {{
  background: linear-gradient(135deg, #1E3A5F 0%, #1E40AF 100%);
  color: #fff;
  border-radius: 12px;
  padding: 28px 28px 20px;
  margin-bottom: 24px;
}}
.header h1 {{ font-size: 22px; font-weight: 700; margin-bottom: 6px; }}
.header .meta {{ font-size: 12px; opacity: .75; margin-bottom: 14px; }}
.badges {{ display: flex; gap: 10px; margin-bottom: 18px; flex-wrap: wrap; }}
.badge {{
  display: inline-flex; align-items: center; gap: 6px;
  padding: 5px 12px; border-radius: 99px; font-size: 12px; font-weight: 600;
}}
.badge-a {{ background: #2563EB; color: #fff; }}
.badge-b {{ background: #16A34A; color: #fff; }}

/* ── Summary bar ── */
.summary-bar {{
  display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px;
}}
@media (max-width: 700px) {{ .summary-bar {{ grid-template-columns: repeat(2, 1fr); }} }}
.summary-metric {{
  background: rgba(255,255,255,.12);
  border-radius: 10px;
  padding: 12px 14px;
  border: 1px solid rgba(255,255,255,.18);
}}
.summary-metric .label {{ font-size: 11px; opacity: .8; text-transform: uppercase; letter-spacing: .04em; margin-bottom: 4px; }}
.summary-metric .vals {{ font-size: 13px; font-weight: 600; margin-bottom: 2px; }}
.summary-metric .delta {{ font-size: 12px; font-weight: 500; }}
.delta-red {{ color: #FCA5A5; }}
.delta-green {{ color: #86EFAC; }}
.delta-gray {{ color: rgba(255,255,255,.6); }}

/* ── Section title ── */
.section-title {{
  font-size: 16px; font-weight: 700; color: #1E293B;
  margin-bottom: 16px; padding-bottom: 10px;
  border-bottom: 2px solid #E2E8F0;
  display: flex; align-items: center; gap: 8px;
}}
.section-num {{
  background: #2563EB; color: #fff;
  width: 24px; height: 24px; border-radius: 50%;
  display: inline-flex; align-items: center; justify-content: center;
  font-size: 12px; font-weight: 700; flex-shrink: 0;
}}

/* ── Table ── */
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
thead th {{
  background: #F8FAFC; color: #64748B; font-weight: 600;
  text-transform: uppercase; font-size: 11px; letter-spacing: .04em;
  padding: 10px 12px; text-align: left; border-bottom: 2px solid #E2E8F0;
}}
tbody tr {{ cursor: pointer; transition: background .15s; }}
tbody tr:hover {{ background: #EFF6FF; }}
tbody tr.selected {{ background: #DBEAFE; }}
tbody td {{ padding: 10px 12px; border-bottom: 1px solid #F1F5F9; }}
.num {{ text-align: right; font-variant-numeric: tabular-nums; }}

/* ── Delta pills ── */
.pill {{
  display: inline-block; padding: 2px 8px; border-radius: 99px;
  font-size: 11px; font-weight: 600; white-space: nowrap;
}}
.pill-red {{ background: #FEE2E2; color: #DC2626; }}
.pill-green {{ background: #DCFCE7; color: #16A34A; }}
.pill-gray {{ background: #F1F5F9; color: #64748B; }}
.pill-large-red {{ background: #FEE2E2; color: #9B1C1C; font-size: 12px; font-weight: 700; }}

/* ── Per-sample detail ── */
#sample-detail {{ display: none; }}
#sample-detail.visible {{ display: block; }}
.sample-heading {{
  font-size: 18px; font-weight: 700; color: #1E293B; margin-bottom: 20px;
  display: flex; align-items: center; gap: 10px;
}}
.sample-heading-tag {{
  background: #EFF6FF; color: #2563EB;
  border-radius: 8px; padding: 4px 12px;
  font-size: 14px; font-weight: 600;
}}
.chart-grid {{
  display: grid; grid-template-columns: 1fr 1fr; gap: 20px;
}}
@media (max-width: 900px) {{ .chart-grid {{ grid-template-columns: 1fr; }} }}
.chart-card {{
  background: #FAFAFA; border-radius: 10px;
  border: 1px solid #E2E8F0; padding: 16px;
}}
.chart-card h4 {{ font-size: 13px; font-weight: 600; color: #475569; margin-bottom: 12px; }}
.chart-wrap {{ position: relative; height: 260px; }}

/* ── Section 3 ── */
.cross-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
@media (max-width: 900px) {{ .cross-grid {{ grid-template-columns: 1fr; }} }}
.cross-note {{
  font-size: 12px; color: #64748B; background: #F8FAFC;
  border-left: 3px solid #2563EB; padding: 8px 12px;
  border-radius: 0 6px 6px 0; margin-top: 10px; line-height: 1.6;
}}

/* ── Legend chips ── */
.legend {{ display: flex; gap: 14px; margin-bottom: 10px; flex-wrap: wrap; }}
.legend-item {{ display: flex; align-items: center; gap: 5px; font-size: 12px; color: #475569; }}
.legend-dot {{ width: 12px; height: 12px; border-radius: 3px; flex-shrink: 0; }}
</style>
</head>
<body>
<div class="container">

<!-- ── HEADER ── -->
<div class="header">
  <h1>{title}</h1>
  <div class="meta">
    Run A: <strong>{file_a}</strong> &nbsp;·&nbsp;
    Run B: <strong>{file_b}</strong> &nbsp;·&nbsp;
    Generated: <strong>{date}</strong>
  </div>
  <div class="badges">
    <span class="badge badge-a">● {label_a}</span>
    <span class="badge badge-b">● {label_b}</span>
  </div>
  <div class="summary-bar" id="summary-bar"></div>
</div>

<!-- ── SECTION 1: Overview Table ── -->
<div class="card">
  <div class="section-title"><span class="section-num">1</span> Sample Overview Comparison</div>
  <div class="legend">
    <div class="legend-item"><div class="legend-dot" style="background:#2563EB"></div> {label_a}</div>
    <div class="legend-item"><div class="legend-dot" style="background:#16A34A"></div> {label_b}</div>
  </div>
  <table id="overview-table">
    <thead>
      <tr>
        <th>Sample</th>
        <th class="num">Reads A</th><th class="num">Reads B</th><th class="num">Δ Reads</th>
        <th class="num">Cells A</th><th class="num">Cells B</th><th class="num">Δ Cells</th>
        <th class="num">Clones A</th><th class="num">Clones B</th><th class="num">Δ Clones</th>
      </tr>
    </thead>
    <tbody id="overview-tbody"></tbody>
  </table>
  <p style="font-size:12px;color:#94A3B8;margin-top:10px;">Click a row or use the dropdown below to view per-sample detail.</p>
</div>

<!-- ── SECTION 2: Per-sample detail ── -->
<div class="card" id="sample-detail">
  <div class="section-title"><span class="section-num">2</span> Per-Sample Detail</div>
  <div style="margin-bottom:16px;">
    <label for="sample-selector" style="font-size:13px;font-weight:600;color:#475569;margin-right:10px;">Sample:</label>
    <select id="sample-selector" onchange="onSelectorChange(this.value)" style="font-family:inherit;font-size:13px;padding:6px 12px;border:1px solid #CBD5E1;border-radius:6px;background:#fff;color:#1E293B;cursor:pointer;">
      <option value="">— select a sample —</option>
    </select>
  </div>
  <span id="selected-sample-name" style="display:none;"></span>
  <div class="legend" style="margin-bottom:16px;">
    <div class="legend-item"><div class="legend-dot" style="background:#2563EB"></div> {label_a}</div>
    <div class="legend-item"><div class="legend-dot" style="background:#16A34A"></div> {label_b}</div>
  </div>
  <div class="chart-grid">
    <div class="chart-card">
      <h4>A) Ranked Clone Abundance (log scale)</h4>
      <div class="chart-wrap"><canvas id="chart-abundance"></canvas></div>
    </div>
    <div class="chart-card">
      <h4>B) Clone Size Distribution</h4>
      <div class="chart-wrap"><canvas id="chart-size-dist"></canvas></div>
    </div>
    <div class="chart-card">
      <h4>C) Top 10 Clones Overlap</h4>
      <div class="chart-wrap"><canvas id="chart-overlap"></canvas></div>
    </div>
    <div class="chart-card">
      <h4>D) Clonality Metrics</h4>
      <div class="chart-wrap"><canvas id="chart-clonality"></canvas></div>
    </div>
  </div>
</div>

<!-- ── SECTION 3: Cross-sample summary ── -->
<div class="card">
  <div class="section-title"><span class="section-num">3</span> Cross-Sample Summary</div>
  <div class="legend" style="margin-bottom:16px;">
    <div class="legend-item"><div class="legend-dot" style="background:#2563EB"></div> {label_a}</div>
    <div class="legend-item"><div class="legend-dot" style="background:#16A34A"></div> {label_b}</div>
  </div>
  <div class="cross-grid">
    <div>
      <div style="font-size:13px;font-weight:600;color:#475569;margin-bottom:10px;">E) Clone Count per Sample</div>
      <div style="position:relative;height:340px;"><canvas id="chart-cross-clones"></canvas></div>
      <div class="cross-note">
        Reference mode uses a complete barcode library whitelist; Discovery mode identifies barcodes de novo from data above a detection threshold.
      </div>
    </div>
    <div>
      <div style="font-size:13px;font-weight:600;color:#475569;margin-bottom:10px;">F) Cell Recovery per Sample</div>
      <div style="position:relative;height:340px;"><canvas id="chart-cross-cells"></canvas></div>
      <div class="cross-note">
        Cell counts are similar across modes (~90% recovery), validating that the core clonal architecture is preserved even in discovery mode.
      </div>
    </div>
  </div>
</div>

</div><!-- /container -->

<script>
// ── Injected data ──
const DATA = {data_json};
const LABEL_A = DATA.label_a;
const LABEL_B = DATA.label_b;
const COLOR_A = '#2563EB';
const COLOR_B = '#16A34A';
const COLOR_A_ALPHA = 'rgba(37,99,235,0.15)';
const COLOR_B_ALPHA = 'rgba(22,163,74,0.15)';

// ── Utility ──
function fmt(n) {{
  if (n === null || n === undefined) return '—';
  return n.toLocaleString();
}}

function deltaClass(d) {{
  if (d === null) return 'gray';
  if (d <= -5) return 'red';
  if (d >= 5) return 'green';
  return 'gray';
}}

function deltaPill(d, large) {{
  if (d === null) return '<span class="pill pill-gray">—</span>';
  const sign = d > 0 ? '+' : '';
  const cls = d <= -5 ? (large ? 'pill-large-red' : 'pill-red') :
              d >= 5  ? 'pill-green' : 'pill-gray';
  return `<span class="pill ${{cls}}">${{sign}}${{d}}%</span>`;
}}

// ── Summary bar ──
(function() {{
  const s = DATA.summary;
  function metricHTML(label, va, vb, delta) {{
    const dc = delta.startsWith('+') ? 'delta-green' :
               delta.startsWith('-') ? 'delta-red' : 'delta-gray';
    return `<div class="summary-metric">
      <div class="label">${{label}}</div>
      <div class="vals">${{va}} vs ${{vb}}</div>
      <div class="delta ${{dc}}">${{delta}}</div>
    </div>`;
  }}
  const bar = document.getElementById('summary-bar');
  bar.innerHTML =
    metricHTML('Total Reads', fmt(s.total_reads_a), fmt(s.total_reads_b), s.delta_reads) +
    metricHTML('Total Cells', fmt(s.total_cells_a), fmt(s.total_cells_b), s.delta_cells) +
    metricHTML('Total Clones', fmt(s.total_clones_a), fmt(s.total_clones_b), s.delta_clones) +
    metricHTML('Samples', s.samples_a, s.samples_b, s.samples_a === s.samples_b ? '=' : '≠');
}})();

// ── Populate sample selector dropdown ──
(function() {{
  const sel = document.getElementById('sample-selector');
  DATA.sample_rows.forEach(row => {{
    const opt = document.createElement('option');
    opt.value = row.sample;
    opt.textContent = row.sample;
    sel.appendChild(opt);
  }});
}})();

function onSelectorChange(sample) {{
  if (!sample) return;
  const tr = document.querySelector(`#overview-tbody tr[data-sample="${{sample}}"]`);
  selectSample(sample, tr);
}}

// ── Overview table ──
(function() {{
  const tbody = document.getElementById('overview-tbody');
  DATA.sample_rows.forEach((row, i) => {{
    const tr = document.createElement('tr');
    tr.dataset.sample = row.sample;
    tr.innerHTML = `
      <td><strong>${{row.sample}}</strong></td>
      <td class="num">${{fmt(row.reads_a)}}</td>
      <td class="num">${{fmt(row.reads_b)}}</td>
      <td class="num">${{deltaPill(row.delta_reads, false)}}</td>
      <td class="num">${{fmt(row.cells_a)}}</td>
      <td class="num">${{fmt(row.cells_b)}}</td>
      <td class="num">${{deltaPill(row.delta_cells, false)}}</td>
      <td class="num">${{fmt(row.clones_a)}}</td>
      <td class="num">${{fmt(row.clones_b)}}</td>
      <td class="num">${{deltaPill(row.delta_clones, true)}}</td>
    `;
    tr.addEventListener('click', () => {{
      document.getElementById('sample-selector').value = row.sample;
      selectSample(row.sample, tr);
    }});
    tbody.appendChild(tr);
  }});
}})();

// Auto-select first sample after page fully loads (ensures Chart.js is ready
// and canvases have non-zero dimensions)
window.addEventListener('load', function() {{
  if (DATA.sample_rows.length > 0) {{
    const first = DATA.sample_rows[0];
    const sel = document.getElementById('sample-selector');
    sel.value = first.sample;
    const firstTr = document.querySelector('#overview-tbody tr');
    // Select without scrolling on initial load
    if (firstTr) firstTr.classList.add('selected');
    document.getElementById('sample-selector').value = first.sample;
    document.getElementById('selected-sample-name').textContent = first.sample;
    
    document.getElementById('sample-detail').classList.add('visible');
    renderSampleCharts(first.sample);
  }}
}});

// ── Chart instances ──
let chartAbundance = null, chartSizeDist = null, chartOverlap = null, chartClonality = null;

function destroyCharts() {{
  [chartAbundance, chartSizeDist, chartOverlap, chartClonality].forEach(c => {{ if (c) c.destroy(); }});
}}

function selectSample(sample, tr, scroll=true) {{
  // Highlight row
  document.querySelectorAll('#overview-tbody tr').forEach(r => r.classList.remove('selected'));
  if (tr) tr.classList.add('selected');

  // Sync dropdown
  document.getElementById('sample-selector').value = sample;

  // Show sample name heading
  document.getElementById('selected-sample-name').textContent = sample;
  

  const detail = document.getElementById('sample-detail');
  detail.classList.add('visible');
  if (scroll) detail.scrollIntoView({{ behavior: 'smooth', block: 'start' }});

  destroyCharts();
  renderSampleCharts(sample);
}}

function renderSampleCharts(sample) {{
  const d = DATA.sample_detail[sample];
  if (!d) return;

  // ── A) Ranked Clone Abundance ──
  const ranksA = d.ranked_a.map((v, i) => ({{ x: i+1, y: v }}));
  const ranksB = d.ranked_b.map((v, i) => ({{ x: i+1, y: v }}));
  chartAbundance = new Chart(document.getElementById('chart-abundance'), {{
    type: 'line',
    data: {{
      datasets: [
        {{
          label: `${{LABEL_A}} (${{d.clones_a}} clones)`,
          data: ranksA,
          borderColor: COLOR_A,
          backgroundColor: COLOR_A_ALPHA,
          pointRadius: 0,
          borderWidth: 2.5,
          tension: 0.1,
          fill: false,
        }},
        {{
          label: `${{LABEL_B}} (${{d.clones_b}} clones)`,
          data: ranksB,
          borderColor: COLOR_B,
          backgroundColor: COLOR_B_ALPHA,
          pointRadius: 0,
          borderWidth: 2.5,
          tension: 0.1,
          fill: false,
        }},
      ]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        x: {{ type: 'linear', title: {{ display: true, text: 'Clone Rank', font: {{ size: 11 }} }} }},
        y: {{
          type: 'logarithmic',
          title: {{ display: true, text: 'Cells (log)', font: {{ size: 11 }} }},
          ticks: {{
            callback: v => Number.isInteger(Math.log10(v)) ? v : ''
          }}
        }}
      }},
      plugins: {{
        legend: {{ position: 'top', labels: {{ font: {{ size: 11 }} }} }},
        tooltip: {{ mode: 'index', intersect: false }},
      }}
    }}
  }});

  // ── B) Clone Size Distribution ──
  const bucketLabels = ['Singleton\n(=1)', 'Small\n(2-5)', 'Medium\n(6-20)', 'Large\n(21-100)', 'Dominant\n(>100)'];
  const bucketKeys = ['singleton', 'small', 'medium', 'large', 'dominant'];
  chartSizeDist = new Chart(document.getElementById('chart-size-dist'), {{
    type: 'bar',
    data: {{
      labels: ['Singleton (=1)', 'Small (2-5)', 'Medium (6-20)', 'Large (21-100)', 'Dominant (>100)'],
      datasets: [
        {{
          label: LABEL_A,
          data: bucketKeys.map(k => d.buckets_a[k] || 0),
          backgroundColor: COLOR_A,
          borderRadius: 4,
        }},
        {{
          label: LABEL_B,
          data: bucketKeys.map(k => d.buckets_b[k] || 0),
          backgroundColor: COLOR_B,
          borderRadius: 4,
        }},
      ]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ position: 'top', labels: {{ font: {{ size: 11 }} }} }} }},
      scales: {{
        x: {{ ticks: {{ font: {{ size: 10 }}, maxRotation: 0 }} }},
        y: {{ title: {{ display: true, text: '# Clones', font: {{ size: 11 }} }} }}
      }}
    }}
  }});

  // ── C) Top Clones Overlap ──
  const overlapLabels = d.overlap.map(o => o.label);
  chartOverlap = new Chart(document.getElementById('chart-overlap'), {{
    type: 'bar',
    data: {{
      labels: overlapLabels,
      datasets: [
        {{
          label: LABEL_A,
          data: d.overlap.map(o => o.cells_a),
          backgroundColor: COLOR_A,
          borderRadius: 4,
        }},
        {{
          label: LABEL_B,
          data: d.overlap.map(o => o.cells_b),
          backgroundColor: COLOR_B,
          borderRadius: 4,
        }},
      ]
    }},
    options: {{
      indexAxis: 'y',
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ position: 'top', labels: {{ font: {{ size: 11 }} }} }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Cells', font: {{ size: 11 }} }} }},
        y: {{ ticks: {{ font: {{ size: 10 }} }} }}
      }}
    }}
  }});

  // ── D) Clonality Metrics ──
  chartClonality = new Chart(document.getElementById('chart-clonality'), {{
    type: 'bar',
    data: {{
      labels: ['Top 1%', 'Top 3%', 'Top 10%'],
      datasets: [
        {{
          label: LABEL_A,
          data: [d.top1_a, d.top3_a, d.top10_a],
          backgroundColor: COLOR_A,
          borderRadius: 4,
        }},
        {{
          label: LABEL_B,
          data: [d.top1_b, d.top3_b, d.top10_b],
          backgroundColor: COLOR_B,
          borderRadius: 4,
        }},
      ]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ position: 'top', labels: {{ font: {{ size: 11 }} }} }} }},
      scales: {{
        x: {{ ticks: {{ font: {{ size: 11 }} }} }},
        y: {{
          title: {{ display: true, text: '% Cells in Top Clones', font: {{ size: 11 }} }},
          max: 100,
        }}
      }}
    }}
  }});
}}

// ── Section 3: Cross-sample charts ──
(function() {{
  const cs = DATA.cross_sorted;
  const labels = cs.map(r => r.sample);

  // E) Clone counts
  new Chart(document.getElementById('chart-cross-clones'), {{
    type: 'bar',
    data: {{
      labels: labels,
      datasets: [
        {{
          label: LABEL_A,
          data: cs.map(r => r.clones_a),
          backgroundColor: COLOR_A,
          borderRadius: 4,
        }},
        {{
          label: LABEL_B,
          data: cs.map(r => r.clones_b),
          backgroundColor: COLOR_B,
          borderRadius: 4,
        }},
      ]
    }},
    options: {{
      indexAxis: 'y',
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ position: 'top', labels: {{ font: {{ size: 11 }} }} }} }},
      scales: {{
        x: {{ title: {{ display: true, text: '# Clones', font: {{ size: 11 }} }} }},
        y: {{ ticks: {{ font: {{ size: 10 }} }} }}
      }}
    }}
  }});

  // F) Cell counts
  new Chart(document.getElementById('chart-cross-cells'), {{
    type: 'bar',
    data: {{
      labels: labels,
      datasets: [
        {{
          label: LABEL_A,
          data: cs.map(r => r.cells_a),
          backgroundColor: COLOR_A,
          borderRadius: 4,
        }},
        {{
          label: LABEL_B,
          data: cs.map(r => r.cells_b),
          backgroundColor: COLOR_B,
          borderRadius: 4,
        }},
      ]
    }},
    options: {{
      indexAxis: 'y',
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ position: 'top', labels: {{ font: {{ size: 11 }} }} }} }},
      scales: {{
        x: {{ title: {{ display: true, text: '# Cells', font: {{ size: 11 }} }} }},
        y: {{ ticks: {{ font: {{ size: 10 }} }} }}
      }}
    }}
  }});
}})();
</script>
</body>
</html>
"""


def generate_html(comparison, title, file_a, file_b):
    data_json = json.dumps(comparison, separators=(',', ':'))
    return HTML_TEMPLATE.format(
        title=title,
        file_a=os.path.basename(file_a),
        file_b=os.path.basename(file_b),
        label_a=comparison["label_a"],
        label_b=comparison["label_b"],
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        data_json=data_json,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="NextClone Comparison Report Generator")
    parser.add_argument("csv_a", help="Input CSV file A (e.g. reference mode)")
    parser.add_argument("csv_b", help="Input CSV file B (e.g. discovery mode)")
    parser.add_argument("--label-a", default="Run A", help="Label for run A")
    parser.add_argument("--label-b", default="Run B", help="Label for run B")
    parser.add_argument("--output", "-o", default="report_comparison.html", help="Output HTML file")
    parser.add_argument("--title", default="NextClone: Run A vs Run B", help="Report title")
    args = parser.parse_args()

    print(f"Loading {args.csv_a} …")
    data_a = load_data(args.csv_a)
    print(f"  → {len(data_a)} samples, {sum(s['reads'] for s in data_a.values()):,} reads")

    print(f"Loading {args.csv_b} …")
    data_b = load_data(args.csv_b)
    print(f"  → {len(data_b)} samples, {sum(s['reads'] for s in data_b.values()):,} reads")

    print("Computing comparison stats …")
    comparison = build_comparison_data(data_a, data_b, args.label_a, args.label_b)

    print("Generating HTML …")
    html = generate_html(comparison, args.title, args.csv_a, args.csv_b)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(html)

    size_kb = os.path.getsize(args.output) / 1024
    print(f"✓ Report written to: {args.output}  ({size_kb:.1f} KB)")


if __name__ == "__main__":
    main()
