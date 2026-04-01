#!/usr/bin/env python3
"""
NextClone Report Generator
Reads clone_barcodes.csv and generates a self-contained HTML dashboard.

Usage:
    python3 generate_report.py <input_csv> [--output report.html] [--title "My Run"]
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
        "flank_edit": defaultdict(int),
        "barcode_edit": defaultdict(int),
    })

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["SourceBAMFile"]
            cell = row["CellBarcode"]
            clone = row["CloneBarcode"]
            try:
                fed = int(row["FlankEditDist"])
            except (ValueError, KeyError):
                fed = -1
            try:
                bed = int(row["BarcodeEditDist"])
            except (ValueError, KeyError):
                bed = -1

            s = samples[sample]
            s["reads"] += 1
            s["cells"].add(cell)
            s["clone_cells"][clone].add(cell)
            if fed >= 0:
                s["flank_edit"][min(fed, 5)] += 1
            if bed >= 0:
                s["barcode_edit"][min(bed, 5)] += 1

    return samples


def compute_stats(samples):
    """Turn raw per-sample data into serialisable stats dicts."""
    result = {}
    for sample, raw in sorted(samples.items()):
        n_reads = raw["reads"]
        n_cells = len(raw["cells"])

        # Clone sizes (by unique cells per clone)
        clone_sizes = {clone: len(cells) for clone, cells in raw["clone_cells"].items()}
        n_clones = len(clone_sizes)

        # Ranked sizes (descending)
        ranked = sorted(clone_sizes.values(), reverse=True)

        # Clone size distribution buckets
        buckets = {"Singleton": 0, "Small (2-5)": 0, "Medium (6-20)": 0,
                   "Large (21-100)": 0, "Dominant (>100)": 0}
        for sz in ranked:
            if sz == 1:
                buckets["Singleton"] += 1
            elif sz <= 5:
                buckets["Small (2-5)"] += 1
            elif sz <= 20:
                buckets["Medium (6-20)"] += 1
            elif sz <= 100:
                buckets["Large (21-100)"] += 1
            else:
                buckets["Dominant (>100)"] += 1

        # Top 20 clones
        top_clones_raw = sorted(clone_sizes.items(), key=lambda x: x[1], reverse=True)[:20]
        top_clones = [
            {
                "barcode": bc[:20],
                "n_cells": cnt,
                "pct": round(cnt / n_cells * 100, 2) if n_cells else 0,
            }
            for bc, cnt in top_clones_raw
        ]

        # Clonality metrics
        def top_n_pct(n):
            if n_cells == 0:
                return 0.0
            top_cells = sum(ranked[:n])
            return round(top_cells / n_cells * 100, 2)

        # Edit distance distributions (keys 0-5)
        def ed_dist(d):
            return [d.get(i, 0) for i in range(6)]

        result[sample] = {
            "reads": n_reads,
            "cells": n_cells,
            "clones": n_clones,
            "ranked_sizes": ranked,
            "clone_size_buckets": buckets,
            "top_clones": top_clones,
            "top1_pct": top_n_pct(1),
            "top3_pct": top_n_pct(3),
            "top10_pct": top_n_pct(10),
            "flank_edit_dist": ed_dist(raw["flank_edit"]),
            "barcode_edit_dist": ed_dist(raw["barcode_edit"]),
        }

    return result


def global_stats(stats):
    total_reads = sum(s["reads"] for s in stats.values())
    total_cells = sum(s["cells"] for s in stats.values())
    total_samples = len(stats)
    # Unique clones across all samples (count clones that appear in each sample independently)
    total_clones = sum(s["clones"] for s in stats.values())
    return {
        "total_reads": total_reads,
        "total_cells": total_cells,
        "total_samples": total_samples,
        "total_clones": total_clones,
    }


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>{{TITLE}}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet"/>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: 'Inter', system-ui, sans-serif; background: #F8FAFC; color: #1E293B; font-size: 14px; line-height: 1.6; }
  a { color: #2563EB; text-decoration: none; }
  a:hover { text-decoration: underline; }

  /* Layout */
  .container { max-width: 1400px; margin: 0 auto; padding: 0 24px; }

  /* Header */
  .header { background: linear-gradient(135deg, #1E3A5F 0%, #2563EB 100%); color: white; padding: 32px 0 28px; }
  .header h1 { font-size: 26px; font-weight: 700; letter-spacing: -0.3px; }
  .header-meta { margin-top: 8px; opacity: 0.8; font-size: 13px; display: flex; gap: 20px; flex-wrap: wrap; }
  .run-mode-badge { background: rgba(255,255,255,0.2); border-radius: 99px; padding: 2px 12px; font-size: 12px; font-weight: 500; }

  /* Summary bar */
  .summary-bar { background: white; border-bottom: 1px solid #E2E8F0; padding: 16px 0; }
  .summary-stats { display: flex; gap: 0; }
  .stat-item { flex: 1; text-align: center; padding: 8px 16px; border-right: 1px solid #E2E8F0; }
  .stat-item:last-child { border-right: none; }
  .stat-value { font-size: 28px; font-weight: 700; color: #2563EB; }
  .stat-label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; color: #64748B; margin-top: 2px; }

  /* Sections */
  .section { padding: 28px 0; }
  .section-title { font-size: 18px; font-weight: 600; color: #1E293B; margin-bottom: 16px; display: flex; align-items: center; gap: 8px; }
  .section-title::before { content: ''; display: block; width: 4px; height: 20px; background: #2563EB; border-radius: 2px; }

  /* Card */
  .card { background: white; border-radius: 12px; box-shadow: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.04); padding: 20px; }

  /* Table */
  .table-wrapper { overflow-x: auto; border-radius: 12px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
  table { width: 100%; border-collapse: collapse; background: white; }
  thead tr { background: #F8FAFC; }
  th { padding: 12px 16px; text-align: left; font-size: 12px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #64748B; border-bottom: 2px solid #E2E8F0; cursor: pointer; user-select: none; white-space: nowrap; }
  th:hover { background: #EFF6FF; color: #2563EB; }
  th .sort-icon { display: inline-block; margin-left: 4px; opacity: 0.4; }
  th.sort-asc .sort-icon::after { content: ' ↑'; opacity: 1; }
  th.sort-desc .sort-icon::after { content: ' ↓'; opacity: 1; }
  tbody tr { border-bottom: 1px solid #F1F5F9; cursor: pointer; transition: background 0.1s; }
  tbody tr:last-child { border-bottom: none; }
  tbody tr:hover { background: #EFF6FF; }
  tbody tr.selected { background: #DBEAFE; }
  td { padding: 12px 16px; }
  .num-cell { text-align: right; font-variant-numeric: tabular-nums; }

  /* Clonality pill */
  .pill { display: inline-block; padding: 2px 10px; border-radius: 99px; font-size: 12px; font-weight: 500; }
  .pill-green { background: #DCFCE7; color: #16A34A; }
  .pill-amber { background: #FEF3C7; color: #D97706; }
  .pill-red { background: #FEE2E2; color: #DC2626; }

  /* Sample detail */
  .detail-header { display: flex; align-items: center; justify-content: space-between; margin-bottom: 16px; flex-wrap: wrap; gap: 12px; }
  .detail-select { padding: 8px 12px; border: 1px solid #CBD5E1; border-radius: 8px; font-family: inherit; font-size: 14px; background: white; cursor: pointer; }
  .detail-select:focus { outline: none; border-color: #2563EB; box-shadow: 0 0 0 3px rgba(37,99,235,0.1); }
  .charts-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
  @media (max-width: 900px) { .charts-grid { grid-template-columns: 1fr; } }
  .chart-card { background: white; border-radius: 12px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); padding: 20px; }
  .chart-title { font-size: 14px; font-weight: 600; color: #374151; margin-bottom: 12px; }
  .chart-container { position: relative; }

  /* Cross-sample */
  .comparison-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
  @media (max-width: 900px) { .comparison-grid { grid-template-columns: 1fr; } }

  /* Footer */
  .footer { background: #1E293B; color: #94A3B8; text-align: center; padding: 20px; font-size: 12px; margin-top: 20px; }
  .footer a { color: #60A5FA; }

  /* Divider */
  .divider { height: 1px; background: #E2E8F0; margin: 0; }

  /* Placeholder */
  .placeholder { text-align: center; color: #94A3B8; padding: 40px; font-size: 14px; }
</style>
</head>
<body>

<!-- Header -->
<div class="header">
  <div class="container">
    <h1>{{TITLE}}</h1>
    <div class="header-meta">
      <span>📄 {{INPUT_FILE}}</span>
      <span>📅 Generated {{TIMESTAMP}}</span>
      <span class="run-mode-badge">{{RUN_MODE}}</span>
    </div>
  </div>
</div>

<!-- Summary bar -->
<div class="summary-bar">
  <div class="container">
    <div class="summary-stats" id="summary-stats"></div>
  </div>
</div>

<!-- Main content -->
<div class="container">

  <!-- Section 1: Sample Overview -->
  <div class="section">
    <div class="section-title">Sample Overview</div>
    <div class="table-wrapper">
      <table id="sample-table">
        <thead>
          <tr>
            <th data-col="sample" data-type="str">Sample<span class="sort-icon"></span></th>
            <th data-col="reads" data-type="num" class="num-cell">Reads<span class="sort-icon"></span></th>
            <th data-col="cells" data-type="num" class="num-cell">Cells<span class="sort-icon"></span></th>
            <th data-col="clones" data-type="num" class="num-cell">Clones<span class="sort-icon"></span></th>
            <th data-col="top1_pct" data-type="num" class="num-cell">Top Clone %<span class="sort-icon"></span></th>
            <th data-col="top3_pct" data-type="num" class="num-cell">Top 3 Clones %<span class="sort-icon"></span></th>
            <th data-col="top1_pct" data-type="num">Clonality<span class="sort-icon"></span></th>
          </tr>
        </thead>
        <tbody id="sample-tbody"></tbody>
      </table>
    </div>
  </div>

  <div class="divider"></div>

  <!-- Section 2: Sample Detail -->
  <div class="section">
    <div class="detail-header">
      <div class="section-title" style="margin-bottom:0">Sample Detail</div>
      <select class="detail-select" id="sample-select">
        <option value="">— Select a sample —</option>
      </select>
    </div>
    <div id="detail-placeholder" class="placeholder card">Click a row in the table above or select a sample from the dropdown to view detailed charts.</div>
    <div id="detail-charts" style="display:none;">
      <div class="charts-grid">
        <div class="chart-card">
          <div class="chart-title">A) Ranked Clone Abundance</div>
          <div class="chart-container" style="height:300px"><canvas id="chartAbundance"></canvas></div>
        </div>
        <div class="chart-card">
          <div class="chart-title">B) Clone Size Distribution</div>
          <div class="chart-container" style="height:300px"><canvas id="chartSizeDist"></canvas></div>
        </div>
        <div class="chart-card">
          <div class="chart-title">C) Top 20 Clones</div>
          <div class="chart-container" style="height:300px"><canvas id="chartTop20"></canvas></div>
        </div>
        <div class="chart-card">
          <div class="chart-title">D) Edit Distance Quality</div>
          <div class="chart-container" style="height:300px"><canvas id="chartEditDist"></canvas></div>
        </div>
      </div>
    </div>
  </div>

  <div class="divider"></div>

  <!-- Section 3: Cross-Sample Comparison -->
  <div class="section">
    <div class="section-title">Cross-Sample Comparison</div>
    <div class="comparison-grid">
      <div class="chart-card">
        <div class="chart-title">E) Cells per Sample</div>
        <div class="chart-container" style="height:320px"><canvas id="chartCellsPerSample"></canvas></div>
      </div>
      <div class="chart-card">
        <div class="chart-title">F) Clonality Comparison</div>
        <div class="chart-container" style="height:320px"><canvas id="chartClonality"></canvas></div>
      </div>
    </div>
  </div>

</div>

<!-- Footer -->
<div class="footer">
  Generated by <a href="https://github.com/phipsonlab/NextClone" target="_blank">NextClone</a> report generator &mdash; {{TIMESTAMP}}
</div>

<script>
// ============================================================
// Embedded data
// ============================================================
const DATA = {{DATA_JSON}};
const GLOBAL = {{GLOBAL_JSON}};
const SAMPLE_NAMES = Object.keys(DATA);

// ============================================================
// Utilities
// ============================================================
function fmt(n) {
  if (n === undefined || n === null) return '—';
  return Number(n).toLocaleString();
}
function pct(v) { return v.toFixed(1) + '%'; }

// ============================================================
// Summary bar
// ============================================================
function renderSummary() {
  const el = document.getElementById('summary-stats');
  const items = [
    { label: 'Total Reads', value: fmt(GLOBAL.total_reads) },
    { label: 'Total Cells', value: fmt(GLOBAL.total_cells) },
    { label: 'Samples', value: fmt(GLOBAL.total_samples) },
    { label: 'Total Clone Assignments', value: fmt(GLOBAL.total_clones) },
  ];
  el.innerHTML = items.map(i =>
    `<div class="stat-item"><div class="stat-value">${i.value}</div><div class="stat-label">${i.label}</div></div>`
  ).join('');
}

// ============================================================
// Sample table
// ============================================================
let sortCol = null, sortDir = 1;

function clonalityPill(v) {
  if (v < 10) return `<span class="pill pill-green">${pct(v)}</span>`;
  if (v < 30) return `<span class="pill pill-amber">${pct(v)}</span>`;
  return `<span class="pill pill-red">${pct(v)}</span>`;
}

function renderTable(names) {
  const tbody = document.getElementById('sample-tbody');
  tbody.innerHTML = names.map(name => {
    const s = DATA[name];
    return `<tr data-sample="${name}">
      <td>${name}</td>
      <td class="num-cell">${fmt(s.reads)}</td>
      <td class="num-cell">${fmt(s.cells)}</td>
      <td class="num-cell">${fmt(s.clones)}</td>
      <td class="num-cell">${pct(s.top1_pct)}</td>
      <td class="num-cell">${pct(s.top3_pct)}</td>
      <td>${clonalityPill(s.top1_pct)}</td>
    </tr>`;
  }).join('');

  tbody.querySelectorAll('tr').forEach(row => {
    row.addEventListener('click', () => selectSample(row.dataset.sample));
  });
}

function sortTable(col, type) {
  if (sortCol === col) sortDir *= -1;
  else { sortCol = col; sortDir = 1; }

  // update header classes
  document.querySelectorAll('th').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
    if (th.dataset.col === col) th.classList.add(sortDir === 1 ? 'sort-asc' : 'sort-desc');
  });

  const sorted = [...SAMPLE_NAMES].sort((a, b) => {
    let va = col === 'sample' ? a : DATA[a][col];
    let vb = col === 'sample' ? b : DATA[b][col];
    if (type === 'num') return (va - vb) * sortDir;
    return va.localeCompare(vb) * sortDir;
  });
  renderTable(sorted);
  // re-highlight selected
  if (currentSample) {
    document.querySelectorAll('#sample-tbody tr').forEach(r => {
      if (r.dataset.sample === currentSample) r.classList.add('selected');
    });
  }
}

document.querySelectorAll('th[data-col]').forEach(th => {
  th.addEventListener('click', () => sortTable(th.dataset.col, th.dataset.type));
});

// ============================================================
// Dropdown
// ============================================================
function populateDropdown() {
  const sel = document.getElementById('sample-select');
  SAMPLE_NAMES.forEach(name => {
    const opt = document.createElement('option');
    opt.value = name; opt.textContent = name;
    sel.appendChild(opt);
  });
  sel.addEventListener('change', e => {
    if (e.target.value) selectSample(e.target.value);
  });
}

// ============================================================
// Chart instances
// ============================================================
let charts = {};
function destroyChart(id) {
  if (charts[id]) { charts[id].destroy(); delete charts[id]; }
}

// ============================================================
// Sample selection & detail charts
// ============================================================
let currentSample = null;

function selectSample(name) {
  currentSample = name;
  // highlight row
  document.querySelectorAll('#sample-tbody tr').forEach(r => {
    r.classList.toggle('selected', r.dataset.sample === name);
  });
  // sync dropdown
  document.getElementById('sample-select').value = name;
  // show charts
  document.getElementById('detail-placeholder').style.display = 'none';
  document.getElementById('detail-charts').style.display = 'block';

  renderAbundance(name);
  renderSizeDist(name);
  renderTop20(name);
  renderEditDist(name);
}

// Chart A: Ranked Clone Abundance
function renderAbundance(name) {
  destroyChart('abundance');
  const s = DATA[name];
  const ranked = s.ranked_sizes;
  const labels = ranked.map((_, i) => i + 1);

  // Annotate top 3 with barcode labels
  const pointLabels = ranked.map((v, i) => {
    if (i < 3 && s.top_clones[i]) return s.top_clones[i].barcode;
    return null;
  });

  const ctx = document.getElementById('chartAbundance').getContext('2d');
  charts['abundance'] = new Chart(ctx, {
    type: 'line',
    data: {
      labels,
      datasets: [{
        label: 'Cells per Clone',
        data: ranked,
        borderColor: '#2563EB',
        backgroundColor: 'rgba(37,99,235,0.05)',
        borderWidth: 1.5,
        pointRadius: 0,
        fill: true,
        tension: 0.1,
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        y: {
          type: 'logarithmic',
          title: { display: true, text: 'Cells (log scale)', font: { size: 11 } },
          ticks: { callback: v => v },
        },
        x: {
          title: { display: true, text: 'Clone Rank', font: { size: 11 } },
          ticks: { maxTicksLimit: 10 },
        }
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            title: ctx => `Rank #${ctx[0].label}`,
            label: ctx => `${fmt(ctx.raw)} cells`,
            afterLabel: ctx => {
              const i = ctx.dataIndex;
              if (i < 3 && s.top_clones[i]) return `Barcode: ${s.top_clones[i].barcode}`;
              return '';
            }
          }
        },
        annotation: undefined,
      }
    },
    plugins: [{
      id: 'topAnnotations',
      afterDatasetsDraw(chart) {
        const { ctx, scales: { x, y } } = chart;
        const ds = chart.data.datasets[0];
        [0, 1, 2].forEach(i => {
          if (!s.top_clones[i] || ranked[i] === undefined) return;
          const xPx = x.getPixelForValue(i + 1);
          const yPx = y.getPixelForValue(ranked[i]);
          ctx.save();
          ctx.fillStyle = '#DC2626';
          ctx.font = '10px Inter, sans-serif';
          ctx.textAlign = 'left';
          ctx.fillText(s.top_clones[i].barcode, xPx + 4, yPx - 4);
          ctx.beginPath();
          ctx.arc(xPx, yPx, 3, 0, 2 * Math.PI);
          ctx.fillStyle = '#DC2626';
          ctx.fill();
          ctx.restore();
        });
      }
    }]
  });
}

// Chart B: Clone Size Distribution
function renderSizeDist(name) {
  destroyChart('sizedist');
  const s = DATA[name];
  const keys = ['Singleton', 'Small (2-5)', 'Medium (6-20)', 'Large (21-100)', 'Dominant (>100)'];
  const vals = keys.map(k => s.clone_size_buckets[k] || 0);
  const colors = ['#94A3B8', '#60A5FA', '#F59E0B', '#EF4444', '#DC2626'];

  const ctx = document.getElementById('chartSizeDist').getContext('2d');
  charts['sizedist'] = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: keys,
      datasets: [{ data: vals, backgroundColor: colors, borderRadius: 4 }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: { callbacks: { label: ctx => `${fmt(ctx.raw)} clones` } }
      },
      scales: {
        y: { title: { display: true, text: 'Number of Clones', font: { size: 11 } } },
        x: { ticks: { font: { size: 11 } } }
      }
    }
  });
}

// Chart C: Top 20 Clones
function renderTop20(name) {
  destroyChart('top20');
  const s = DATA[name];
  const top = s.top_clones;
  const labels = top.map(c => c.barcode).reverse();
  const values = top.map(c => c.n_cells).reverse();
  const pcts = top.map(c => c.pct).reverse();
  const colors = top.map((_, i) => {
    const ri = top.length - 1 - i; // reversed index
    if (ri < 3) return '#DC2626';
    if (ri < 10) return '#D97706';
    return '#2563EB';
  }).reverse();

  const ctx = document.getElementById('chartTop20').getContext('2d');
  charts['top20'] = new Chart(ctx, {
    type: 'bar',
    data: {
      labels,
      datasets: [{ data: values, backgroundColor: colors, borderRadius: 3 }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      indexAxis: 'y',
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: ctx => {
              const i = labels.length - 1 - ctx.dataIndex;
              return `${fmt(ctx.raw)} cells (${pcts[ctx.dataIndex]}%)`;
            }
          }
        },
        datalabels: undefined,
      },
      scales: {
        x: { title: { display: true, text: 'Number of Cells', font: { size: 11 } } },
        y: { ticks: { font: { size: 10 } } }
      }
    },
    plugins: [{
      id: 'barPctLabels',
      afterDatasetsDraw(chart) {
        const { ctx: c, scales: { x } } = chart;
        chart.data.datasets[0].data.forEach((val, i) => {
          const meta = chart.getDatasetMeta(0);
          const bar = meta.data[i];
          const pctVal = pcts[i];
          c.save();
          c.font = '10px Inter, sans-serif';
          c.fillStyle = '#374151';
          c.textAlign = 'left';
          c.textBaseline = 'middle';
          c.fillText(`${pctVal}%`, bar.x + 4, bar.y);
          c.restore();
        });
      }
    }]
  });
}

// Chart D: Edit Distance Quality
function renderEditDist(name) {
  destroyChart('editdist');
  const s = DATA[name];
  const labels = ['0', '1', '2', '3', '4', '5+'];

  const ctx = document.getElementById('chartEditDist').getContext('2d');
  charts['editdist'] = new Chart(ctx, {
    type: 'bar',
    data: {
      labels,
      datasets: [
        {
          label: 'FlankEditDist',
          data: s.flank_edit_dist,
          backgroundColor: 'rgba(37,99,235,0.7)',
          borderRadius: 3,
        },
        {
          label: 'BarcodeEditDist',
          data: s.barcode_edit_dist,
          backgroundColor: 'rgba(220,38,38,0.6)',
          borderRadius: 3,
        }
      ]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { position: 'top', labels: { font: { size: 11 } } },
        tooltip: { callbacks: { label: ctx => `${fmt(ctx.raw)} reads` } }
      },
      scales: {
        y: { title: { display: true, text: 'Number of Reads', font: { size: 11 } } },
        x: { title: { display: true, text: 'Edit Distance', font: { size: 11 } } }
      }
    }
  });
}

// ============================================================
// Cross-sample charts
// ============================================================
function renderCrossCharts() {
  // Sort by cells descending for Chart E
  const sorted = [...SAMPLE_NAMES].sort((a, b) => DATA[b].cells - DATA[a].cells);

  // Chart E: Cells per sample
  {
    const ctx = document.getElementById('chartCellsPerSample').getContext('2d');
    charts['cellsPerSample'] = new Chart(ctx, {
      type: 'bar',
      data: {
        labels: sorted,
        datasets: [{
          label: 'Unique Cells',
          data: sorted.map(n => DATA[n].cells),
          backgroundColor: '#2563EB',
          borderRadius: 4,
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        indexAxis: 'y',
        plugins: {
          legend: { display: false },
          tooltip: { callbacks: { label: ctx => `${fmt(ctx.raw)} cells` } }
        },
        scales: {
          x: { title: { display: true, text: 'Unique Cells', font: { size: 11 } } },
          y: { ticks: { font: { size: 11 } } }
        }
      }
    });
  }

  // Chart F: Clonality comparison
  {
    const ctx = document.getElementById('chartClonality').getContext('2d');
    charts['clonality'] = new Chart(ctx, {
      type: 'bar',
      data: {
        labels: SAMPLE_NAMES,
        datasets: [
          {
            label: 'Top 1%',
            data: SAMPLE_NAMES.map(n => DATA[n].top1_pct),
            backgroundColor: '#DC2626',
            borderRadius: 3,
          },
          {
            label: 'Top 3%',
            data: SAMPLE_NAMES.map(n => DATA[n].top3_pct),
            backgroundColor: '#D97706',
            borderRadius: 3,
          },
          {
            label: 'Top 10%',
            data: SAMPLE_NAMES.map(n => DATA[n].top10_pct),
            backgroundColor: '#16A34A',
            borderRadius: 3,
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: { position: 'top', labels: { font: { size: 11 } } },
          tooltip: { callbacks: { label: ctx => `${ctx.dataset.label}: ${ctx.raw.toFixed(1)}%` } }
        },
        scales: {
          y: {
            title: { display: true, text: '% of Cells', font: { size: 11 } },
            max: 100,
          },
          x: { ticks: { font: { size: 10 }, maxRotation: 30 } }
        }
      }
    });
  }
}

// ============================================================
// Init
// ============================================================
renderSummary();
renderTable(SAMPLE_NAMES);
populateDropdown();
renderCrossCharts();

// Auto-select first sample
if (SAMPLE_NAMES.length > 0) selectSample(SAMPLE_NAMES[0]);
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def detect_run_mode(stats):
    """Heuristic: if all clone barcodes look random (no common prefix/pattern), call it Discovery."""
    # We can't reliably detect reference barcodes from this CSV alone.
    # For now, default to Discovery mode unless user passes a flag.
    return "Discovery Mode"


def generate_report(csv_path, output_path, title):
    print(f"[1/4] Loading data from {csv_path}...")
    raw = load_data(csv_path)

    print(f"[2/4] Computing stats for {len(raw)} samples...")
    stats = compute_stats(raw)
    glob = global_stats(stats)

    run_mode = detect_run_mode(stats)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    input_filename = os.path.basename(csv_path)

    print(f"[3/4] Building HTML report...")
    data_json = json.dumps(stats, separators=(",", ":"))
    global_json = json.dumps(glob, separators=(",", ":"))

    html = HTML_TEMPLATE
    html = html.replace("{{TITLE}}", title)
    html = html.replace("{{INPUT_FILE}}", input_filename)
    html = html.replace("{{TIMESTAMP}}", timestamp)
    html = html.replace("{{RUN_MODE}}", run_mode)
    html = html.replace("{{DATA_JSON}}", data_json)
    html = html.replace("{{GLOBAL_JSON}}", global_json)

    print(f"[4/4] Writing to {output_path}...")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    size_kb = os.path.getsize(output_path) / 1024
    print(f"\n✅ Report generated: {output_path} ({size_kb:.1f} KB)")
    print(f"   Samples: {glob['total_samples']}")
    print(f"   Reads:   {glob['total_reads']:,}")
    print(f"   Cells:   {glob['total_cells']:,}")
    print(f"   Clones:  {glob['total_clones']:,}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a NextClone HTML report from clone_barcodes.csv"
    )
    parser.add_argument("input_csv", help="Path to clone_barcodes.csv")
    parser.add_argument("--output", default="report.html", help="Output HTML file (default: report.html)")
    parser.add_argument("--title", default="NextClone Report", help="Report title")
    args = parser.parse_args()

    if not os.path.isfile(args.input_csv):
        print(f"Error: input file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)

    generate_report(args.input_csv, args.output, args.title)


if __name__ == "__main__":
    main()
