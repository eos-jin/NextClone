#!/usr/bin/env python3
"""
NextClone Report Generator v2
Reads clone_barcodes.csv and generates a self-contained HTML dashboard.

New features (v2):
- Clone overlap table (shared clones across samples at different thresholds)
- Heterogeneity metrics (Gini coefficient, Shannon index)
- Clone size density plot (KDE curve)
- Reversed top 20 clones (largest at top)

Usage:
    python3 generate_report.py <input_csv> [--output report.html] [--title "My Run"]
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from datetime import datetime


# ---------------------------------------------------------------------------
# Data loading & stats computation
# ---------------------------------------------------------------------------

def load_data(csv_path):
    """
    Parse the CSV and return a dict of per-sample data structures.
    Also extracts run information from header lines starting with #.
    """
    samples = defaultdict(lambda: {
        "reads": 0,
        "cells": set(),
        "clone_cells": defaultdict(set),  # clone_barcode -> set of cell barcodes
        "flank_edit": defaultdict(int),
        "barcode_edit": defaultdict(int),
    })
    
    run_info = {
        "mode": None,
        "command": None,
        "parameters": {},
    }

    with open(csv_path, newline="", encoding="utf-8") as f:
        # First pass: read header comments for run information
        header_lines = []
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.strip())
            else:
                # Found first data line, reset file pointer
                f.seek(0)
                break
        
        # Parse header comments
        for line in header_lines:
            line = line.lstrip('#').strip()
            if line.startswith('mode:'):
                run_info["mode"] = line.split(':', 1)[1].strip()
            elif line.startswith('command:'):
                run_info["command"] = line.split(':', 1)[1].strip()
            elif ':' in line:
                key, val = line.split(':', 1)
                run_info["parameters"][key.strip()] = val.strip()
        
        # Second pass: read CSV data
        f.seek(0)
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

    return samples, run_info


def compute_gini(values):
    """
    Calculate Gini coefficient for clone size distribution.
    0 = perfect equality (all clones same size)
    1 = perfect inequality (one clone has all cells)
    
    Formula: G = sum(|xi - xj|) / (2 * n * sum(x))
    """
    if not values or sum(values) == 0:
        return 0.0
    
    n = len(values)
    if n == 1:
        return 0.0
    
    # Sort values
    sorted_vals = sorted(values)
    
    # Calculate Gini using the efficient formula
    # G = (2 * sum(i * x_i) - (n + 1) * sum(x_i)) / (n * sum(x_i))
    total = sum(sorted_vals)
    weighted_sum = sum((i + 1) * val for i, val in enumerate(sorted_vals))
    
    gini = (2 * weighted_sum - (n + 1) * total) / (n * total)
    return round(gini, 4)


def compute_shannon(values):
    """
    Calculate Shannon diversity index for clone distribution.
    Higher = more diverse (many clones with similar sizes)
    Lower = less diverse (few dominant clones)
    
    Formula: H = -sum(pi * ln(pi))
    """
    if not values or sum(values) == 0:
        return 0.0
    
    total = sum(values)
    h = 0.0
    
    for val in values:
        if val > 0:
            pi = val / total
            h -= pi * math.log(pi)
    
    return round(h, 4)


def compute_clone_overlap(samples):
    """
    Compute clone overlap across samples at different cell thresholds.
    Returns a dict with thresholds as keys and per-sample counts + "in_all" count.
    """
    thresholds = [5, 10, 15, 20, 50, 100]
    sample_names = sorted(samples.keys())
    
    # For each sample, get clones with >= threshold cells
    sample_clone_sets = {}
    for sample, raw in samples.items():
        clone_sizes = {clone: len(cells) for clone, cells in raw["clone_cells"].items()}
        sample_clone_sets[sample] = clone_sizes
    
    # Compute overlap for each threshold
    overlap_data = {}
    for thresh in thresholds:
        overlap_data[thresh] = {
            "per_sample": {},
            "in_all": 0
        }
        
        # Get clones meeting threshold for each sample
        clones_above_thresh = {}
        for sample in sample_names:
            clones_above = [
                clone for clone, size in sample_clone_sets[sample].items()
                if size >= thresh
            ]
            clones_above_thresh[sample] = set(clones_above)
            overlap_data[thresh]["per_sample"][sample] = len(clones_above)
        
        # Clones present in ALL samples above threshold
        if len(sample_names) > 0:
            common_clones = set.intersection(*clones_above_thresh.values())
            overlap_data[thresh]["in_all"] = len(common_clones)
    
    return overlap_data


def compute_stats(samples):
    """Turn raw per-sample data into serialisable stats dicts."""
    result = {}
    for sample, raw in sorted(samples.items()):
        n_reads = raw["reads"]
        n_cells = len(raw["cells"])

        # Clone sizes (by unique cells per clone)
        clone_sizes = {clone: len(cells) for clone, cells in raw["clone_cells"].items()}
        n_clones = len(clone_sizes)
        clone_size_values = list(clone_sizes.values())

        # Ranked sizes (descending)
        ranked = sorted(clone_size_values, reverse=True)

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

        # Clone size density (for KDE plot)
        # Create binned density for log-transformed clone sizes
        density_data = []
        if clone_size_values:
            # Use log scale for better visualization
            log_sizes = [math.log10(sz) for sz in clone_size_values if sz > 0]
            if log_sizes:
                min_log = min(log_sizes)
                max_log = max(log_sizes)
                n_bins = min(30, len(log_sizes))
                if n_bins > 1:
                    bin_width = (max_log - min_log) / n_bins
                    for i in range(n_bins):
                        bin_start = min_log + i * bin_width
                        bin_count = sum(1 for ls in log_sizes if bin_start <= ls < bin_start + bin_width)
                        density_data.append({
                            "x": round(10 ** bin_start, 2),
                            "y": bin_count
                        })

        # Top 20 clones (already sorted descending - largest first)
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

        # Heterogeneity metrics
        gini = compute_gini(clone_size_values)
        shannon = compute_shannon(clone_size_values)

        result[sample] = {
            "reads": n_reads,
            "cells": n_cells,
            "clones": n_clones,
            "ranked_sizes": ranked,
            "clone_size_buckets": buckets,
            "clone_size_density": density_data,
            "top_clones": top_clones,
            "top1_pct": top_n_pct(1),
            "top3_pct": top_n_pct(3),
            "top10_pct": top_n_pct(10),
            "flank_edit_dist": ed_dist(raw["flank_edit"]),
            "barcode_edit_dist": ed_dist(raw["barcode_edit"]),
            "gini": gini,
            "shannon": shannon,
        }

    return result


def compute_global_overlap(samples):
    """Compute clone overlap data for all samples."""
    return compute_clone_overlap(samples)


def global_stats(stats):
    total_reads = sum(s["reads"] for s in stats.values())
    total_cells = sum(s["cells"] for s in stats.values())
    total_samples = len(stats)
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

  /* Heterogeneity badge */
  .het-badge { display: inline-block; padding: 2px 8px; border-radius: 6px; font-size: 11px; font-weight: 600; margin-left: 6px; }
  .het-low { background: #DCFCE7; color: #16A34A; }
  .het-med { background: #FEF3C7; color: #D97706; }
  .het-high { background: #FEE2E2; color: #DC2626; }

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
  
  /* Overlap table */
  .overlap-table-wrapper { overflow-x: auto; border-radius: 12px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); margin-top: 16px; }
  .overlap-table { width: 100%; border-collapse: collapse; }
  .overlap-table th { background: #F1F5F9; font-size: 11px; padding: 10px 12px; text-align: right; border: 1px solid #E2E8F0; }
  .overlap-table th:first-child { text-align: left; background: #F8FAFC; }
  .overlap-table td { padding: 10px 12px; text-align: right; border: 1px solid #E2E8F0; font-variant-numeric: tabular-nums; }
  .overlap-table tr:hover { background: #F8FAFC; }
  .overlap-table th:first-child, .overlap-table td:first-child { text-align: left; background: #F8FAFC; font-weight: 600; color: #1E293B; }
  .overlap-table .in-all-col { background: #DBEAFE; font-weight: 600; color: #1E40AF; }

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
            <th data-col="gini" data-type="num" class="num-cell">Gini<span class="sort-icon"></span></th>
            <th data-col="shannon" data-type="num" class="num-cell">Shannon<span class="sort-icon"></span></th>
          </tr>
        </thead>
        <tbody id="sample-tbody"></tbody>
      </table>
    </div>
  </div>

  <div class="divider"></div>

  <!-- Section 2: Clone Overlap Across Samples -->
  <div class="section">
    <div class="section-title">Clone Overlap Across Samples</div>
    <div class="card">
      <p style="color:#64748B;font-size:13px;margin-bottom:12px;">
        Number of clones detected in each sample (and in ALL samples) at different cell count thresholds.
        Higher overlap indicates consistent clone detection across samples.
      </p>
      <div class="overlap-table-wrapper">
        <table class="overlap-table" id="overlap-table">
          <thead>
            <tr id="overlap-header"></tr>
          </thead>
          <tbody id="overlap-tbody"></tbody>
        </table>
      </div>
    </div>
  </div>

  <div class="divider"></div>

  <!-- Section 3: Sample Detail -->
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
          <div class="chart-title">B) Clone Size Density</div>
          <div class="chart-container" style="height:300px"><canvas id="chartSizeDensity"></canvas></div>
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

  <!-- Section 5: Cross-Sample Comparison -->
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
const OVERLAP = {{OVERLAP_JSON}};
const SAMPLE_NAMES = Object.keys(DATA);

// ============================================================
// Utilities
// ============================================================
function fmt(n) {
  if (n === undefined || n === null) return '—';
  return Number(n).toLocaleString();
}
function pct(v) { return v.toFixed(1) + '%'; }
function fmt4(v) { return v.toFixed(4); }

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

function giniBadge(v) {
  if (v < 0.3) return `<span class="het-badge het-low">Low</span>`;
  if (v < 0.6) return `<span class="het-badge het-med">Med</span>`;
  return `<span class="het-badge het-high">High</span>`;
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
      <td class="num-cell">${fmt4(s.gini)} ${giniBadge(s.gini)}</td>
      <td class="num-cell">${fmt4(s.shannon)}</td>
    </tr>`;
  }).join('');

  tbody.querySelectorAll('tr').forEach(row => {
    row.addEventListener('click', () => selectSample(row.dataset.sample));
  });
}

function sortTable(col, type) {
  if (sortCol === col) sortDir *= -1;
  else { sortCol = col; sortDir = 1; }

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
// Overlap table
// ============================================================
function renderOverlapTable() {
  const header = document.getElementById('overlap-header');
  const tbody = document.getElementById('overlap-tbody');
  
  // Header row
  const thresholds = Object.keys(OVERLAP).map(Number);
  header.innerHTML = '<th>Threshold</th>' + 
    SAMPLE_NAMES.map(s => `<th>${s}</th>`).join('') +
    '<th class="in-all-col">In ALL Samples</th>';
  
  // Data rows
  tbody.innerHTML = thresholds.map(thresh => {
    const row = OVERLAP[thresh];
    return '<tr>' +
      `<td>≥${thresh} cells</td>` +
      SAMPLE_NAMES.map(s => `<td>${fmt(row.per_sample[s])}</td>`).join('') +
      `<td class="in-all-col">${fmt(row.in_all)}</td>` +
      '</tr>';
  }).join('');
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
  document.querySelectorAll('#sample-tbody tr').forEach(r => {
    r.classList.toggle('selected', r.dataset.sample === name);
  });
  document.getElementById('sample-select').value = name;
  document.getElementById('detail-placeholder').style.display = 'none';
  document.getElementById('detail-charts').style.display = 'block';

  renderAbundance(name);
  renderSizeDensity(name);
  renderTop20(name);
  renderEditDist(name);
}

// Chart A: Ranked Clone Abundance
function renderAbundance(name) {
  destroyChart('abundance');
  const s = DATA[name];
  const ranked = s.ranked_sizes;
  const labels = ranked.map((_, i) => i + 1);

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

// Chart B: Clone Size Density (KDE-style)
function renderSizeDensity(name) {
  destroyChart('sizedensity');
  const s = DATA[name];
  const densityData = s.clone_size_density;
  
  if (!densityData || densityData.length === 0) {
    const ctx = document.getElementById('chartSizeDensity').getContext('2d');
    ctx.canvas.parentNode.innerHTML = '<div class="placeholder">No density data available</div>';
    return;
  }
  
  const labels = densityData.map(d => fmt(d.x));
  const values = densityData.map(d => d.y);

  const ctx = document.getElementById('chartSizeDensity').getContext('2d');
  charts['sizedensity'] = new Chart(ctx, {
    type: 'line',
    data: {
      labels,
      datasets: [{
        label: 'Clone Count',
        data: values,
        borderColor: '#16A34A',
        backgroundColor: 'rgba(22,163,74,0.1)',
        borderWidth: 2,
        pointRadius: 0,
        fill: true,
        tension: 0.3,
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        y: {
          title: { display: true, text: 'Number of Clones', font: { size: 11 } },
          beginAtZero: true,
        },
        x: {
          title: { display: true, text: 'Clone Size (cells, log scale)', font: { size: 11 } },
          min: 0,
          ticks: { maxTicksLimit: 10 }
        }
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            title: ctx => `Clone size: ~${ctx[0].label} cells`,
            label: ctx => `${fmt(ctx.raw)} clones`,
          }
        }
      }
    }
  });
}

// Chart C: Top 20 Clones (reversed - largest at top)
function renderTop20(name) {
  destroyChart('top20');
  const s = DATA[name];
  const top = s.top_clones;
  
  // Reverse so largest is at top (index 0)
  const labels = top.map(c => c.barcode);
  const values = top.map(c => c.n_cells);
  const pcts = top.map(c => c.pct);
  const colors = top.map((_, i) => {
    if (i < 3) return '#DC2626';
    if (i < 10) return '#D97706';
    return '#2563EB';
  });

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
            label: ctx => `${fmt(ctx.raw)} cells (${pcts[ctx.dataIndex]}%)`
          }
        },
      },
      scales: {
        x: { title: { display: true, text: 'Number of Cells', font: { size: 11 } } },
        y: { 
          ticks: { font: { size: 10 } },
          reverse: false  // Largest at top
        }
      }
    },
    plugins: [{
      id: 'barPctLabels',
      afterDatasetsDraw(chart) {
        const { ctx: c } = chart;
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
renderOverlapTable();
renderCrossCharts();

if (SAMPLE_NAMES.length > 0) selectSample(SAMPLE_NAMES[0]);
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def detect_run_mode(run_info):
    """Determine run mode from run_info or return 'Unknown'."""
    if run_info.get("mode"):
        mode = run_info["mode"]
        # Capitalize appropriately
        if mode.lower() == "discovery":
            return "Discovery Mode"
        elif mode.lower() == "whitelist" or mode.lower() == "reference":
            return "Whitelist Mode"
        else:
            return mode
    return "Run Mode Unknown"


def generate_report(csv_path, output_path, title):
    print(f"[1/6] Loading data from {csv_path}...")
    raw, run_info = load_data(csv_path)

    print(f"[2/6] Computing stats for {len(raw)} samples...")
    stats = compute_stats(raw)
    glob = global_stats(stats)
    overlap = compute_global_overlap(raw)

    run_mode = detect_run_mode(run_info)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    input_filename = os.path.basename(csv_path)

    print(f"[3/6] Building HTML report...")
    data_json = json.dumps(stats, separators=(",", ":"))
    global_json = json.dumps(glob, separators=(",", ":"))
    overlap_json = json.dumps(overlap, separators=(",", ":"))

    html = HTML_TEMPLATE
    html = html.replace("{{TITLE}}", title)
    html = html.replace("{{INPUT_FILE}}", input_filename)
    html = html.replace("{{TIMESTAMP}}", timestamp)
    html = html.replace("{{RUN_MODE}}", run_mode)
    html = html.replace("{{DATA_JSON}}", data_json)
    html = html.replace("{{GLOBAL_JSON}}", global_json)
    html = html.replace("{{OVERLAP_JSON}}", overlap_json)

    print(f"[4/6] Writing to {output_path}...")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"[5/6] Complete!")
    print(f"[6/6] Summary:")
    size_kb = os.path.getsize(output_path) / 1024
    print(f"\n✅ Report generated: {output_path} ({size_kb:.1f} KB)")
    print(f"   Samples: {glob['total_samples']}")
    print(f"   Reads:   {glob['total_reads']:,}")
    print(f"   Cells:   {glob['total_cells']:,}")
    print(f"   Clones:  {glob['total_clones']:,}")
    if run_info.get("mode"):
        print(f"   Mode:    {run_mode}")
    if run_info.get("command"):
        print(f"   Command: {run_info['command'][:80]}...")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a NextClone HTML report from clone_barcodes.csv (v2)"
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
