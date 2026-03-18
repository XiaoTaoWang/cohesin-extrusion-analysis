#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enhanced K562 Micro-C subtraction heatmaps for simulation.
  - CTCF_degron - WT
  - RAD21_degron - WT
Region: chr4:37302000-38102000
Uses diagonal normalization to align sequencing depth before subtraction.
"""

# ── Path configuration ─────────────────────────────────────────────
import os
_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
# scripts/pairwise/region_visualization -> scripts/pairwise -> scripts -> K562_chromatin_folding
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import math
import cooler
import pyBigWig
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap, Normalize

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

# =============================================================================
#  Helper Functions
# =============================================================================
def plot_y_axis(ax, max_val, label_size=10):
    max_str = str(int(max_val)) if max_val % 1 == 0 else f"{max_val:.1f}"
    ax_pos = ax.get_position().bounds
    gap, width = 0.015, 0.005
    y_ax = ax.figure.add_axes([ax_pos[0] - gap, ax_pos[1], width, ax_pos[3]])
    y_ax.plot([1, 1], [0, 1], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.plot([0, 1], [0, 0], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.plot([0, 1], [1, 1], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.text(-0.5, 0, "0", ha='right', va='center', fontsize=label_size, transform=y_ax.transAxes)
    y_ax.text(-0.5, 1, max_str, ha='right', va='center', fontsize=label_size, transform=y_ax.transAxes)
    y_ax.axis('off')

def normalize_diag(m_wt, m_ko):
    m_wt_clean = np.nan_to_num(m_wt)
    m_ko_clean = np.nan_to_num(m_ko)
    diag_wt = np.sum(np.diagonal(m_wt_clean))
    diag_ko = np.sum(np.diagonal(m_ko_clean))
    coeff = (diag_wt / diag_ko) if diag_ko > 0 else 1.0
    return m_ko_clean * coeff, coeff

def pad_matrix(m, target_n):
    if m.shape[0] < target_n:
        return np.pad(m, ((0, target_n - m.shape[0]), (0, target_n - m.shape[0])), 'constant')
    return m

def find_chrom(clr, chrom):
    if chrom in clr.chromnames: return chrom
    if 'chr' + chrom in clr.chromnames: return 'chr' + chrom
    if chrom.startswith('chr') and chrom[3:] in clr.chromnames: return chrom[3:]
    return chrom

def get_raw_max_value(bw_path, chrom, start, end):
    try:
        bw = pyBigWig.open(bw_path)
        chrom_in_bw = next((v for v in [chrom, chrom.lstrip('chr'), f"chr{chrom.lstrip('chr')}"] if v in bw.chroms()), None)
        if chrom_in_bw is None: return 0
        max_val = bw.stats(chrom_in_bw, start, end, type="max")[0]
        bw.close()
        return max_val if max_val is not None else 0
    except Exception:
        return 0

# =============================================================================
#  Plotting Class
# =============================================================================
class SubtractionPlot:

    def __init__(self, chrom, start, end, n_bins, res, figsize=(12, 14),
                 n_rows=13, track_partition=None, space=0.08):
        self.fig = plt.figure(figsize=figsize, facecolor='white')
        self.grid = GridSpec(n_rows, 1, figure=self.fig, left=0.15, right=0.92,
                             bottom=0.03, top=0.97, hspace=space, height_ratios=track_partition)
        self.track_count = 0
        self.chrom = chrom
        self.start = start
        self.end = end
        self.n_bins = n_bins
        self.res = res
        self.cmap_diff = plt.get_cmap('RdYlBu_r')
        self.cmap_diff.set_bad(color='white')

    def plot_difference_matrix(self, diff, vmax=None, title=None):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        if vmax is None:
            flat = np.abs(diff).flatten()
            flat_valid = flat[flat > 0]
            vmax = np.percentile(flat_valid, 95) if len(flat_valid) > 0 else 1

        im = h_ax.imshow(diff, cmap=self.cmap_diff, vmax=vmax, vmin=-vmax,
                         interpolation='none', origin='upper', aspect='equal')
        h_ax.axis('off')

        hb = h_ax.get_position().bounds
        w, h = hb[2], hb[3]
        size = min(w, h)
        h_ax.set_position([hb[0] + (w - size) / 2, hb[1] + (h - size) / 2, size, size])
        self.heatmap_pos = h_ax.get_position().bounds

        if title:
            h_ax.set_title(title, fontsize=14, pad=10)

        cbar_ax = self.fig.add_axes([self.heatmap_pos[0] + size + 0.01,
                                     self.heatmap_pos[1], 0.02, size * 0.2])
        self.fig.colorbar(im, cax=cbar_ax)
        cbar_ax.tick_params(labelsize=9)

    def plot_coordinates(self, labelsize=10):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        if hasattr(self, 'heatmap_pos'):
            pos = ax.get_position().bounds
            ax.set_position([self.heatmap_pos[0], pos[1], self.heatmap_pos[2], pos[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_xticks([0, self.n_bins])
        ax.set_xticklabels([f"{self.start:,}", f"{self.end:,}"], fontsize=labelsize)
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False, length=3)
        ax.set_yticks([])
        ax.text(-0.02, 0.5, self.chrom, fontsize=labelsize, va='center', ha='right', transform=ax.transAxes)

    def plot_motif(self, input_fil, track_name, strand='+', label_size=10, marker_size=60):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        loci = []
        chrom_variants = {self.chrom, self.chrom.lstrip('chr'), f"chr{self.chrom.lstrip('chr')}"}
        try:
            with open(input_fil, 'r') as source:
                for line in source:
                    parts = line.rstrip().split()
                    if len(parts) < 3: continue
                    c, s, e = parts[0], int(parts[1]), int(parts[2])
                    file_strand = parts[3] if len(parts) >= 4 else '+'
                    if (c in chrom_variants) and (s >= self.start) and (e <= self.end) and (file_strand == strand):
                        relative_pos = (((s + e) / 2) - self.start) / self.res
                        loci.append(relative_pos)
        except Exception:
            pass
        if loci:
            if strand == '+':
                ax.scatter(loci, [0] * len(loci), s=marker_size, c='red', marker='>', edgecolors='none', clip_on=True)
            else:
                ax.scatter(loci, [0] * len(loci), s=marker_size, c='blue', marker='<', edgecolors='none', clip_on=True)
        if hasattr(self, 'heatmap_pos'):
            pos = ax.get_position().bounds
            ax.set_position([self.heatmap_pos[0], pos[1], self.heatmap_pos[2], pos[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_ylim(-1, 1)
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(track_name, rotation=0, fontsize=label_size, ha='right', va='center', labelpad=15)

    def plot_signal(self, track_name, bw_fil, color='#666666', max_value=None, label_size=10, nBins=1000):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        arr = np.array([])
        try:
            with pyBigWig.open(bw_fil) as db:
                chrom_to_fetch = next((v for v in [self.chrom, self.chrom.lstrip('chr'),
                                                   f"chr{self.chrom.lstrip('chr')}"] if v in db.chroms()), None)
                if chrom_to_fetch:
                    stats = db.stats(chrom_to_fetch, self.start, self.end, nBins=nBins, type="max")
                    arr = np.array(stats).astype(float)
        except Exception:
            pass
        arr[np.isnan(arr)] = 0
        if max_value is None:
            max_value = np.max(arr) if arr.size > 0 else 1
        x = np.linspace(0, self.n_bins, len(arr))
        ax.fill_between(x, arr, color=color, edgecolor='none', linewidth=0)
        if hasattr(self, 'heatmap_pos'):
            pos = ax.get_position().bounds
            ax.set_position([self.heatmap_pos[0], pos[1], self.heatmap_pos[2], pos[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_ylim(0, max_value)
        ax.axis('off')
        ax.text(-0.02, 0.5, track_name, transform=ax.transAxes, fontsize=label_size, ha='right', va='center')
        plot_y_axis(ax, max_value, label_size=label_size - 2)

    def outfig(self, outfile, dpi=300):
        self.fig.savefig(outfile, dpi=dpi, bbox_inches='tight')
        print(f"  Saved: {outfile}")


# =============================================================================
#  Main
# =============================================================================
if __name__ == '__main__':

    # ── Region ──
    CHROM = 'chr4'
    START = 37302000
    END   = 38102000
    REGION_NAME = f"{CHROM}_{START}_{END}"

    # ── Output directory (parent of region_visualization, ending with for_simulation) ──
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    OUT_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "region_visualization_for_simulation")
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Enhanced cool files (K562 Micro-C, 5kb resolution) ──
    ENHANCE_BASE = f"{DATA_ROOT}/data_enhancement_1kb_cool_LY"
    cool_wt    = os.path.join(ENHANCE_BASE, "K562_MC_WT_1kb.cool_5000", "HiCFoundation_enhanced.cool")
    cool_ctcf  = os.path.join(ENHANCE_BASE, "K562_MC_CTCF_6h_9h_merge_007_026_1kb.cool_5000", "HiCFoundation_enhanced.cool")
    cool_rad21 = os.path.join(ENHANCE_BASE, "K562_MC_RAD21_6h_9h_merge_024_031_1kb.cool_5000", "HiCFoundation_enhanced.cool")

    # ── Motif & signal tracks (matching K562_example_2kb.py) ──
    MOTIF_FILE = f"{PROJECT_ROOT}/example_for_simulation/K562.CTCF_merged.bed"

    BW_BASE = f"{DATA_ROOT}/mly/K562_CUT_TAG/4.IgG_signal_blocking"
    bw_files = {
        'CTCF':    f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_CTCF_h019_153g_deep.RPGC.bw",
        'Cohesin': f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_RAD21_h020_160g_deep.RPGC.bw",
        'NIPBL':   os.path.join(BW_BASE, "CTCF_0h_NIPBL_positive.bw"),
        'WAPL':    os.path.join(BW_BASE, "CTCF_0h_WAPL_positive.bw"),
        'RNAPII':  os.path.join(BW_BASE, "CTCF_0h_Pol2_positive.bw"),
        'H3K27ac': f"{DATA_ROOT}/mly/K562_CUT_TAG/3.loading_site_heatmap/K562_H3K27ac_ChIP-seq_ENCFF465GBD.bigWig",
        'H3K4me1': f"{DATA_ROOT}/mly/K562_CUT_TAG/3.loading_site_heatmap/K562_H3K4me1_ChIP-seq_ENCFF783QIW.bigWig",
    }
    track_cfg = [
        ('CTCF',    '#0000B2'),
        ('Cohesin', '#008000'),
        ('NIPBL',   '#E69F00'),
        ('WAPL',    '#8B0000'),
        ('RNAPII',  '#800080'),
        ('H3K27ac', '#CD3700'),
        ('H3K4me1', '#CD0000'),
    ]

    # ── Load matrices ──
    print(f"Loading enhanced cool files for {REGION_NAME} ...")
    clr_wt    = cooler.Cooler(cool_wt)
    clr_ctcf  = cooler.Cooler(cool_ctcf)
    clr_rad21 = cooler.Cooler(cool_rad21)
    res = clr_wt.binsize

    c_wt    = find_chrom(clr_wt, CHROM)
    c_ctcf  = find_chrom(clr_ctcf, CHROM)
    c_rad21 = find_chrom(clr_rad21, CHROM)

    m_wt    = clr_wt.matrix(balance=False, sparse=False).fetch((c_wt, START, END))
    m_ctcf  = clr_ctcf.matrix(balance=False, sparse=False).fetch((c_ctcf, START, END))
    m_rad21 = clr_rad21.matrix(balance=False, sparse=False).fetch((c_rad21, START, END))

    m_wt[np.isnan(m_wt)] = 0
    m_ctcf[np.isnan(m_ctcf)] = 0
    m_rad21[np.isnan(m_rad21)] = 0

    # Pad to same size
    n_max = max(m_wt.shape[0], m_ctcf.shape[0], m_rad21.shape[0])
    m_wt    = pad_matrix(m_wt, n_max)
    m_ctcf  = pad_matrix(m_ctcf, n_max)
    m_rad21 = pad_matrix(m_rad21, n_max)

    # ── Diagonal normalization ──
    m_ctcf_norm, coef_ctcf   = normalize_diag(m_wt, m_ctcf)
    m_rad21_norm, coef_rad21 = normalize_diag(m_wt, m_rad21)
    print(f"  Diagonal normalization coefficients: CTCF={coef_ctcf:.4f}, RAD21={coef_rad21:.4f}")

    # ── Compute difference matrices ──
    diff_ctcf  = m_ctcf_norm - m_wt
    diff_rad21 = m_rad21_norm - m_wt

    # Shared vmax across both diffs
    all_diff_vals = np.concatenate([np.abs(diff_ctcf)[np.abs(diff_ctcf) > 0],
                                    np.abs(diff_rad21)[np.abs(diff_rad21) > 0]])
    shared_diff_vmax = np.percentile(all_diff_vals, 95) if len(all_diff_vals) > 0 else 1

    # ── Pre-calculate signal track max values ──
    max_vals = {}
    for protein, path in bw_files.items():
        raw = get_raw_max_value(path, CHROM, START, END)
        max_vals[protein] = math.ceil(raw / 10.0) * 10 if raw > 0 else 10

    # ── Plot setup ──
    n_tracks = len(track_cfg)
    track_ratios = [30, 1, 1.5, 1.5] + [2] * n_tracks
    n_rows = len(track_ratios)

    # ===== CTCF degron - WT =====
    print("Plotting: CTCF_degron - WT ...")
    vis_ctcf = SubtractionPlot(CHROM, START, END, n_max, res,
                                figsize=(12, 14), n_rows=n_rows, track_partition=track_ratios)
    vis_ctcf.plot_difference_matrix(diff_ctcf, vmax=shared_diff_vmax,
                                     title=f"Enhanced MC: CTCF degron - WT [diag x{coef_ctcf:.2f}]")
    vis_ctcf.plot_coordinates()
    vis_ctcf.plot_motif(MOTIF_FILE, 'Motif (+)', strand='+')
    vis_ctcf.plot_motif(MOTIF_FILE, 'Motif (-)', strand='-')
    for track_label, color in track_cfg:
        vis_ctcf.plot_signal(track_label, bw_files[track_label], color=color,
                              max_value=max_vals.get(track_label, 10))
    vis_ctcf.outfig(os.path.join(OUT_DIR, f"enhanced_MC_diff_CTCF-degron_vs_WT_{res}bp_{REGION_NAME}.pdf"))
    plt.close(vis_ctcf.fig)

    # ===== RAD21 degron - WT =====
    print("Plotting: RAD21_degron - WT ...")
    vis_rad21 = SubtractionPlot(CHROM, START, END, n_max, res,
                                 figsize=(12, 14), n_rows=n_rows, track_partition=track_ratios)
    vis_rad21.plot_difference_matrix(diff_rad21, vmax=shared_diff_vmax,
                                      title=f"Enhanced MC: RAD21 degron - WT [diag x{coef_rad21:.2f}]")
    vis_rad21.plot_coordinates()
    vis_rad21.plot_motif(MOTIF_FILE, 'Motif (+)', strand='+')
    vis_rad21.plot_motif(MOTIF_FILE, 'Motif (-)', strand='-')
    for track_label, color in track_cfg:
        vis_rad21.plot_signal(track_label, bw_files[track_label], color=color,
                               max_value=max_vals.get(track_label, 10))
    vis_rad21.outfig(os.path.join(OUT_DIR, f"enhanced_MC_diff_RAD21-degron_vs_WT_{res}bp_{REGION_NAME}.pdf"))
    plt.close(vis_rad21.fig)

    print(f"\nAll done! Output directory: {OUT_DIR}")
