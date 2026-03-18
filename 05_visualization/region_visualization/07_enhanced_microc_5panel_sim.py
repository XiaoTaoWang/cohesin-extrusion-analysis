#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import re
import argparse
import cooler
import pyBigWig
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap, Normalize

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR))))
DATA_ROOT = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

def parse_region(region_str):
    region_str = region_str.replace(',', '')
    match = re.match(r'(chr)?(\d+|[XYM]):(\d+)-(\d+)', region_str)
    if not match:
        raise ValueError(
            f"Invalid region format: {region_str}. "
            "Use format: 8:104200000-104900000 or chr8:104200000-104900000"
        )
    chrom = 'chr' + match.group(2) if not match.group(2).startswith('chr') else match.group(2)
    return chrom, int(match.group(3)), int(match.group(4))

def parse_resolution_arg(resolution_arg):
    if resolution_arg is None:
        return 2000
    s = str(resolution_arg).strip().lower().replace(',', '')
    if s.endswith('kb'):
        value = int(float(s[:-2]) * 1000)
    elif s.endswith('k'):
        value = int(float(s[:-1]) * 1000)
    elif s.endswith('mb'):
        value = int(float(s[:-2]) * 1000000)
    elif s.endswith('m'):
        value = int(float(s[:-1]) * 1000000)
    elif s.endswith('bp'):
        value = int(float(s[:-2]))
    else:
        value = int(float(s))

    if value <= 0:
        raise ValueError(f"Resolution must be positive, got: {resolution_arg}")
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description="Enhanced Micro-C 5-panel plot with optional resolution parameter"
    )
    parser.add_argument(
        "region",
        help="Genomic region, e.g. 1:236093000-236501000 or chr1:236,093,000-236,501,000"
    )
    parser.add_argument(
        "resolution",
        nargs='?',
        default="2000",
        help="Optional resolution, e.g. 2000 / 2kb / 5kb (default: 2000)"
    )
    parser.add_argument(
        "-r", "--res",
        dest="res_opt",
        default=None,
        help="Optional resolution override, same formats as positional resolution"
    )
    args = parser.parse_args()

    region_str = args.region
    resolution_str = args.res_opt if args.res_opt is not None else args.resolution
    resolution_bp = parse_resolution_arg(resolution_str)
    return region_str, resolution_bp

def plot_y_axis(ax, max_val, label_size=8):
    max_str = str(int(max_val)) if max_val % 1 == 0 else f"{max_val:.1f}"
    ax_pos = ax.get_position().bounds
    y_ax = ax.figure.add_axes([ax_pos[0] - 0.02, ax_pos[1], 0.005, ax_pos[3]])
    for spine in ['left', 'top', 'bottom', 'right']: y_ax.spines[spine].set_visible(False)
    y_ax.plot([1, 1], [0, 1], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.plot([0, 1], [0, 0], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.plot([0, 1], [1, 1], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.text(-0.8, 0, "0", ha='right', va='center', fontsize=label_size, transform=y_ax.transAxes)
    y_ax.text(-0.8, 1, max_str, ha='right', va='center', fontsize=label_size, transform=y_ax.transAxes)
    y_ax.axis('off')

def pad_matrix(m, target_n):
    return np.pad(m, ((0, target_n - m.shape[0]), (0, target_n - m.shape[0])), 'constant') if m.shape[0] < target_n else m

def find_chrom(clr, chrom):
    if chrom in clr.chromnames: return chrom
    if 'chr' + chrom in clr.chromnames: return 'chr' + chrom
    if chrom.startswith('chr') and chrom[3:] in clr.chromnames: return chrom[3:]
    return chrom

def calc_shared_raw_vmax(matrices, percentile=93):
    all_vals = [m[m > 0] for m in matrices if len(m[m > 0]) > 0]
    return np.percentile(np.concatenate(all_vals), percentile) if all_vals else 1

def calc_shared_diff_vmax(diff_matrices, percentile=95):
    abs_vals = [np.abs(d)[np.abs(d) > 0] for d in diff_matrices if len(np.abs(d)[np.abs(d) > 0]) > 0]
    return np.percentile(np.concatenate(abs_vals), percentile) if abs_vals else 1

class FivePanelPlot:
    def __init__(self, chrom, start, end, n_bins, res, figsize=(12, 32), n_rows=14, track_partition=None, space=0.01):
        self.fig = plt.figure(figsize=figsize, facecolor='white')
        # simple comment
        self.grid = GridSpec(n_rows, 1, figure=self.fig, left=0.1, right=0.9, bottom=0.06, top=0.97, hspace=space, height_ratios=track_partition)
        self.track_count, self.chrom, self.start, self.end, self.n_bins, self.res = 0, chrom, start, end, n_bins, res
        self.raw_cmap = LinearSegmentedColormap.from_list('red_white', ['#FFFFFF', '#FFDFDF', '#FF7575', '#FF2626', '#F70000'], N=256)
        self.raw_cmap.set_bad(color='white')
        self.diff_cmap = plt.get_cmap('RdYlBu_r').copy()
        self.diff_cmap.set_bad(color='white')

    def _finalize_heatmap_axis(self, h_ax, title=None):
        h_ax.axis('off')
        hb, (fig_w, fig_h) = h_ax.get_position().bounds, self.fig.get_size_inches()
        size_inches = min(hb[2] * fig_w, hb[3] * fig_h)
        rel_w, rel_h = size_inches / fig_w, size_inches / fig_h
        new_x, new_y = hb[0] + (hb[2] - rel_w) / 2, hb[1] + (hb[3] - rel_h) / 2
        h_ax.set_position([new_x, new_y, rel_w, rel_h])
        self.heatmap_pos = [new_x, new_y, rel_w, rel_h]
        if title: h_ax.set_title(title, fontsize=13, pad=10)
        return rel_w

    def _add_colorbar(self, im, rel_w):
        cbar_ax = self.fig.add_axes([self.heatmap_pos[0] + rel_w + 0.015, self.heatmap_pos[1] + self.heatmap_pos[3] * 0.2, 0.015, self.heatmap_pos[3] * 0.4])
        self.fig.colorbar(im, cax=cbar_ax).ax.tick_params(labelsize=8)

    def plot_split_matrix(self, m_wt, m_aid, vmax, title=None):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        m_plot = np.triu(m_wt.astype(float), k=1) + np.tril(m_aid.astype(float), k=0)
        m_plot[m_plot <= 0] = np.nan
        im = h_ax.imshow(m_plot, cmap=self.raw_cmap, norm=Normalize(vmin=0, vmax=vmax), interpolation='none', origin='upper', aspect='equal', rasterized=True)
        self._add_colorbar(im, self._finalize_heatmap_axis(h_ax, title=title))

    def plot_diff_matrix(self, diff, vmax, title=None):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        im = h_ax.imshow(diff, cmap=self.diff_cmap, vmax=vmax, vmin=-vmax, interpolation='none', origin='upper', aspect='equal', rasterized=True)
        self._add_colorbar(im, self._finalize_heatmap_axis(h_ax, title=title))

    def plot_diff_split_matrix(self, diff_upper, diff_lower, vmax, title=None):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        m_plot = np.triu(diff_upper.astype(float), k=1) + np.tril(diff_lower.astype(float), k=0)
        im = h_ax.imshow(m_plot, cmap=self.diff_cmap, vmax=vmax, vmin=-vmax, interpolation='none', origin='upper', aspect='equal', rasterized=True)
        self._add_colorbar(im, self._finalize_heatmap_axis(h_ax, title=title))

    def plot_coordinates(self, labelsize=10):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        if hasattr(self, 'heatmap_pos'): ax.set_position([self.heatmap_pos[0], ax.get_position().bounds[1], self.heatmap_pos[2], ax.get_position().bounds[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_xticks([0, self.n_bins])
        ax.set_xticklabels([f"{self.start:,}", f"{self.end:,}"], fontsize=labelsize)
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False, length=3)
        ax.set_yticks([])
        ax.text(-0.06, 0.5, self.chrom, fontsize=labelsize, va='center', ha='right', transform=ax.transAxes)

    def plot_motif_single_line(self, input_fil, track_name='Motif', label_size=10, marker_size=52):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        loci_pos, loci_neg = [], []
        chrom_variants = {self.chrom, self.chrom.lstrip('chr'), f"chr{self.chrom.lstrip('chr')}"}
        prefixes = tuple([c + '\t' for c in chrom_variants])

        try:
            with open(input_fil, 'r') as source:
                for line in source:
                    if not line.startswith(prefixes): continue 
                    parts = line.strip().split('\t')
                    if len(parts) < 3: continue
                    s, e = int(parts[1]), int(parts[2])
                    if s >= self.start and e <= self.end:
                        strand = parts[3] if len(parts) >= 4 else '+'
                        rel_pos = (((s + e) / 2) - self.start) / self.res
                        if strand == '+': loci_pos.append(rel_pos)
                        elif strand == '-': loci_neg.append(rel_pos)
        except Exception: pass

        if loci_pos: ax.scatter(loci_pos, [0]*len(loci_pos), s=marker_size, c='red', marker='>', edgecolors='none', clip_on=True)
        if loci_neg: ax.scatter(loci_neg, [0]*len(loci_neg), s=marker_size, c='blue', marker='<', edgecolors='none', clip_on=True)

        if hasattr(self, 'heatmap_pos'): ax.set_position([self.heatmap_pos[0], ax.get_position().bounds[1], self.heatmap_pos[2], ax.get_position().bounds[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_ylim(-1, 1)
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(-0.06, 0.5, track_name, transform=ax.transAxes, fontsize=label_size, ha='right', va='center')

    def plot_signal(self, track_name, bw_fil, color='#666666', label_size=10, nBins=1000):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        arr = np.array([])
        max_value = 10
        try:
            with pyBigWig.open(bw_fil) as db:
                c_fetch = next((v for v in [self.chrom, self.chrom.lstrip('chr'), f"chr{self.chrom.lstrip('chr')}"] if v in db.chroms()), None)
                if c_fetch: 
                    arr = np.array(db.stats(c_fetch, self.start, self.end, nBins=nBins, type="max")).astype(float)
                    raw_max = db.stats(c_fetch, self.start, self.end, type="max")[0]
                    max_value = math.ceil(raw_max / 10.0) * 10 if raw_max and raw_max > 0 else 10
        except Exception: pass

        arr[np.isnan(arr)] = 0
        ax.fill_between(np.linspace(0, self.n_bins, len(arr)), arr, color=color, edgecolor='none', linewidth=0)

        if hasattr(self, 'heatmap_pos'): ax.set_position([self.heatmap_pos[0], ax.get_position().bounds[1], self.heatmap_pos[2], ax.get_position().bounds[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_ylim(0, max_value)
        ax.axis('off')
        ax.text(-0.06, 0.5, track_name, transform=ax.transAxes, fontsize=label_size, ha='right', va='center')
        plot_y_axis(ax, max_value, label_size=label_size - 2)

    def outfig(self, outfile, dpi=300):
        self.fig.savefig(outfile, dpi=dpi, bbox_inches='tight')

if __name__ == '__main__':
    region_str, resolution_bp = parse_args()

    CHROM, START, END = parse_region(region_str)
    REGION_NAME = f"{CHROM}_{START}_{END}"
    OUT_DIR = os.path.join(PROJECT_ROOT, "enhanced_microc_fivepanel_for_simulation")
    os.makedirs(OUT_DIR, exist_ok=True)

    ENHANCE_BASE = f"{DATA_ROOT}/data_enhancement_1kb_cool_LY"
    cool_wt = os.path.join(ENHANCE_BASE, f"K562_MC_WT_1kb.cool_{resolution_bp}", "HiCFoundation_enhanced.cool")
    cool_ctcf = os.path.join(ENHANCE_BASE, f"K562_MC_CTCF_6h_9h_merge_007_026_1kb.cool_{resolution_bp}", "HiCFoundation_enhanced.cool")
    cool_rad21 = os.path.join(ENHANCE_BASE, f"K562_MC_RAD21_6h_9h_merge_024_031_1kb.cool_{resolution_bp}", "HiCFoundation_enhanced.cool")

    MOTIF_FILE = "/data/home/ruanlab/ly/K562_chromatin_folding/ChIA-PET/2.CTCF_motif_filter/K562_CTCF_ENCFF221SKA.motifs_p0.05.bed"
    BW_BASE = f"{DATA_ROOT}/mly/K562_CUT_TAG/4.IgG_signal_blocking"

    bw_files = {
        'CTCF': f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_CTCF_h019_153g_deep.RPGC.bw",
        'Cohesin': f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_RAD21_h020_160g_deep.RPGC.bw",
        'NIPBL': os.path.join(BW_BASE, "CTCF_0h_NIPBL_positive.bw"),
        'RNAPII': os.path.join(BW_BASE, "CTCF_0h_Pol2_positive.bw"),
        'H3K27ac': f"{DATA_ROOT}/mly/K562_CUT_TAG/3.loading_site_heatmap/K562_H3K27ac_ChIP-seq_ENCFF465GBD.bigWig",
        'H3K4me1': f"{DATA_ROOT}/mly/K562_CUT_TAG/3.loading_site_heatmap/K562_H3K4me1_ChIP-seq_ENCFF783QIW.bigWig",
    }
    track_cfg = [('CTCF', '#0000B2'), ('Cohesin', '#008000'), ('NIPBL', '#E69F00'), ('RNAPII', '#800080'), ('H3K27ac', '#CD3700'), ('H3K4me1', '#CD0000')]

    clr_wt, clr_ctcf, clr_rad21 = cooler.Cooler(cool_wt), cooler.Cooler(cool_ctcf), cooler.Cooler(cool_rad21)
    res = clr_wt.binsize
    if res != resolution_bp:
        print(f"[Warning] Requested resolution={resolution_bp}, but loaded binsize={res}")

    c_wt, c_ctcf, c_rad21 = find_chrom(clr_wt, CHROM), find_chrom(clr_ctcf, CHROM), find_chrom(clr_rad21, CHROM)
    m_wt = clr_wt.matrix(balance=False, sparse=False).fetch((c_wt, START, END))
    m_ctcf = clr_ctcf.matrix(balance=False, sparse=False).fetch((c_ctcf, START, END))
    m_rad21 = clr_rad21.matrix(balance=False, sparse=False).fetch((c_rad21, START, END))

    n_max = max(m_wt.shape[0], m_ctcf.shape[0], m_rad21.shape[0])
    m_wt, m_ctcf, m_rad21 = pad_matrix(np.nan_to_num(m_wt), n_max), pad_matrix(np.nan_to_num(m_ctcf), n_max), pad_matrix(np.nan_to_num(m_rad21), n_max)

    diff_rad21, diff_ctcf = m_rad21 - m_wt, m_ctcf - m_wt

    shared_raw_vmax = calc_shared_raw_vmax([m_wt, m_rad21, m_ctcf])
    shared_diff_vmax = calc_shared_diff_vmax([diff_rad21, diff_ctcf])

    track_ratios = [
        50, 4,
        50, 4,
        50, 4,
        50, 4,
        50, 6,
        3, 5
    ] + [5] * len(track_cfg)

    vis = FivePanelPlot(CHROM, START, END, n_max, res, figsize=(12, 32), n_rows=len(track_ratios), track_partition=track_ratios, space=0.01)

    vis.plot_split_matrix(m_wt, m_rad21, vmax=shared_raw_vmax, title="Enhanced Micro-C Split: WT (Upper) vs RAD21-AID (Lower)")
    vis.track_count += 1

    vis.plot_split_matrix(m_wt, m_ctcf, vmax=shared_raw_vmax, title="Enhanced Micro-C Split: WT (Upper) vs CTCF-AID (Lower)")
    vis.track_count += 1

    vis.plot_diff_matrix(diff_rad21, vmax=shared_diff_vmax, title="Enhanced Micro-C Difference: RAD21-AID - WT")
    vis.track_count += 1

    vis.plot_diff_matrix(diff_ctcf, vmax=shared_diff_vmax, title="Enhanced Micro-C Difference: CTCF-AID - WT")
    vis.track_count += 1

    vis.plot_diff_split_matrix(diff_upper=diff_rad21, diff_lower=diff_ctcf, vmax=shared_diff_vmax, title="Enhanced Micro-C Difference Split: RAD21-AID-WT (Upper) vs CTCF-AID-WT (Lower)")
    vis.track_count += 1
    
    vis.plot_coordinates(labelsize=10)
    vis.plot_motif_single_line(MOTIF_FILE, track_name='Motif', label_size=10)
    
    for track_label, color in track_cfg:
        vis.plot_signal(track_label, bw_files[track_label], color=color, label_size=10)

    vis.outfig(os.path.join(OUT_DIR, f"enhanced_MC_5panel_{res}bp_{REGION_NAME}.pdf"))