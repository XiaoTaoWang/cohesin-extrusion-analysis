#!/usr/bin/env python
# -*- coding: utf-8 -*-


# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
import os
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import os
import math
import cooler
import pyBigWig
import numpy as np
import collections
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
import concurrent.futures
from functools import partial

# =============================================================================
#  0. Global Settings
# =============================================================================
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

# =============================================================================
#  1. Helper Functions
# =============================================================================
def plot_y_axis(ax, max_val, label_size=10):
    max_str = str(int(max_val)) if max_val % 1 == 0 else f"{max_val:.1f}"
    ax_pos = ax.get_position().bounds
    gap = 0.015 
    width = 0.005
    y_ax = ax.figure.add_axes([ax_pos[0] - gap, ax_pos[1], width, ax_pos[3]])
    y_ax.plot([1, 1], [0, 1], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.plot([0, 1], [0, 0], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.plot([0, 1], [1, 1], color='black', linewidth=1.0, transform=y_ax.transAxes, clip_on=False)
    y_ax.text(-0.5, 0, "0", ha='right', va='center', fontsize=label_size, transform=y_ax.transAxes)
    y_ax.text(-0.5, 1, max_str, ha='right', va='center', fontsize=label_size, transform=y_ax.transAxes)
    y_ax.axis('off')

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

def calculate_unified_max_values(loading_sites, bw_files):
    print("="*50 + "\nPre-calculating Y-axis max values...")
    unified_maxes = {}
    for site in loading_sites:
        name = site['name']
        unified_maxes[name] = {}
        for protein, path in bw_files.items():
            raw = get_raw_max_value(path, site['chrom'], site['start'], site['end'])
            unified_maxes[name][protein] = math.ceil(raw / 10.0) * 10 if raw > 0 else 10
    print("Pre-calculation done!\n" + "="*50)
    return unified_maxes

def normalize_diag(m_wt, m_ko):
    """
    通过对比0h(wt)和6h(ko)主对角线的总和，计算系数并对ko矩阵进行缩放对齐。
    返回: 缩放后的矩阵, 缩放系数
    """
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

# =============================================================================
#  2. Plotting Class
# =============================================================================
class intraChrom(object):
    
    def __init__(self, uri_cool, chrom, start, end, uri_treated=None, figsize=(12, 14), 
                n_rows=13, track_partition=None, space=0.08, matrix=None, matrix_treated=None):
        
        self.clr = cooler.Cooler(uri_cool)
        self.res = self.clr.binsize
        self.fig = plt.figure(figsize=figsize, facecolor='white')
        self.grid = GridSpec(n_rows, 1, figure=self.fig, left=0.15, right=0.92, 
                    bottom=0.03, top=0.97, hspace=space, height_ratios=track_partition)
        self.track_count = 0
        self.start = start
        self.end = end
        self.chrom = chrom
        
        if matrix is not None:
            self.matrix = matrix
            self.n_bins = matrix.shape[0]
            
        if matrix_treated is not None:
            self.matrix_treated = matrix_treated
            self.n_bins = max(self.n_bins, matrix_treated.shape[0])
            
        colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
        self.cmap_black = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
        self.cmap_black.set_bad(color='white')
        self.cmap_diff = plt.get_cmap('RdYlBu_r')
        self.cmap_diff.set_bad(color='white')

    def square_split_matrix_plot(self, vmax=None, vmin=None, title=None, mode='linear'):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        M_wt = self.matrix.astype(float)
        M_ko = self.matrix_treated.astype(float)

        M_wt_part = np.triu(M_wt, k=1)
        M_ko_part = np.tril(M_ko, k=0)
        M_to_plot = M_wt_part + M_ko_part
        M_to_plot[M_to_plot <= 0] = np.nan

        if vmax is None:
            valid_vals = M_to_plot[~np.isnan(M_to_plot)]
            vmax = np.percentile(valid_vals, 95) if len(valid_vals) > 0 else 1

        if mode == 'log':
            if vmin is None: vmin = vmax / 50.0 if vmax > 0 else 0.01
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            if vmin is None: vmin = 0
            norm = Normalize(vmin=vmin, vmax=vmax)

        im = h_ax.imshow(M_to_plot, cmap=self.cmap_black, norm=norm,
                         interpolation='none', origin='upper', aspect='equal')
        h_ax.axis('off')
        
        hb = h_ax.get_position().bounds
        w, h = hb[2], hb[3]
        size = min(w, h)
        h_ax.set_position([hb[0] + (w-size)/2, hb[1] + (h-size)/2, size, size])
        self.heatmap_pos = h_ax.get_position().bounds 

        if title:
            h_ax.set_title(title, fontsize=14, pad=10)

        cbar_ax = self.fig.add_axes([self.heatmap_pos[0] + size + 0.01, self.heatmap_pos[1], 0.02, size * 0.2])
        self.fig.colorbar(im, cax=cbar_ax)
        cbar_ax.tick_params(labelsize=9)

    def plot_difference_matrix(self, vmax=None, title=None):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        M_control = self.matrix.astype(float)
        M_treated = self.matrix_treated.astype(float)
        diff = M_treated - M_control
        
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
        h_ax.set_position([hb[0] + (w-size)/2, hb[1] + (h-size)/2, size, size])
        self.heatmap_pos = h_ax.get_position().bounds 

        if title:
            h_ax.set_title(title, fontsize=14, pad=10)

        cbar_ax = self.fig.add_axes([self.heatmap_pos[0] + size + 0.01, self.heatmap_pos[1], 0.02, size * 0.2])
        self.fig.colorbar(im, cax=cbar_ax)
        cbar_ax.tick_params(labelsize=9)

    def square_split_diff_matrix_plot(self, diff_lower, diff_upper, vmax=None, title=None):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        M_plot = np.tril(diff_lower, k=0) + np.triu(diff_upper, k=1)

        if vmax is None:
            flat = np.abs(M_plot).flatten()
            flat_valid = flat[flat > 0]
            vmax = np.percentile(flat_valid, 95) if len(flat_valid) > 0 else 1

        im = h_ax.imshow(M_plot, cmap=self.cmap_diff, vmax=vmax, vmin=-vmax,
                         interpolation='none', origin='upper', aspect='equal')
        h_ax.axis('off')
        
        hb = h_ax.get_position().bounds
        w, h = hb[2], hb[3]
        size = min(w, h)
        h_ax.set_position([hb[0] + (w-size)/2, hb[1] + (h-size)/2, size, size])
        self.heatmap_pos = h_ax.get_position().bounds 

        if title:
            h_ax.set_title(title, fontsize=14, pad=10)

        cbar_ax = self.fig.add_axes([self.heatmap_pos[0] + size + 0.01, self.heatmap_pos[1], 0.02, size * 0.2])
        self.fig.colorbar(im, cax=cbar_ax)
        cbar_ax.tick_params(labelsize=9)

    def plot_coordinates(self, labelsize=10):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        if hasattr(self, 'heatmap_pos'):
            pos = ax.get_position().bounds
            ax.set_position([self.heatmap_pos[0], pos[1], self.heatmap_pos[2], pos[3]])
        ax.set_xlim(0, self.n_bins)
        labels = [f"{self.start:,}", f"{self.end:,}"]
        ax.set_xticks([0, self.n_bins])
        ax.set_xticklabels(labels, fontsize=labelsize)
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False, length=3)
        ax.set_yticks([])
        ax.text(-0.02, 0.5, self.chrom, fontsize=labelsize, va='center', ha='right', transform=ax.transAxes)

    def plot_elements(self, input_fil, track_name, label_size=10, marker_size=60):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        loci_pos, loci_neg = [], []
        chrom_variants = {self.chrom, self.chrom.lstrip('chr'), f"chr{self.chrom.lstrip('chr')}"}
        try:
            with open(input_fil, 'r') as source:
                for line in source:
                    parts = line.rstrip().split()
                    if len(parts) < 3: continue
                    c, s, e = parts[0], int(parts[1]), int(parts[2])
                    strand = parts[3] if len(parts) >= 4 else '+'
                    if (c in chrom_variants) and (s >= self.start) and (e <= self.end):
                        relative_pos = (((s + e) / 2) - self.start) / self.res
                        if strand == '+': loci_pos.append(relative_pos)
                        elif strand == '-': loci_neg.append(relative_pos)
        except Exception:
            pass
        if loci_pos: ax.scatter(loci_pos, [0]*len(loci_pos), s=marker_size, c='red', marker='>')
        if loci_neg: ax.scatter(loci_neg, [0]*len(loci_neg), s=marker_size, c='blue', marker='<')
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
                chrom_to_fetch = next((v for v in [self.chrom, self.chrom.lstrip('chr'), f"chr{self.chrom.lstrip('chr')}"] if v in db.chroms()), None)
                if chrom_to_fetch:
                    stats = db.stats(chrom_to_fetch, self.start, self.end, nBins=nBins, type="max")
                    arr = np.array(stats).astype(float)
        except Exception:
            pass
        arr[np.isnan(arr)] = 0
        if max_value is None: max_value = np.max(arr) if arr.size > 0 else 1
        x = np.linspace(0, self.n_bins, len(arr))
        ax.fill_between(x, arr, color=color, edgecolor='none', linewidth=0)
        if hasattr(self, 'heatmap_pos'):
            pos = ax.get_position().bounds
            ax.set_position([self.heatmap_pos[0], pos[1], self.heatmap_pos[2], pos[3]])
        ax.set_xlim(0, self.n_bins)
        ax.set_ylim(0, max_value)
        ax.axis('off')
        ax.text(-0.02, 0.5, track_name, transform=ax.transAxes, fontsize=label_size, ha='right', va='center')
        plot_y_axis(ax, max_value, label_size=label_size-2)

    def outfig(self, outfile, dpi=300):
        self.fig.savefig(outfile, dpi=dpi, bbox_inches='tight')

# =============================================================================
#  3. Task Process Function (Multiprocessing safe)
# =============================================================================
def process_region_task(region, hic_files, bw_files, track_cfg, unified_max_values, out_dir, motif_file):
    max_vals = unified_max_values[region['name']]
    track_ratios = [30, 1, 1.5] + [2] * len(track_cfg)
    res = 1000

    try:
        c_ctcf_0h = cooler.Cooler(hic_files["CTCF_0h"])
        c_ctcf_6h = cooler.Cooler(hic_files["CTCF_6h"])
        c_rad21_0h = cooler.Cooler(hic_files["RAD21_0h"])
        c_rad21_6h = cooler.Cooler(hic_files["RAD21_6h"])
        
        m_ctcf_0h = c_ctcf_0h.matrix(balance=False, sparse=False).fetch((region['chrom'], region['start'], region['end']))
        m_ctcf_6h_raw = c_ctcf_6h.matrix(balance=False, sparse=False).fetch((region['chrom'], region['start'], region['end']))
        m_rad21_0h = c_rad21_0h.matrix(balance=False, sparse=False).fetch((region['chrom'], region['start'], region['end']))
        m_rad21_6h_raw = c_rad21_6h.matrix(balance=False, sparse=False).fetch((region['chrom'], region['start'], region['end']))
        
        m_ctcf_0h[np.isnan(m_ctcf_0h)] = 0
        m_ctcf_6h_raw[np.isnan(m_ctcf_6h_raw)] = 0
        m_rad21_0h[np.isnan(m_rad21_0h)] = 0
        m_rad21_6h_raw[np.isnan(m_rad21_6h_raw)] = 0
        
        n_max = max(m_ctcf_0h.shape[0], m_ctcf_6h_raw.shape[0], m_rad21_0h.shape[0], m_rad21_6h_raw.shape[0])
        m_ctcf_0h = pad_matrix(m_ctcf_0h, n_max)
        m_ctcf_6h_raw = pad_matrix(m_ctcf_6h_raw, n_max)
        m_rad21_0h = pad_matrix(m_rad21_0h, n_max)
        m_rad21_6h_raw = pad_matrix(m_rad21_6h_raw, n_max)

        # simple comment
        m_ctcf_6h_norm, coef_ctcf = normalize_diag(m_ctcf_0h, m_ctcf_6h_raw)
        m_rad21_6h_norm, coef_rad21 = normalize_diag(m_rad21_0h, m_rad21_6h_raw)
        
        # simple comment
        vals_ctcf_all = np.concatenate([m_ctcf_0h[m_ctcf_0h>0], m_ctcf_6h_norm[m_ctcf_6h_norm>0]])
        vmax_split_ctcf = np.percentile(vals_ctcf_all, 95) if len(vals_ctcf_all) > 0 else 1

        # simple comment
        vals_rad21_all = np.concatenate([m_rad21_0h[m_rad21_0h>0], m_rad21_6h_norm[m_rad21_6h_norm>0]])
        vmax_split_rad21 = np.percentile(vals_rad21_all, 95) if len(vals_rad21_all) > 0 else 1
        
        # simple comment
        diff_ctcf = m_ctcf_6h_norm - m_ctcf_0h
        diff_rad21 = m_rad21_6h_norm - m_rad21_0h
        
        all_diff_vals = np.concatenate([np.abs(diff_ctcf)[np.abs(diff_ctcf)>0], 
                                        np.abs(diff_rad21)[np.abs(diff_rad21)>0]])
        shared_diff_vmax = np.percentile(all_diff_vals, 95) if len(all_diff_vals) > 0 else 1
        
    except Exception as e:
        return f"Skipped {region['name']} at {res}bp: {e}"

    # simple comment
    # simple comment
    vis_split_ctcf = intraChrom(uri_cool=hic_files["CTCF_0h"], chrom=region['chrom'], start=region['start'], end=region['end'],
                                uri_treated=hic_files["CTCF_6h"], figsize=(12, 14), n_rows=len(track_ratios), track_partition=track_ratios,
                                matrix=m_ctcf_0h, matrix_treated=m_ctcf_6h_norm)
    vis_split_ctcf.square_split_matrix_plot(vmax=vmax_split_ctcf, title=f"Split enhanced 1kb: WT (Up) vs CTCF 6h (Dn) [x{coef_ctcf:.2f}]", mode="linear")
    vis_split_ctcf.plot_coordinates()
    vis_split_ctcf.plot_elements(motif_file, 'Motif')
    for track_label, color in track_cfg:
        if track_label in bw_files: vis_split_ctcf.plot_signal(track_label, bw_files[track_label], color=color, max_value=max_vals.get(track_label, 10))
        else: vis_split_ctcf.track_count += 1
    vis_split_ctcf.outfig(os.path.join(out_dir, f"enhanced_split_CTCF_{res}bp_{region['name']}.pdf"))
    plt.close(vis_split_ctcf.fig)

    vis_diff_ctcf = intraChrom(uri_cool=hic_files["CTCF_0h"], chrom=region['chrom'], start=region['start'], end=region['end'],
                               uri_treated=hic_files["CTCF_6h"], figsize=(12, 14), n_rows=len(track_ratios), track_partition=track_ratios,
                               matrix=m_ctcf_0h, matrix_treated=m_ctcf_6h_norm)
    vis_diff_ctcf.plot_difference_matrix(vmax=shared_diff_vmax, title=f"Diff enhanced 1kb: CTCF 6h - 0h")
    vis_diff_ctcf.plot_coordinates()
    vis_diff_ctcf.plot_elements(motif_file, 'Motif')
    for track_label, color in track_cfg:
        if track_label in bw_files: vis_diff_ctcf.plot_signal(track_label, bw_files[track_label], color=color, max_value=max_vals.get(track_label, 10))
        else: vis_diff_ctcf.track_count += 1
    vis_diff_ctcf.outfig(os.path.join(out_dir, f"enhanced_diff_CTCF_{res}bp_{region['name']}.pdf"))
    plt.close(vis_diff_ctcf.fig)

    # simple comment
    # simple comment
    vis_split_rad = intraChrom(uri_cool=hic_files["RAD21_0h"], chrom=region['chrom'], start=region['start'], end=region['end'],
                               uri_treated=hic_files["RAD21_6h"], figsize=(12, 14), n_rows=len(track_ratios), track_partition=track_ratios,
                               matrix=m_rad21_0h, matrix_treated=m_rad21_6h_norm)
    vis_split_rad.square_split_matrix_plot(vmax=vmax_split_rad21, title=f"Split enhanced 1kb: WT (Up) vs RAD21 6h (Dn) [x{coef_rad21:.2f}]", mode="linear")
    vis_split_rad.plot_coordinates()
    vis_split_rad.plot_elements(motif_file, 'Motif')
    for track_label, color in track_cfg:
        if track_label in bw_files: vis_split_rad.plot_signal(track_label, bw_files[track_label], color=color, max_value=max_vals.get(track_label, 10))
        else: vis_split_rad.track_count += 1
    vis_split_rad.outfig(os.path.join(out_dir, f"enhanced_split_RAD21_{res}bp_{region['name']}.pdf"))
    plt.close(vis_split_rad.fig)

    vis_diff_rad = intraChrom(uri_cool=hic_files["RAD21_0h"], chrom=region['chrom'], start=region['start'], end=region['end'],
                              uri_treated=hic_files["RAD21_6h"], figsize=(12, 14), n_rows=len(track_ratios), track_partition=track_ratios,
                              matrix=m_rad21_0h, matrix_treated=m_rad21_6h_norm)
    vis_diff_rad.plot_difference_matrix(vmax=shared_diff_vmax, title=f"Diff enhanced 1kb: RAD21 6h - 0h")
    vis_diff_rad.plot_coordinates()
    vis_diff_rad.plot_elements(motif_file, 'Motif')
    for track_label, color in track_cfg:
        if track_label in bw_files: vis_diff_rad.plot_signal(track_label, bw_files[track_label], color=color, max_value=max_vals.get(track_label, 10))
        else: vis_diff_rad.track_count += 1
    vis_diff_rad.outfig(os.path.join(out_dir, f"enhanced_diff_RAD21_{res}bp_{region['name']}.pdf"))
    plt.close(vis_diff_rad.fig)

    # simple comment
    vis_split_diff = intraChrom(uri_cool=hic_files["RAD21_0h"], chrom=region['chrom'], start=region['start'], end=region['end'],
                                figsize=(12, 14), n_rows=len(track_ratios), track_partition=track_ratios, matrix=m_rad21_0h)
    title_split_diff = f"Split Diff enhanced 1kb: RAD21 6h-0h (Up) vs CTCF 6h-0h (Dn)"
    vis_split_diff.square_split_diff_matrix_plot(diff_ctcf, diff_rad21, vmax=shared_diff_vmax, title=title_split_diff)
    vis_split_diff.plot_coordinates()
    vis_split_diff.plot_elements(motif_file, 'Motif')
    for track_label, color in track_cfg:
        if track_label in bw_files: vis_split_diff.plot_signal(track_label, bw_files[track_label], color=color, max_value=max_vals.get(track_label, 10))
        else: vis_split_diff.track_count += 1
    vis_split_diff.outfig(os.path.join(out_dir, f"enhanced_split_diff_CTCF_vs_RAD21_{res}bp_{region['name']}.pdf"))
    plt.close(vis_split_diff.fig)
        
    return f"Successfully processed {region['name']} at {res}bp"

# =============================================================================
#  4. Main Pipeline
# =============================================================================
if __name__ == '__main__':
    
    OUT_DIR = f"{PROJECT_ROOT}/diag_norm_HCT116_version2_split"
    os.makedirs(OUT_DIR, exist_ok=True)
    
    hic_base = f"{DATA_ROOT}/single_molecule/HCT116_intact_hic_enhance"
    hic_files = {
        "CTCF_0h": os.path.join(hic_base, "CTCF_0h_1kb/HiCFoundation_enhanced.cool"),
        "CTCF_6h": os.path.join(hic_base, "CTCF_6h_1kb/HiCFoundation_enhanced.cool"),
        "RAD21_0h": os.path.join(hic_base, "RAD21_0h_1kb/HiCFoundation_enhanced.cool"),
        "RAD21_6h": os.path.join(hic_base, "RAD21_6h_1kb/HiCFoundation_enhanced.cool")
    }

    bed_file = f"{PROJECT_ROOT}/script_minji/20260213_from_minji/HCT116-cohesin-specific-regions_20230710_filt_25-75percentile_20231030.bed"
    MOTIF_FILE = f'{PROJECT_ROOT}/ChIA-PET/2.CTCF_motif_filter/K562_CTCF_ENCFF221SKA.motifs_p0.05.bed'
    EXTENSION = 10000 
    NUM_WORKERS = 10
    
    bw_files_0h = {
        'CTCF':    f"{DATA_ROOT}/ChIA-PET/HCT116-CTCF-ChIA-PET/HCT116-CTCF-ChIA-PET.forBASIC.RPGC.bw",
        'RAD21':   f"{DATA_ROOT}/ChIA-PET/HCT116-RAD21-ChIA-PET/HCT116-RAD21-ChIA-PET.forBASIC.RPGC.bw",
        'NIPBL':   f"{PROJECT_ROOT}/HCT116/4DNFI84R4CIL.bw",
        'RNAPII':  f"{DATA_ROOT}/ChIA-PET/HCT116-RNAPII-ChIA-PET/HCT116-RNAPII-ChIA-PET.forBASIC.RPGC.bw",
        'H3K27ac': f"{PROJECT_ROOT}/HCT116/public_data/HCT116_H3K27ac_chip-seq_4DNFII8WPH9V.bw",
        'H3K4me1': f"{PROJECT_ROOT}/HCT116/public_data/HCT116_H3K4me1_chip-seq_4DNFIRNLICB9.bw"
    }

    track_cfg = [
        ('CTCF', '#0000B2'), 
        ('RAD21', '#008000'),
        ('NIPBL', '#E69F00'), 
        ('RNAPII', '#800080'), 
        ('H3K27ac', '#CD3700'),
        ('H3K4me1', '#CD0000')
    ]

    regions_dict = collections.defaultdict(dict)
    
    if os.path.exists(bed_file):
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    chrom, s, e, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    
                    if '_' not in name: continue
                    base_name, suffix = name.split('_', 1)
                    
                    regions_dict[base_name]['chrom'] = chrom
                    if suffix == 'S':
                        regions_dict[base_name]['start'] = s
                    elif suffix == 'E':
                        regions_dict[base_name]['end'] = e
    else:
        print(f"Error: {bed_file} not found.")
        exit(1)

    regions = []
    for base_name, data in regions_dict.items():
        if 'start' in data and 'end' in data:
            regions.append({
                'name': base_name,
                'chrom': data['chrom'],
                'start': max(0, data['start'] - EXTENSION),
                'end': data['end'] + EXTENSION
            })
            
    if not regions:
        print("Error: No valid S and E pairs found in BED file.")
        exit(1)

    unified_max_values = calculate_unified_max_values(regions, bw_files_0h)

    print(f"Starting multiprocessing with {NUM_WORKERS} workers...")
    
    worker = partial(
        process_region_task,
        hic_files=hic_files,
        bw_files=bw_files_0h,
        track_cfg=track_cfg,
        unified_max_values=unified_max_values,
        out_dir=OUT_DIR,
        motif_file=MOTIF_FILE
    )

    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        results = executor.map(worker, regions)
        for result in results:
            if result:
                print(result)

    print("All tasks completed.")