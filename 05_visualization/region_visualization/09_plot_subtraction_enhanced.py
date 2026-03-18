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
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize, LogNorm

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
    except Exception as e:
        print(f"BW Read Error {os.path.basename(bw_path)}: {e}")
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

# =============================================================================
#  2. Plotting Class
# =============================================================================
class intraChrom(object):
    
    def __init__(self, uri_cool, chrom, start, end, uri_treated=None, figsize=(12, 14), 
                n_rows=13, track_partition=None, space=0.08):
        
        self.clr = cooler.Cooler(uri_cool)
        self.res = self.clr.binsize

        self.fig = plt.figure(figsize=figsize, facecolor='white')
        self.grid = GridSpec(n_rows, 1, figure=self.fig, left=0.15, right=0.92, 
                    bottom=0.03, top=0.97, hspace=space, height_ratios=track_partition)
        self.track_count = 0

        self.start = start
        self.end = end
        self.chrom = chrom
        
        c_name = self._find_chrom(self.clr, chrom)
        
        # Using balance=False as per original script
        try:
            M = self.clr.matrix(balance=False, sparse=False).fetch((c_name, start, end))
        except:
            M = self.clr.matrix(balance=False, sparse=False).fetch((c_name, start, end))
            
        M[np.isnan(M)] = 0
        self.matrix = M
        self.n_bins = M.shape[0]
        
        if uri_treated:
            self.clr_treated = cooler.Cooler(uri_treated)
            c_name_t = self._find_chrom(self.clr_treated, chrom)
            M_t = self.clr_treated.matrix(balance=False, sparse=False).fetch((c_name_t, start, end))
            M_t[np.isnan(M_t)] = 0
            
            n = max(self.n_bins, M_t.shape[0])
            if self.n_bins < n:
                self.matrix = np.pad(self.matrix, ((0,n-self.n_bins),(0,n-self.n_bins)), 'constant')
            if M_t.shape[0] < n:
                M_t = np.pad(M_t, ((0,n-M_t.shape[0]),(0,n-M_t.shape[0])), 'constant')
            
            self.matrix_treated = M_t
            self.n_bins = n

    def _find_chrom(self, clr, chrom):
        if chrom in clr.chromnames: return chrom
        if 'chr'+chrom in clr.chromnames: return 'chr'+chrom
        if chrom.startswith('chr') and chrom[3:] in clr.chromnames: return chrom[3:]
        return chrom

    def plot_difference_matrix(self, vmax=None, title=None, use_log=False):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        # Treated (matrix_treated) - Control (matrix)
        M_control = self.matrix.astype(float)
        M_treated = self.matrix_treated.astype(float)
        
    except AttributeError:
        raise AttributeError("`matrix_treated` is not available. Ensure `uri_treated` is provided during `intraChrom` initialization.")

        diff = M_treated - M_control

        if use_log:
            # Apply log transform, handling non-positive values
            diff_abs = np.abs(diff)
            diff_log = np.log1p(diff_abs)  # log(x+1) to handle 0 and small values
            
            # Determine scale based on logged values
            flat = diff_log.flatten()
            if len(flat) > 0:
                vmax = np.percentile(flat, 98) if vmax is None else vmax
            else:
                vmax = 1
            vmin = 0  # Log scale starts from 0
            
            cmap = plt.get_cmap('RdBu_r')  # You can adjust the colormap as needed
            cmap.set_bad('white')

            im = h_ax.imshow(diff_log, cmap=cmap, vmax=vmax, vmin=vmin,
                             interpolation='none', origin='upper', aspect='equal')
        else:
            # Determine scale
            if vmax is None:
                # Use 98th percentile of absolute values for robust scaling
                flat = np.abs(diff).flatten()
                if len(flat) > 0:
                    vmax = np.percentile(flat, 98)
                else:
                    vmax = 1
            
            vmin = -vmax
            
            cmap = plt.get_cmap('RdBu_r')
            cmap.set_bad('white')

            im = h_ax.imshow(diff, cmap=cmap, vmax=vmax, vmin=vmin,
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
        xticks = [0, self.n_bins]
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels, fontsize=labelsize)
        
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False, length=3)
        ax.set_yticks([])
        
        ax.text(-0.02, 0.5, self.chrom, fontsize=labelsize, weight='normal', va='center', ha='right', transform=ax.transAxes)

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
                        mid_bp = (s + e) / 2
                        relative_pos = (mid_bp - self.start) / self.res
                        
                        if strand == '+': loci_pos.append(relative_pos)
                        elif strand == '-': loci_neg.append(relative_pos)
        except Exception as e:
            print(f"Motif Read Warning: {e}")

        if loci_pos: 
            ax.scatter(loci_pos, [0]*len(loci_pos), s=marker_size, c='red', marker='>', edgecolors='none', clip_on=True)
        if loci_neg: 
            ax.scatter(loci_neg, [0]*len(loci_neg), s=marker_size, c='blue', marker='<', edgecolors='none', clip_on=True)
        
        if hasattr(self, 'heatmap_pos'):
            pos = ax.get_position().bounds
            ax.set_position([self.heatmap_pos[0], pos[1], self.heatmap_pos[2], pos[3]])

        ax.set_xlim(0, self.n_bins)
        ax.set_ylim(-1, 1)
        
        for spine in ax.spines: ax.spines[spine].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(track_name, rotation=0, fontsize=label_size, ha='right', va='center', labelpad=15)

    def plot_signal(self, track_name, bw_fil, color='#666666', max_value=None, 
                    label_size=10, nBins=1000):
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        arr = np.array([])
        try:
            with pyBigWig.open(bw_fil) as db:
                chrom_to_fetch = next((v for v in [self.chrom, self.chrom.lstrip('chr'), f"chr{self.chrom.lstrip('chr')}"] if v in db.chroms()), None)
                if chrom_to_fetch:
                    stats = db.stats(chrom_to_fetch, self.start, self.end, nBins=nBins, type="max")
                    arr = np.array(stats).astype(float)
        except Exception as e:
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
        ax.text(-0.02, 0.5, track_name, transform=ax.transAxes, 
                fontsize=label_size, ha='right', va='center', rotation=0)

        plot_y_axis(ax, max_value, label_size=label_size-2)

    def outfig(self, outfile, dpi=300):
        self.fig.savefig(outfile, dpi=dpi, bbox_inches='tight')
 # ============================================================================
 # 3. Main
 # ============================================================================
if __name__ == '__main__':
    
    OUT_DIR = f"{PROJECT_ROOT}/examples"
    # region = {
    #     'name': 'chr4_3753000_3787000',
    #     'chrom': 'chr4',
    #     'start': 37533000,
    #     'end': 37870000
    # }
    region = {
        'name': 'chr4_37580000_37820000',
        'chrom': 'chr4',
        'start': 37580000,
        'end': 37820000
    }
    cool_files = [
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_MC_WT_1kb.mcool::/resolutions/2000", "Micro-C_WT"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_MC_CTCF_6h_9h_merge_007_026_1kb.mcool::/resolutions/2000", "Micro-C_CTCF_degron"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_MC_RAD21_6h_9h_merge_024_031_1kb.mcool::/resolutions/2000", "Micro-C_RAD21_degron"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_CP_WT_RAD21_merge_1kb.mcool::/resolutions/2000", "ChIA-PET_WT"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_CP_RAD21_6h_9h_merge_RAD21_1kb.mcool::/resolutions/2000", "ChIA-PET_RAD21_degron")
    ]


    BASE_PATH = f'{DATA_ROOT}/mly/K562_CUT_TAG/4.IgG_signal_blocking'
    MOTIF_FILE = f'{PROJECT_ROOT}/ChIA-PET/2.CTCF_motif_filter/K562_CTCF_ENCFF221SKA.motifs_p0.05.bed'

    file_map = {label: path for path, label in cool_files}
    comparison_pairs = [
        ("Micro-C_WT", "Micro-C_CTCF_degron"),
        ("Micro-C_WT", "Micro-C_RAD21_degron"),
    ]

    for label_wt, label_ko in comparison_pairs:
        print(f"Processing Subtraction: {label_ko} - {label_wt} ...")
        vis = intraChrom(
            uri_cool=file_map[label_wt],
            chrom=region['chrom'],
            start=region['start'],
            end=region['end'],
        )
        vis.plot_difference_matrix(title=f"Subtraction: {label_ko} - {label_wt} (Linear)", use_log=False)
        vis.plot_coordinates()
        vis.plot_elements(MOTIF_FILE, 'Motif')
