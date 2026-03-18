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
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize

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
                n_rows=13, track_partition=None, space=0.08, cmap_theme='red'):
        
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
        
        if cmap_theme == 'black':
            colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
            self.cmap = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
        else:
            self.cmap = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
            
        self.cmap.set_bad(color='white')

    def _find_chrom(self, clr, chrom):
        if chrom in clr.chromnames: return chrom
        if 'chr'+chrom in clr.chromnames: return 'chr'+chrom
        if chrom.startswith('chr') and chrom[3:] in clr.chromnames: return chrom[3:]
        return chrom

    def plot_matrix(self, vmax=None, vmin=None, title=None, mode='log'):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        M_to_plot = self.matrix.astype(float)
        
        if mode == 'log':
            M_to_plot[M_to_plot <= 0] = np.nan
            flat = M_to_plot[~np.isnan(M_to_plot)].flatten()
            if len(flat) > 0:
                if vmax is None: vmax = np.percentile(flat, 95)
                if vmin is None: vmin = vmax / 100.0
            else:
                vmax, vmin = 1, 0.01
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            if vmax is None: vmax = np.percentile(M_to_plot, 98)
            if vmin is None: vmin = 0
            norm = Normalize(vmin=vmin, vmax=vmax)

        im = h_ax.imshow(M_to_plot, cmap=self.cmap, norm=norm,
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

    def square_split_matrix_plot(self, vmax=None, vmin=None, title=None, mode='log'):
        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        M_wt = self.matrix.astype(float)
        M_ko = self.matrix_treated.astype(float)

        M_wt_part = np.triu(M_wt, k=1)
        M_ko_part = np.tril(M_ko, k=0)
        M_to_plot = M_wt_part + M_ko_part
        
        if mode == 'log':
            M_to_plot[M_to_plot <= 0] = np.nan
            if vmax is None:
                combined = np.concatenate((M_wt[~np.isnan(M_wt)], M_ko[~np.isnan(M_ko)]))
                if len(combined) > 0: vmax = np.percentile(combined, 95)
                else: vmax = 1
            if vmin is None: vmin = vmax / 100.0
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            if vmax is None:
                combined = np.concatenate((M_wt.flatten(), M_ko.flatten()))
                vmax = np.percentile(combined, 98)
            if vmin is None: vmin = 0
            norm = Normalize(vmin=vmin, vmax=vmax)

        im = h_ax.imshow(M_to_plot, cmap=self.cmap, norm=norm,
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

# =============================================================================
#  3. Main
# =============================================================================
if __name__ == '__main__':
    
    OUT_DIR = f"{PROJECT_ROOT}/examples"
    os.makedirs(OUT_DIR, exist_ok=True)

    region = {
        'name': 'chr8_104400000_104700000',
        'chrom': 'chr8',
        'start': 104480000,
        'end': 104720000
    }
    
    cool_files = [
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/K562_MC_WT_1kb.cool_2000/HiCFoundation_enhanced.cool", "Micro-C_WT"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/K562_MC_CTCF_6h_9h_merge_007_026_1kb.cool_2000/HiCFoundation_enhanced.cool", "Micro-C_CTCF_degron"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/K562_MC_RAD21_6h_9h_merge_024_031_1kb.cool_2000/HiCFoundation_enhanced.cool", "Micro-C_RAD21_degron"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/K562_CP_WT_RAD21_merge_1kb.cool_2000/HiCFoundation_enhanced.cool", "ChIA-PET_WT"),
        (f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/K562_CP_RAD21_6h_9h_merge_RAD21_1kb.cool_2000/HiCFoundation_enhanced.cool", "ChIA-PET_RAD21_degron")
    ]

    BASE_PATH = f'{DATA_ROOT}/mly/K562_CUT_TAG/4.IgG_signal_blocking'
    MOTIF_FILE = f'{PROJECT_ROOT}/ChIA-PET/2.CTCF_motif_filter/K562_CTCF_ENCFF221SKA.motifs_p0.05.bed'

    bw_files = {
        'CTCF':  f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_CTCF_h019_153g_deep.RPGC.bw",
        'RAD21': f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_RAD21_h020_160g_deep.RPGC.bw",
        'Cohesin': f"{DATA_ROOT}/mly/ChIA-PET/bw_all/WT_RAD21_h020_160g_deep.RPGC.bw",
        'SMC1':  os.path.join(BASE_PATH, 'CTCF_0h_SMC1_positive.bw'),
        'Pol2':  os.path.join(BASE_PATH, 'CTCF_0h_Pol2_positive.bw'),
        'RNAPII': os.path.join(BASE_PATH, 'CTCF_0h_Pol2_positive.bw'),
        'NIPBL': os.path.join(BASE_PATH, 'CTCF_0h_NIPBL_positive.bw'),
        'WAPL':  os.path.join(BASE_PATH, 'CTCF_0h_WAPL_positive.bw'),
        'YY1':   os.path.join(BASE_PATH, 'CTCF_0h_YY1_positive.bw'),
        'P300':  os.path.join(BASE_PATH, 'CTCF_0h_P300_positive.bw'),
        'BRD4':  os.path.join(BASE_PATH, 'CTCF_0h_BRD4_positive.bw'),
        'H3K27ac': f"{DATA_ROOT}/mly/K562_CUT_TAG/3.loading_site_heatmap/K562_H3K27ac_ChIP-seq_ENCFF465GBD.bigWig",
        'H3K4me1': f"{DATA_ROOT}/mly/K562_CUT_TAG/3.loading_site_heatmap/K562_H3K4me1_ChIP-seq_ENCFF783QIW.bigWig"
    }

    track_cfg = [
        ('CTCF', '#0000B2'), 
        ('Cohesin', '#008000'),
        ('NIPBL', '#E69F00'), 
        ('WAPL', '#8B0000'), 
        ('RNAPII', '#800080'), 
        ('H3K27ac', '#CD3700'),
        ('H3K4me1', '#CD0000')
    ]

    PLOT_MODE = 'linear' # Options: 'log', 'linear'
    COLOR_THEME = 'black'  # Options: 'red', 'black'

    unified_max_values = calculate_unified_max_values([region], bw_files)
    max_vals = unified_max_values[region['name']]

    color_scales = {}
    def get_group_name(label):
        if "Micro-C" in label: return "Micro-C"
        if "ChIA-PET" in label: return "ChIA-PET"
        return "Other"

    for cool_path, label in cool_files:
        print(f"Processing {label}...")
        
        n_tracks = len(track_cfg)
        track_ratios = [30, 1, 1.5] + [2] * n_tracks
        
        vis = intraChrom(
            uri_cool=cool_path,
            chrom=region['chrom'], 
            start=region['start'], 
            end=region['end'],
            figsize=(12, 14),
            n_rows=len(track_ratios), 
            track_partition=track_ratios,
            cmap_theme=COLOR_THEME
        )
        
        group = get_group_name(label)
        vmax_val, vmin_val = None, None

        if "WT" in label:
            M_tmp = vis.matrix.astype(float)
            if PLOT_MODE == 'log':
                M_tmp[M_tmp <= 0] = np.nan
                flat = M_tmp[~np.isnan(M_tmp)].flatten()
                if len(flat) > 0:
                    vmax_val = np.percentile(flat, 93)
                    vmin_val = vmax_val / 100.0
            else:
                flat = M_tmp.flatten()
                vmax_val = np.percentile(flat, 93)
                vmin_val = 0
            
            color_scales[group] = (vmax_val, vmin_val)

            print(f"  [Scale] Fixed {group} scale to vmax={vmax_val:.4f}, vmin={vmin_val:.4f}")
        
        if group in color_scales:
            vmax_val, vmin_val = color_scales[group]

        vis.plot_matrix(title=label, vmax=vmax_val, vmin=vmin_val, mode=PLOT_MODE)
        vis.plot_coordinates()
        vis.plot_elements(MOTIF_FILE, 'Motif')
        
        for track_label, color in track_cfg:
            key = track_label
            if key in bw_files:
                path = bw_files[key]
                limit = max_vals.get(key, 10)
                vis.plot_signal(track_label, path, color=color, max_value=limit)
            else:
                vis.track_count += 1
        
        out_f = os.path.join(OUT_DIR, f"log_style_{label}_{region['name']}.pdf")
        vis.outfig(out_f)
        plt.close(vis.fig)
        print(f"Saved: {out_f}")

    # Split Heatmaps
    file_map = {label: path for path, label in cool_files}
    comparison_pairs = [
        ("Micro-C_WT", "Micro-C_CTCF_degron"),
        ("Micro-C_WT", "Micro-C_RAD21_degron"),
        ("ChIA-PET_WT", "ChIA-PET_RAD21_degron")
    ]
    
    n_tracks = len(track_cfg)
    track_ratios = [30, 1, 1.5] + [2] * n_tracks

    for label_wt, label_ko in comparison_pairs:
        if label_wt not in file_map or label_ko not in file_map: continue
        
        print(f"Processing Split: {label_wt} vs {label_ko}...")
        
        vis = intraChrom(
            uri_cool=file_map[label_wt],
            chrom=region['chrom'], 
            start=region['start'], 
            end=region['end'],
            uri_treated=file_map[label_ko],
            figsize=(12, 14),
            n_rows=len(track_ratios), 
            track_partition=track_ratios,
            cmap_theme=COLOR_THEME
        )
        
        group = get_group_name(label_wt)
        vmax_val, vmin_val = color_scales.get(group, (None, None))
        
        vis.square_split_matrix_plot(title=f"{label_wt} (Upper) vs {label_ko} (Lower)", vmax=vmax_val, vmin=vmin_val, mode=PLOT_MODE)
        vis.plot_coordinates()
        vis.plot_elements(MOTIF_FILE, 'Motif')
        
        for track_label, color in track_cfg:
            key = track_label
            if key in bw_files:
                path = bw_files[key]
                limit = max_vals.get(key, 10)
                vis.plot_signal(track_label, path, color=color, max_value=limit)
            else:
                vis.track_count += 1
        
        out_f = os.path.join(OUT_DIR, f"split_log_{label_wt}_vs_{label_ko}_{region['name']}.pdf")
        vis.outfig(out_f)
        plt.close(vis.fig)
        print(f"Saved: {out_f}")

    print("All done.")