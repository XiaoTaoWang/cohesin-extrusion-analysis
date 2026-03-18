
# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
import os
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
from matplotlib.colors import LogNorm, LinearSegmentedColormap
import copy

# simple comment
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'

BED_FILE = f'{PROJECT_ROOT}/ADA/K562-loading-sites-within-convergent-loops.bed'
COOL_FILES = {
    'WT_0h': f'{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_CP_WT_RAD21_merge_1kb.mcool::/resolutions/2000',
    'Treated': f'{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool/K562_CP_RAD21_6h_9h_merge_RAD21_1kb.mcool::/resolutions/2000'
}

def process_bed(bed_path):
    bed = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name'])
    bed['loop_id'] = bed['name'].str.split('_').str[0]
    bed['type'] = bed['name'].str.split('_').str[1]
    
    triplets = []
    for loop_id, group in bed.groupby('loop_id'):
        s_row = group[group['type'] == 'S']
        e_row = group[group['type'] == 'E']
        m_rows = group[group['type'].str.startswith('M')]
        
        if s_row.empty or e_row.empty or m_rows.empty:
            continue
            
        s_pos = (s_row['start'].values[0] + s_row['end'].values[0]) // 2
        e_pos = (e_row['start'].values[0] + e_row['end'].values[0]) // 2
        chrom = s_row['chrom'].values[0]
        start_pos, end_pos = min(s_pos, e_pos), max(s_pos, e_pos)
        
        for _, m_row in m_rows.iterrows():
            m_pos = (m_row['start'] + m_row['end']) // 2
            if start_pos < m_pos < end_pos:
                triplets.append((chrom, start_pos, m_pos, end_pos))
    return triplets

def extract_and_scale_raw(clr, chrom, s_pos, m_pos, e_pos, K=30, flank_ratio=0.5):
    res = clr.binsize
    len_left = m_pos - s_pos
    len_right = e_pos - m_pos
    
    R0 = s_pos - len_left * flank_ratio
    R1 = s_pos
    R2 = m_pos
    R3 = e_pos
    R4 = e_pos + len_right * flank_ratio
    
    start_bin = int(R0 // res)
    end_bin = int(R4 // res) + 1
    
    try:
        raw_mat = clr.matrix(balance=True, sparse=False).fetch(f"{chrom}:{start_bin*res}-{end_bin*res}")
        raw_mat = np.nan_to_num(raw_mat)
    except Exception:
        return None
        
    flank_bins = int(K * flank_ratio)
    target_size = 2 * flank_bins + 2 * K
    
    target_edges = [0, flank_bins, flank_bins + K, flank_bins + 2 * K, target_size - 1]
    raw_edges = [(R0 - start_bin*res)/res, 
                 (R1 - start_bin*res)/res, 
                 (R2 - start_bin*res)/res, 
                 (R3 - start_bin*res)/res, 
                 (R4 - start_bin*res)/res]
                 
    target_coords = np.arange(target_size)
    map_1d = np.interp(target_coords, target_edges, raw_edges)
    X, Y = np.meshgrid(map_1d, map_1d, indexing='ij')
    
    scaled_mat = map_coordinates(raw_mat, [X, Y], order=1, mode='nearest')
    return scaled_mat

# --- Core workflow ---
triplets = process_bed(BED_FILE)
K_bins = 30
flank_r = 0.5

# 1. Compute and store aggregate matrices
agg_matrices = {}
for cond, cool_path in COOL_FILES.items():
    clr = cooler.Cooler(cool_path)
    mats = []
    for trip in triplets:
        mat = extract_and_scale_raw(clr, trip[0], trip[1], trip[2], trip[3], K=K_bins, flank_ratio=flank_r)
        if mat is not None:
            mats.append(mat)
    if mats:
        agg_matrices[cond] = np.nanmean(mats, axis=0)

# 2. Compute shared vmin/vmax (key fix: mask diagonals)
off_diag_values = []
for mat in agg_matrices.values():
    n = mat.shape[0]
    y, x = np.mgrid[0:n, 0:n]
    # Mask off-diagonal beyond |i-j| > 2 to remove dominant main diagonal
    mask = np.abs(y - x) > 2 
    off_diag_values.extend(mat[mask])

valid_values = np.array(off_diag_values)
valid_values = valid_values[valid_values > 0]

# After removing high-frequency noise, set vmin to 15th percentile to show weak stripes
global_vmin = np.percentile(valid_values, 15) 
global_vmax = np.percentile(valid_values, 99)

# ====== Core change: custom pure Hi-C red-white palette ======
# simple comment
custom_colors = ['white', '#FF0000'] 
cmap = LinearSegmentedColormap.from_list('juicebox_red', custom_colors)
cmap.set_under('white')
cmap.set_bad('white')
# =======================================================

# simple comment
fig, axes = plt.subplots(2, 1, figsize=(5, 8))
title_map = {'WT_0h': 'Wild Type', 'Treated': 'RAD21 depleted'}

for ax, (cond, agg_mat) in zip(axes, agg_matrices.items()):
    # simple comment
    im = ax.imshow(agg_mat, cmap=cmap, norm=LogNorm(vmin=global_vmin, vmax=global_vmax, clip=True), interpolation='none')
    ax.set_title(title_map.get(cond, cond), fontsize=16)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(2)

cbar_ax = fig.add_axes([0.05, 0.72, 0.03, 0.15]) 
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('Normalized\ncontact signals', fontsize=12, labelpad=10)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('right')

plt.subplots_adjust(hspace=0.2)
plt.savefig('Loading_Site_ADA_Clean_2kb_JuiceboxStyle.pdf', bbox_inches='tight')
plt.show()