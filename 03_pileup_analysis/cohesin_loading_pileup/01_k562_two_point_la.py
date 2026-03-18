# -*- coding: utf-8 -*-

# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
import os
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import os
import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from scipy.ndimage import zoom

from coolpuppy.coolpup import CoordCreator, PileUpper

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'


def prep_2point_dfs(bed_file, min_dist=10000, add_chr=True):
    df = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0,1,2,3], names=['chrom', 'start', 'end', 'name'])
    
    df['chrom'] = df['chrom'].astype(str)
    if add_chr:
        df['chrom'] = df['chrom'].apply(lambda x: x if x.startswith('chr') else 'chr' + x)
    else:
        df['chrom'] = df['chrom'].str.replace('^chr', '', regex=True)
        
    df['name'] = df['name'].str.rstrip('\r')
    df[['loop_id', 'site_type']] = df['name'].str.split('_', n=1, expand=True)
    
    s_dict = df[df['site_type'] == 'S'].set_index('loop_id')[['chrom', 'start']].to_dict('index')
    e_dict = df[df['site_type'] == 'E'].set_index('loop_id')[['chrom', 'end']].to_dict('index')
    
    res = {'S-E': [], 'S-L': [], 'L-E': []}
    l_df = df[df['site_type'].str.startswith('M-')]
    
    for _, row in l_df.iterrows():
        l_id = row['loop_id']
        if l_id in s_dict and l_id in e_dict:
            chrom = s_dict[l_id]['chrom']
            if chrom != e_dict[l_id]['chrom'] or chrom != row['chrom']: 
                continue
            
            s_pos = s_dict[l_id]['start']
            e_pos = e_dict[l_id]['end']
            l_pos = (row['start'] + row['end']) // 2 
            
            if s_pos > e_pos:
                s_pos, e_pos = e_pos, s_pos
            
            if (l_pos - s_pos) >= min_dist and (e_pos - l_pos) >= min_dist:
                res['S-E'].append({'chrom': chrom, 'start': s_pos, 'end': e_pos})
                res['S-L'].append({'chrom': chrom, 'start': s_pos, 'end': l_pos})
                res['L-E'].append({'chrom': chrom, 'start': l_pos, 'end': e_pos})
                
    return {k: pd.DataFrame(v) for k, v in res.items()}


class PileUpper2P(PileUpper):
    def _rescale_snip(self, snip):
        if snip["data"].size == 0 or np.all(np.isnan(snip["data"])):
            snip["data"] = np.full((self.rescale_size, self.rescale_size), np.nan)
            return snip

        raw_data = snip["data"]
        valid_mask = (~np.isnan(raw_data)).astype(float)
        data = np.nan_to_num(raw_data)
        res = self.resolution
        
        s_bp = snip["start1"]
        e_bp = snip["end1"]
        exp_start = snip["exp_start1"]

        bin_1 = int((s_bp - exp_start) // res)
        bin_2 = int((e_bp - exp_start) // res)
        bin_end = data.shape[0]

        bin_1 = max(0, min(bin_1, bin_end))
        bin_2 = max(0, min(bin_2, bin_end))
        
        safe_1 = min(bin_1, valid_mask.shape[0]-1)
        safe_2 = min(bin_2, valid_mask.shape[0]-1)
        if valid_mask[safe_1, safe_1] == 0 or valid_mask[safe_2, safe_2] == 0:
            snip["data"] = np.full((self.rescale_size, self.rescale_size), np.nan)
            return snip

        bins = [0, bin_1, bin_2, bin_end]
        
        N = self.rescale_size 
        pad_frac = self.rescale_flank 
        total_len = 2 * pad_frac + 1.0
        
        pix_1 = int(N * pad_frac / total_len)
        pix_2 = int(N * (pad_frac + 1.0) / total_len)
        pixs = [0, pix_1, pix_2, N]
        
        row_blocks = []
        valid_row_blocks = []
        for i in range(3):
            col_blocks = []
            valid_col_blocks = []
            h_target = pixs[i+1] - pixs[i]
            
            for j in range(3):
                w_target = pixs[j+1] - pixs[j]
                target_shape = (h_target, w_target)
                
                q = data[bins[i]:bins[i+1], bins[j]:bins[j+1]]
                q_valid = valid_mask[bins[i]:bins[i+1], bins[j]:bins[j+1]]
                
                if h_target == 0 or w_target == 0 or q.shape[0] == 0 or q.shape[1] == 0:
                    r = np.zeros(target_shape)
                    r_val = np.zeros(target_shape)
                else:
                    zoom_factors = (h_target / q.shape[0], w_target / q.shape[1])
                    r = zoom(q, zoom_factors, order=1)
                    r_val = zoom(q_valid, zoom_factors, order=1)
                
                col_blocks.append(r)
                valid_col_blocks.append(r_val)
                
            row_blocks.append(np.hstack(col_blocks))
            valid_row_blocks.append(np.hstack(valid_col_blocks))
            
        final_data = np.vstack(row_blocks)
        final_valid = np.vstack(valid_row_blocks)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            final_data = final_data / (final_valid + 1e-12)
        
        final_data[final_valid < 0.5] = np.nan
        snip["data"] = final_data

        if self.coverage_norm:
            snip["cov_start"] = np.zeros(N)
            snip["cov_end"] = np.zeros(N)

        return snip


def run_2point_pileup():
    bed_file = f"{PROJECT_ROOT}/ADA/K562-loading-sites-within-convergent-loops.bed"
    MCOOL_DIR = f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool"
    
    mcool_files = [
        os.path.join(MCOOL_DIR, "K562_MC_WT_1kb.mcool::/resolutions/2000"),
        os.path.join(MCOOL_DIR, "K562_MC_CTCF_6h_9h_merge_007_026_1kb.mcool::/resolutions/2000"),
        os.path.join(MCOOL_DIR, "K562_MC_RAD21_6h_9h_merge_024_031_1kb.mcool::/resolutions/2000"),
    ]
    
    procs = 10
    flank_ratio = 0.5
    min_dist_threshold = 20000 
    cmap_theme = 'black' 
    
    # custom_vmin = 0.0001
    # custom_vmax = 0.005
    custom_vmin = 0.5
    custom_vmax = 6
    

    if cmap_theme == 'black':
        colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
        custom_cmap = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
    else:
        custom_cmap = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
    custom_cmap.set_bad('#FFFFFF')

    def plot_and_save(matrix, out_png, out_txt, labels, title_text):
        np.savetxt(out_txt, matrix, delimiter='\t', fmt='%.6e')
        fig, ax = plt.subplots(figsize=(6, 5))
        valid_vals = matrix[~np.isnan(matrix) & (matrix > 0)]
        
        vmax = custom_vmax if custom_vmax is not None else (np.percentile(valid_vals, 80) if len(valid_vals) > 0 else 1.0)
        vmin = custom_vmin if custom_vmin is not None else (np.percentile(valid_vals, 5) if len(valid_vals) > 0 else 0.0)

        im = ax.imshow(matrix, cmap=custom_cmap, norm=Normalize(vmin=vmin, vmax=vmax))
        
        N = matrix.shape[0]
        total_len = 2 * flank_ratio + 1.0
        pix_1 = int(N * flank_ratio / total_len)         
        pix_2 = int(N * (flank_ratio + 1.0) / total_len) 
        
        ax.set_xticks([pix_1, pix_2])
        ax.set_xticklabels(labels)
        ax.set_yticks([pix_1, pix_2])
        ax.set_yticklabels(labels)
        
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        ax.set_title(title_text, fontsize=10)
        plt.colorbar(im, label='Raw Contact Frequency')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()

    for cool_file in mcool_files:
        base_name = os.path.basename(cool_file.split('::')[0]).replace('.mcool', '')
        clr = cooler.Cooler(cool_file)
        has_chr_prefix = any(chrom.startswith('chr') for chrom in clr.chromnames)
        
        dfs = prep_2point_dfs(bed_file, min_dist=min_dist_threshold, add_chr=has_chr_prefix)
        matrices = {}
        
        for pair_type, df in dfs.items():
            retained_count = len(df)
            print(f"[{base_name}] {pair_type} interactions: {retained_count}")
            if retained_count == 0:
                continue
                
            CC = CoordCreator(
                features=df, resolution=clr.binsize, features_format="bed",
                local=True, nshifts=0, rescale_flank=flank_ratio, 
            )

            PU = PileUpper2P(
                clr=clr, CC=CC, clr_weight_name=None, rescale=True,
                rescale_size=99, ignore_diags=0, nproc=procs,
            )
            
            pups_df = PU.pileupsWithControl(nproc=procs)
            matrix = pups_df['data'].values[0]
            matrices[pair_type] = matrix
            
            out_png = f"coolpuppy_{pair_type}_{base_name}_2K.pdf"
            out_txt = f"coolpuppy_{pair_type}_{base_name}_2K.txt"
            labels = pair_type.split('-')
            plot_and_save(matrix, out_png, out_txt, labels, f"{base_name} ({pair_type})\nRetained: {retained_count}")

        if 'S-L' in matrices and 'L-E' in matrices:
            print(f"[{base_name}] Merging S-L and mirrored L-E...")
            matrix_sl = matrices['S-L']
            matrix_le = matrices['L-E']
            
            matrix_le_mirrored = np.flip(matrix_le).T
            
            with np.errstate(invalid='ignore'):
                merged_matrix = np.nanmean([matrix_sl, matrix_le_mirrored], axis=0)
            
            out_png_merged = f"coolpuppy_Merge_SL_LE_{base_name}_2K.pdf"
            out_txt_merged = f"coolpuppy_Merge_SL_LE_{base_name}_2K.txt"
            
            labels_merged = ['S/E', 'L']
            title_merged = f"{base_name} (Merged S-L & mirrored L-E)\nRetained: 160"
            plot_and_save(merged_matrix, out_png_merged, out_txt_merged, labels_merged, title_merged)

if __name__ == "__main__":
    run_2point_pileup()