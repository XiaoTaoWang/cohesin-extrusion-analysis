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


def prep_3point_df(bed_file, min_la_dist=10000, add_chr=True):
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
    
    res = []
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
            
            dist_SL = l_pos - s_pos
            dist_LE = e_pos - l_pos
            
            if dist_SL >= min_la_dist and dist_LE >= min_la_dist:
                res.append({
                    'chrom': chrom,
                    'start': s_pos,
                    'end': e_pos,
                    'pos_L': l_pos
                })
    return pd.DataFrame(res)


class PileUpper3P(PileUpper):
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
        l_bp = snip["pos_L1"]
        exp_start = snip["exp_start1"]

        bin_S = int((s_bp - exp_start) // res)
        bin_L = int((l_bp - exp_start) // res)
        bin_E = int((e_bp - exp_start) // res)
        bin_end = data.shape[0]

        bin_S = max(0, min(bin_S, bin_end))
        bin_L = max(0, min(bin_L, bin_end))
        bin_E = max(0, min(bin_E, bin_end))
        
        safe_S = min(bin_S, valid_mask.shape[0]-1)
        safe_L = min(bin_L, valid_mask.shape[0]-1)
        safe_E = min(bin_E, valid_mask.shape[0]-1)
        if valid_mask[safe_S, safe_S] == 0 or valid_mask[safe_L, safe_L] == 0 or valid_mask[safe_E, safe_E] == 0:
            snip["data"] = np.full((self.rescale_size, self.rescale_size), np.nan)
            return snip

        bins = [0, bin_S, bin_L, bin_E, bin_end]
        
        N = self.rescale_size 
        pad_frac = self.rescale_flank 
        total_len = 2 * pad_frac + 1.0
        
        pix_S = int(N * pad_frac / total_len)
        pix_L = int(N * (pad_frac + 0.5) / total_len)
        pix_E = int(N * (pad_frac + 1.0) / total_len)
        pixs = [0, pix_S, pix_L, pix_E, N]
        
        row_blocks = []
        valid_row_blocks = []
        for i in range(4):
            col_blocks = []
            valid_col_blocks = []
            h_target = pixs[i+1] - pixs[i]
            
            for j in range(4):
                w_target = pixs[j+1] - pixs[j]
                target_shape = (h_target, w_target)
                
                q = data[bins[i]:bins[i+1], bins[j]:bins[j+1]]
                q_valid = valid_mask[bins[i]:bins[i+1], bins[j]:bins[j+1]]
                
                if h_target == 0 or w_target == 0:
                    r = np.zeros(target_shape)
                    r_val = np.zeros(target_shape)
                elif q.shape[0] == 0 or q.shape[1] == 0:
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


def run_3point_pileup():
    bed_file = f"{PROJECT_ROOT}/ADA/K562-loading-sites-within-convergent-loops.bed"
    MCOOL_DIR = f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool"
    
    mcool_files = [
        # os.path.join(MCOOL_DIR, "K562_MC_WT_1kb.mcool::/resolutions/2000"),
        # os.path.join(MCOOL_DIR, "K562_CP_WT_CTCF_merge_1kb.mcool::/resolutions/2000"),
        # os.path.join(MCOOL_DIR, "K562_CP_WT_RAD21_merge_1kb.mcool::/resolutions/2000"),
        # os.path.join(MCOOL_DIR, "K562_CP_RAD21_6h_9h_merge_RAD21_1kb.mcool::/resolutions/2000"),
        # os.path.join(MCOOL_DIR, "K562_MC_CTCF_6h_9h_merge_007_026_1kb.mcool::/resolutions/2000"),
        # os.path.join(MCOOL_DIR, "K562_MC_RAD21_6h_9h_merge_024_031_1kb.mcool::/resolutions/2000"),
        # f"{PROJECT_ROOT}/GM12878/GM12878_in-situ-hic_4DNFI1UEG1HD.mcool::/resolutions/2000"
        # f"{DATA_ROOT}/mly/ChIA-PET/mcool_all/CTCF_0h_CTCF_045_deep.mcool::/resolutions/2000",
        # f"{DATA_ROOT}/mly/ChIA-PET/mcool_all/CTCF_6h_CTCF_044_deep.mcool::/resolutions/2000",
        # f"{DATA_ROOT}/mly/ChIA-PET/mcool_all/WT_CTCF_019_2kb.cool::/resolutions/2000",
        # f"{DATA_ROOT}/mly/ChIA-PET/mcool_all/CTCF_9h_CTCF_047_deep.mcool::/resolutions/2000"
        # f"{PROJECT_ROOT}/HCT116/public_hic/ENCFF135MUT.cool",
        # f"{PROJECT_ROOT}/HCT116/public_hic/ENCFF259YUS.cool",
        # f"{PROJECT_ROOT}/HCT116/public_hic/ENCFF317OIA.cool",
        # f"{PROJECT_ROOT}/HCT116/public_hic/ENCFF528XGK.cool"
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.allele.M.G1.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.allele.M.G2M.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.allele.M.S.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.allele.P.G1.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.allele.P.G2M.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.allele.P.S.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.G1.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.G2M.e500.juice.gz.mcool::/resolutions/2000",
        f"{PROJECT_ROOT}/cell_cycle/k562.PET.S.e500.juice.gz.mcool::/resolutions/2000"
    ]
    
    # ======== Custom parameters ========
    procs = 16
    flank_ratio = 0.5
    min_dist_threshold = 20000 
    cmap_theme = 'black' 
    
    # Set explicit range to force vmin/vmax (e.g., 0.0005), or None to auto percentile
    # custom_vmin = 0.0005 
    # custom_vmax = 0.01

    custom_vmin = None
    custom_vmax = None
    # ==============================

    if cmap_theme == 'black':
        colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
        custom_cmap = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
    else:
        custom_cmap = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
    custom_cmap.set_bad('#FFFFFF')

    for cool_file in mcool_files:
        base_name = os.path.basename(cool_file.split('::')[0]).replace('.mcool', '')
        out_png = f"K562_cell_cycle/coolpuppy_rescaled_3point_{base_name}_2K.pdf"
        out_txt = f"K562_cell_cycle/coolpuppy_rescaled_3point_{base_name}_2K.txt"
        
        print(f"\n{'='*50}")
        print(f"Processing: {base_name}")
        
        clr = cooler.Cooler(cool_file)
        has_chr_prefix = any(chrom.startswith('chr') for chrom in clr.chromnames)
        
        df_3p = prep_3point_df(bed_file, min_la_dist=min_dist_threshold, add_chr=has_chr_prefix)
        retained_count = len(df_3p)
        print(f"Extracted {retained_count} valid S-L-E interactions after distance filtering.")

        if retained_count == 0:
            print(f"Skipping {base_name}: No valid regions left after filtering.")
            continue
        
        CC = CoordCreator(
            features=df_3p,
            resolution=clr.binsize,
            features_format="bed",
            local=True,
            nshifts=0,       
            rescale_flank=flank_ratio, 
        )

        PU = PileUpper3P(
            clr=clr,
            CC=CC,
            rescale=True,
            clr_weight_name=None,
            rescale_size=99,
            ignore_diags=0,         
            nproc=procs,
        )
        
        pups_df = PU.pileupsWithControl(nproc=procs)
        matrix = pups_df['data'].values[0]

        np.savetxt(out_txt, matrix, delimiter='\t', fmt='%.6e')
        
        fig, ax = plt.subplots(figsize=(6, 5))
        
        valid_vals = matrix[~np.isnan(matrix) & (matrix > 0)]
        
        if custom_vmax is not None:
            vmax = custom_vmax
        else:
            vmax = np.percentile(valid_vals, 95) if len(valid_vals) > 0 else 1.0
            
        if custom_vmin is not None:
            vmin = custom_vmin
        else:
            vmin = np.percentile(valid_vals, 5) if len(valid_vals) > 0 else 0.0

        im = ax.imshow(matrix, cmap=custom_cmap, norm=Normalize(vmin=vmin, vmax=vmax))
        
        N = matrix.shape[0]
        total_len = 2 * flank_ratio + 1.0
        pix_S = int(N * flank_ratio / total_len)         
        pix_L = int(N * (flank_ratio + 0.5) / total_len) 
        pix_E = int(N * (flank_ratio + 1.0) / total_len) 
        
        ax.set_xticks([pix_S, pix_L, pix_E])
        ax.set_xticklabels(['S', 'L', 'E'])
        ax.set_yticks([pix_S, pix_L, pix_E])
        ax.set_yticklabels(['S', 'L', 'E'])
        
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # Add plot title showing source and retained count
        ax.set_title(f"Source: {base_name}\nRetained Regions: 160", fontsize=10)
        
        plt.colorbar(im, label='Balanced Contact Frequency')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Finished {base_name}.")

if __name__ == "__main__":
    run_3point_pileup()