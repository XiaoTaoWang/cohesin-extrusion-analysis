
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
import pickle
import cooler
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.ndimage import zoom
from coolpuppy.coolpup import CoordCreator, PileUpper
from matplotlib.colors import Normalize, LinearSegmentedColormap

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

        bin_S = max(0, min(int((s_bp - exp_start) // res), data.shape[0]))
        bin_L = max(0, min(int((l_bp - exp_start) // res), data.shape[0]))
        bin_E = max(0, min(int((e_bp - exp_start) // res), data.shape[0]))
        bin_end = data.shape[0]

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
        
        pixs = [0, 
                int(N * pad_frac / total_len), 
                int(N * (pad_frac + 0.5) / total_len), 
                int(N * (pad_frac + 1.0) / total_len), 
                N]
        
        row_blocks, valid_row_blocks = [], []
        for i in range(4):
            col_blocks, valid_col_blocks = [], []
            h_target = pixs[i+1] - pixs[i]
            
            for j in range(4):
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

def get_pileup_matrix(cool_path, bed_file, min_dist_threshold, flank_ratio, procs):
    if not os.path.exists(cool_path.split('::')[0]):
        return None
        
    clr = cooler.Cooler(cool_path)
    has_chr_prefix = any(chrom.startswith('chr') for chrom in clr.chromnames)
    
    df_3p = prep_3point_df(bed_file, min_la_dist=min_dist_threshold, add_chr=has_chr_prefix)
    if len(df_3p) == 0: return None

    CC = CoordCreator(features=df_3p, resolution=clr.binsize, features_format="bed", local=True, nshifts=0, rescale_flank=flank_ratio)
    PU = PileUpper3P(clr=clr, CC=CC, clr_weight_name=None, rescale=True, rescale_size=99, ignore_diags=0, nproc=procs)
    
    return PU.pileupsWithControl(nproc=procs)['data'].values[0]

def get_or_calc_matrix(cool_path, bed_file, min_dist, flank, procs, cache_dir):
    if not cool_path or not os.path.exists(cool_path.split('::')[0]):
        return None
    safe_name = os.path.basename(cool_path.split('::')[0]).replace('.mcool', '').replace('.cool', '')
    pkl_path = os.path.join(cache_dir, f"matrix_{safe_name}_3point.pkl")
    
    if os.path.exists(pkl_path):
        with open(pkl_path, 'rb') as f:
            return pickle.load(f)
    else:
        mat = get_pileup_matrix(cool_path, bed_file, min_dist, flank, procs)
        if mat is not None:
            with open(pkl_path, 'wb') as f:
                pickle.dump(mat, f)
        return mat

def plot_single_matrix(matrix, out_file, title, flank, cmap, is_log2=False, is_linear=False):
    fig, ax = plt.subplots(figsize=(6, 5))
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad('#FFFFFF')
    
    if is_log2:
        # vmax = np.nanpercentile(np.abs(matrix), 98)
        vmax = 1.5
        if vmax == 0 or np.isnan(vmax): vmax = 1.0
        vmin = -vmax
        label = 'Log2 Fold Change'
    elif is_linear:
        valid_vals = matrix[~np.isnan(matrix)]
        if len(valid_vals) > 0:
            # max_abs = np.percentile(np.abs(valid_vals), 98)
            max_abs = 1.0
            vmin, vmax = -max_abs, max_abs
        else:
            vmin, vmax = -0.001, 0.001
        label = 'Delta Contact Frequency'
    else:
        # vmin = np.nanpercentile(matrix, 5)
        # vmax = np.nanpercentile(matrix, 95)
        vmin = 0.2
        vmax = 5.0
        if vmin == vmax or np.isnan(vmax): vmax = vmin + 1.0
        label = 'Contact frequency'
        
    im = ax.imshow(matrix, cmap=cmap_obj, norm=Normalize(vmin=vmin, vmax=vmax))
    
    N = matrix.shape[0]
    total_len = 2 * flank + 1.0
    pixs = [int(N * flank / total_len), int(N * (flank + 0.5) / total_len), int(N * (flank + 1.0) / total_len)]
    
    ax.set_xticks(pixs)
    ax.set_xticklabels(['S', 'L', 'E'])
    ax.set_yticks(pixs)
    ax.set_yticklabels(['S', 'L', 'E'])
    ax.grid(False)
    for spine in ax.spines.values(): spine.set_visible(False)
        
    ax.set_title(title.replace('_', ' '), fontsize=12)
    plt.colorbar(im, label=label)
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
# bed_file = f"{PROJECT_ROOT}/script_minji/20260213_from_minji/HCT116-cohesin-specific-regions_20230710_filt_25-75percentile_20231030.bed"
 
def run_hct116_analysis():
    bed_file = f"{PROJECT_ROOT}/K562/minji_loading_site_define/K562-cohesin-specific-regions_filt_25-75percentile.bed"
    base_out_dir = f"{PROJECT_ROOT}/aggregation_plot/rescale_3point/K562_linear_norm_new_define_minji"
    log2_dir = os.path.join(base_out_dir, "log2_plot")
    linear_dir = os.path.join(base_out_dir, "linear_minus")
    cache_dir = os.path.join(base_out_dir, "cache")
    
    for d in [log2_dir, linear_dir, cache_dir]: os.makedirs(d, exist_ok=True)
    
    procs, flank_ratio, min_dist = 32, 0.5, 20000 
    MCOOL_DIR = f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool"

    hic_dir = f"{PROJECT_ROOT}/HCT116/public_hic"
    ctcf_0h = os.path.join(MCOOL_DIR, "K562_MC_WT_1kb.mcool::/resolutions/2000")
    ctcf_6h = os.path.join(MCOOL_DIR, "K562_MC_CTCF_6h_9h_merge_007_026_1kb.mcool::/resolutions/2000")
    rad21_0h = os.path.join(MCOOL_DIR, "K562_MC_WT_1kb.mcool::/resolutions/2000")
    rad21_6h = os.path.join(MCOOL_DIR, "K562_MC_RAD21_6h_9h_merge_024_031_1kb.mcool::/resolutions/2000")
    
    # simple comment
    comparisons = [
        ("RAD21_6h_vs_0h", rad21_6h, "RAD21_6h", rad21_0h, "RAD21_0h"),
        ("CTCF_6h_vs_0h", ctcf_6h, "CTCF_6h", ctcf_0h, "CTCF_0h")
    ]
    
    for comp_name, path_A, name_A, path_B, name_B in comparisons:
        mat_A = get_or_calc_matrix(path_A, bed_file, min_dist, flank_ratio, procs, cache_dir)
        mat_B = get_or_calc_matrix(path_B, bed_file, min_dist, flank_ratio, procs, cache_dir)
        
        if mat_A is None or mat_B is None: continue
        
        colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
        custom_cmap = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
        plot_single_matrix(mat_A, os.path.join(log2_dir, f"{name_A}_raw.pdf"), name_A, flank_ratio, custom_cmap)
        plot_single_matrix(mat_B, os.path.join(log2_dir, f"{name_B}_raw.pdf"), name_B, flank_ratio, custom_cmap)
        
        # simple comment
        mat_A_norm = mat_A / np.nanmean(mat_A)
        mat_B_norm = mat_B / np.nanmean(mat_B)
        pseudo = 1e-5
        
        diff_log2 = np.log2((mat_A_norm + pseudo) / (mat_B_norm + pseudo))
        diff_log2 = np.nan_to_num(diff_log2, nan=0.0, posinf=0.0, neginf=0.0)
        np.savetxt(os.path.join(log2_dir, f"log2_{comp_name}.txt"), diff_log2, delimiter='\t', fmt='%.6e')
        plot_single_matrix(diff_log2, os.path.join(log2_dir, f"log2_{comp_name}.pdf"), f"Log2 {comp_name}", flank_ratio, 'RdYlBu_r', is_log2=True)
        
        # simple comment
        diag_mean_A = np.nanmean(np.diag(mat_A))
        diag_mean_B = np.nanmean(np.diag(mat_B))
        scale_factor = diag_mean_B / diag_mean_A  # simple comment
        mat_A_scaled = mat_A * scale_factor
        
        diff_linear = mat_A_scaled - mat_B
        np.savetxt(os.path.join(linear_dir, f"linear_{comp_name}.txt"), diff_linear, delimiter='\t', fmt='%.6e')
        plot_single_matrix(diff_linear, os.path.join(linear_dir, f"linear_{comp_name}.pdf"), f"Linear {comp_name} (Scaled)", flank_ratio, 'RdYlBu_r', is_linear=True)

    # simple comment
    mat_rad21_6h = get_or_calc_matrix(rad21_6h, bed_file, min_dist, flank_ratio, procs, cache_dir)
    mat_rad21_0h = get_or_calc_matrix(rad21_0h, bed_file, min_dist, flank_ratio, procs, cache_dir)
    mat_ctcf_6h = get_or_calc_matrix(ctcf_6h, bed_file, min_dist, flank_ratio, procs, cache_dir)
    mat_ctcf_0h = get_or_calc_matrix(ctcf_0h, bed_file, min_dist, flank_ratio, procs, cache_dir)
    
    if all(m is not None for m in [mat_rad21_6h, mat_rad21_0h, mat_ctcf_6h, mat_ctcf_0h]):
        comp_name_double = "DoubleDiff_RAD21_vs_CTCF_6h_over_0h"
        pseudo = 1e-5
        
        # simple comment
        rad21_norm_ratio = (mat_rad21_6h / np.nanmean(mat_rad21_6h) + pseudo) / (mat_rad21_0h / np.nanmean(mat_rad21_0h) + pseudo)
        ctcf_norm_ratio = (mat_ctcf_6h / np.nanmean(mat_ctcf_6h) + pseudo) / (mat_ctcf_0h / np.nanmean(mat_ctcf_0h) + pseudo)
        diff_log2_double = np.log2(rad21_norm_ratio) - np.log2(ctcf_norm_ratio)
        diff_log2_double = np.nan_to_num(diff_log2_double, nan=0.0, posinf=0.0, neginf=0.0)
        
        np.savetxt(os.path.join(log2_dir, f"log2_{comp_name_double}.txt"), diff_log2_double, delimiter='\t', fmt='%.6e')
        plot_single_matrix(diff_log2_double, os.path.join(log2_dir, f"log2_{comp_name_double}.pdf"), "Log2 Delta (RAD21 - CTCF)", flank_ratio, 'RdYlBu_r', is_log2=True)
        
        # simple comment
        # simple comment
        scale_rad21 = np.nanmean(np.diag(mat_rad21_0h)) / np.nanmean(np.diag(mat_rad21_6h))
        mat_rad21_6h_scaled = mat_rad21_6h * scale_rad21
        diff_rad21 = mat_rad21_6h_scaled - mat_rad21_0h
        
        # simple comment
        scale_ctcf = np.nanmean(np.diag(mat_ctcf_0h)) / np.nanmean(np.diag(mat_ctcf_6h))
        mat_ctcf_6h_scaled = mat_ctcf_6h * scale_ctcf
        diff_ctcf = mat_ctcf_6h_scaled - mat_ctcf_0h
        
        # simple comment
        diff_linear_double = diff_rad21 - diff_ctcf
        
        np.savetxt(os.path.join(linear_dir, f"linear_{comp_name_double}.txt"), diff_linear_double, delimiter='\t', fmt='%.6e')
        plot_single_matrix(diff_linear_double, os.path.join(linear_dir, f"linear_{comp_name_double}.pdf"), "Linear Delta (RAD21 - CTCF)", flank_ratio, 'RdYlBu_r', is_linear=True)

    if all(m is not None for m in [mat_rad21_6h, mat_rad21_0h, mat_ctcf_6h, mat_ctcf_0h]):
        comp_name_double = "DoubleDiff_CTCF_vs_RAD21_6h_over_0h"
        pseudo = 1e-5
        
        # simple comment
        rad21_norm_ratio = (mat_rad21_6h / np.nanmean(mat_rad21_6h) + pseudo) / (mat_rad21_0h / np.nanmean(mat_rad21_0h) + pseudo)
        ctcf_norm_ratio = (mat_ctcf_6h / np.nanmean(mat_ctcf_6h) + pseudo) / (mat_ctcf_0h / np.nanmean(mat_ctcf_0h) + pseudo)
        diff_log2_double = np.log2(ctcf_norm_ratio) - np.log2(rad21_norm_ratio)
        diff_log2_double = np.nan_to_num(diff_log2_double, nan=0.0, posinf=0.0, neginf=0.0)
        
        np.savetxt(os.path.join(log2_dir, f"log2_{comp_name_double}.txt"), diff_log2_double, delimiter='\t', fmt='%.6e')
        plot_single_matrix(diff_log2_double, os.path.join(log2_dir, f"log2_{comp_name_double}.pdf"), "Log2 Delta (RAD21 - CTCF)", flank_ratio, 'RdYlBu_r', is_log2=True)
        
        # simple comment
        # simple comment
        scale_rad21 = np.nanmean(np.diag(mat_rad21_0h)) / np.nanmean(np.diag(mat_rad21_6h))
        mat_rad21_6h_scaled = mat_rad21_6h * scale_rad21
        diff_rad21 = mat_rad21_6h_scaled - mat_rad21_0h
        
        # simple comment
        scale_ctcf = np.nanmean(np.diag(mat_ctcf_0h)) / np.nanmean(np.diag(mat_ctcf_6h))
        mat_ctcf_6h_scaled = mat_ctcf_6h * scale_ctcf
        diff_ctcf = mat_ctcf_6h_scaled - mat_ctcf_0h
        
        # simple comment
        diff_linear_double = diff_ctcf - diff_rad21
        
        np.savetxt(os.path.join(linear_dir, f"linear_{comp_name_double}.txt"), diff_linear_double, delimiter='\t', fmt='%.6e')
        plot_single_matrix(diff_linear_double, os.path.join(linear_dir, f"linear_{comp_name_double}.pdf"), "Linear Delta (CTCF - RAD21)", flank_ratio, 'RdYlBu_r', is_linear=True)

if __name__ == "__main__":
    run_hct116_analysis()