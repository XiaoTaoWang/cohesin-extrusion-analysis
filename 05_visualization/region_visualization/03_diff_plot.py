
# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
import os
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import os
import cooler, matplotlib, joblib
matplotlib.use('Agg')
from multiprocess import Pool
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from skimage import transform

def parse_regions(fil, remove_chr=True):
    anchors = {}
    loadings = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e, label = line.rstrip().split()
            if remove_chr:
                c = c.lstrip('chr')
            s, e = int(s), int(e)
            p = (s + e) // 2
            key = label.split('_')[0]
            if not label.split('_')[1].startswith('M'):
                if not key in anchors:
                    anchors[key] = []
                anchors[key].append((c, p))
            else:
                if not key in loadings:
                    loadings[key] = []
                loadings[key].append((c, p))
    
    triples = {} # by chrom
    for key in anchors:
        S, E = anchors[key]
        c = S[0]
        for M in loadings[key]:
            if not c in triples:
                triples[c] = []
            triples[c].append((c, S[1], M[1], E[1]))
    
    return triples

def collect_submatrices(uri, triples, balance=False, minsize=10):
    lib = cooler.Cooler(uri)
    res = lib.binsize
    fixed_length = 80
    fixed_shape = (fixed_length//2, fixed_length//2)

    arrs = []
    filtered_triples = []
    for c in triples:
        if not c in lib.chromnames:
            continue

        M = lib.matrix(balance=balance, sparse=True).fetch(c).tocsr()
        for _, s, m, e in triples[c]:
            s = s // res
            m = m // res
            e = e // res
            el = (m - s) // 2
            er = (e - m) // 2
            
            if (m - s < minsize) or (e - m < minsize):
                continue

            if (s - el >= 0) and (e + er + 1 <= M.shape[0]):
                sub1 = M[s-el:m, s-el:m].toarray().astype(float)
                sub1[np.isnan(sub1)] = 0
                sub2 = M[s-el:m, m:e+er+1].toarray().astype(float)
                sub2[np.isnan(sub2)] = 0
                sub3 = M[m:e+er+1, s-el:m].toarray().astype(float)
                sub3[np.isnan(sub3)] = 0
                sub4 = M[m:e+er+1, m:e+er+1].toarray().astype(float)
                sub4[np.isnan(sub4)] = 0
                
                resized_sub1 = transform.resize(sub1, output_shape=fixed_shape)
                resized_sub2 = transform.resize(sub2, output_shape=fixed_shape)
                resized_sub3 = transform.resize(sub3, output_shape=fixed_shape)
                resized_sub4 = transform.resize(sub4, output_shape=fixed_shape)
                
                big = np.zeros((fixed_length, fixed_length))
                big[:fixed_shape[0], :fixed_shape[0]] = resized_sub1
                big[:fixed_shape[0], fixed_shape[0]:] = resized_sub2
                big[fixed_shape[0]:, :fixed_shape[0]] = resized_sub3
                big[fixed_shape[0]:, fixed_shape[0]:] = resized_sub4
                
                mean_val = big.mean()
                if mean_val > 0:
                    big = big / mean_val

                arrs.append(big)
                filtered_triples.append((c, s, m, e))
    
    return arrs, filtered_triples

def pileup(arrs):
    if len(arrs) == 0:
        raise ValueError("No valid submatrices found. Check input parameters or regions.")
    pool = np.r_[arrs]
    avg = pool.mean(axis=0)
    return avg

# ===============================
# PARAMETERS & CONFIGURATION
# ===============================
bed_file = f"{PROJECT_ROOT}/ADA/K562-loading-sites-within-convergent-loops.bed"
out_dir = f"{PROJECT_ROOT}/aggregation_plot/rescale_3point/Xiaotao_script"  # Renamed directory to reflect linear differences
procs = 10
flank_ratio = 0.5
min_dist = 20000 
resolution = 2000

MCOOL_DIR = f"{DATA_ROOT}/data_enhancement_1kb_cool_LY/mcool"

mc_wt = os.path.join(MCOOL_DIR, "K562_MC_WT_1kb.mcool::/resolutions/2000")
mc_ctcf = os.path.join(MCOOL_DIR, "K562_MC_CTCF_6h_9h_merge_007_026_1kb.mcool::/resolutions/2000")
mc_rad21 = os.path.join(MCOOL_DIR, "K562_MC_RAD21_6h_9h_merge_024_031_1kb.mcool::/resolutions/2000")

cp_wt_rad21 = os.path.join(MCOOL_DIR, "K562_CP_WT_RAD21_merge_1kb.mcool::/resolutions/2000")
cp_rad21 = os.path.join(MCOOL_DIR, "K562_CP_RAD21_6h_9h_merge_RAD21_1kb.mcool::/resolutions/2000")

comparisons = [
    ("MC_diff_RAD21-AID_vs_CTCF-AID", mc_rad21, "MC_RAD21", mc_ctcf, "MC_CTCF"),
    ("MC_diff_RAD21-AID_vs_WT", mc_rad21, "MC_RAD21", mc_wt, "MC_WT"),
    ("MC_diff_CTCF-AID_vs_WT", mc_ctcf, "MC_CTCF", mc_wt, "MC_WT"),
    ("CP_diff_RAD21-AID_vs_WT", cp_rad21, "CP_RAD21", cp_wt_rad21, "CP_WT_RAD21")
]

remove_chr = False
balance = False

# Create output directory if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

# Parse regions once
print(f"Parsing regions from {bed_file}...")
triples = parse_regions(bed_file, remove_chr=remove_chr)
min_bins = min_dist // resolution

# Loop through each comparison
for comp_name, uri_trt, name_trt, uri_ctl, name_ctl in comparisons:
    print(f"\nProcessing {comp_name}...")
    
    # 1. Collect Treatment
    print(f"  Fetching treatment: {name_trt}")
    arrs_trt, _ = collect_submatrices(uri_trt, triples, balance=balance, minsize=min_bins)
    avg_trt = pileup(arrs_trt)
    
    # 2. Collect Control
    print(f"  Fetching control: {name_ctl}")
    arrs_ctl, _ = collect_submatrices(uri_ctl, triples, balance=balance, minsize=min_bins)
    avg_ctl = pileup(arrs_ctl)
    
    # 3. Calculate Linear Difference (Treatment - Control)
    print("  Calculating linear difference...")
    diff_map = avg_trt - avg_ctl
    
    # Save Pickled array
    outpkl = os.path.join(out_dir, f"{comp_name}.pkl")
    joblib.dump(diff_map, outpkl)
    
    # 4. Plotting
    outfig = os.path.join(out_dir, f"{comp_name}.svg")
    
    # Adjust vmax_val for linear differences. You might need to tune this 
    # depending on the absolute magnitude of your interaction differences.
    vmax_val = 0.2 
    vmin_val = -vmax_val
    
    size = (2.2, 2)
    fig = plt.figure(figsize=size)
    width = 0.7; Left = 0.1
    HB = 0.1; HH = width * size[0] / size[1]
    ax = fig.add_axes([Left, HB, width, HH])
    
    # Using 'bwr' (Blue-White-Red) divergent colormap
    sc = ax.imshow(diff_map, cmap='bwr', vmin=vmin_val, vmax=vmax_val)
    
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in ['right', 'top', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(0.9)
        
    cax = fig.add_axes([Left+width+0.04, 0.72, 0.03, 0.15])
    fig.colorbar(sc, cax=cax, ticks=[vmin_val, 0, vmax_val], format='%.2g')
    
    plt.savefig(outfig, dpi=500, bbox_inches='tight')
    plt.close()
    print(f"  Saved {outfig}")

print("\nDone processing all comparisons.")