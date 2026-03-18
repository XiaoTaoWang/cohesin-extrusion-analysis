import joblib
import numpy as np
import matplotlib.pyplot as plt
import sys
import cooler
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import matplotlib.cm as cm

# --- 1. 参数与路径配置 ---
# 假设用法: python diff_plot.py ./mut_dir ./wild_dir title
mut_dir = sys.argv[1]
wild_dir = sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else "Mutant - Wild"

paramfile = "./data/all_ctcf_loading/cr139.pkl" 
loop_id, CHROM, START, END, RES, CTCF_left_positions, CTCF_right_positions, site_types = joblib.load(paramfile)

start_bin = (37302000 - START) // RES
end_bin = (38102000 - START) // RES

idx = loop_id
mut_cool = f'{mut_dir}/cool_outputs/{idx}.2000.cool'
wild_cool = f'{wild_dir}/cool_outputs/{idx}.2000.cool'

# --- 2. 读取并归一化数据 ---
def get_matrix(path):
    c = cooler.Cooler(path)
    mat = c.matrix(balance=False, sparse=False).fetch(c.chromnames[0]).astype(float)
    # 归一化：基于对角线中位数，确保两张图在同一量级
    return mat

mat_mut = get_matrix(mut_cool)
mat_wild = get_matrix(wild_cool)

# 计算差异矩阵
mat_diff = mat_mut - mat_wild
#mat_diff = mat_diff[:400, :400]
mat_diff = mat_diff[start_bin:end_bin, start_bin:end_bin]    
#mat /= np.median(np.diag(mat, 2))

# --- 3. 绘图设置 ---
fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 2], 'hspace': 0.1}, figsize=(8, 10))
ax_diff = axes[0]
ax_track = axes[1]

# 设置差异图的 Colormap: 类似图中右侧的蓝-白-红
# 'bwr' (blue-white-red) 或 'seismic' 是最常用的
CMAP_DIFF = cm.bwr 

# 动态计算 vmin/vmax，确保 0 对应白色
v_abs = np.nanmax(np.abs(mat_diff)) * 0.7 # 取最大绝对值的70%以增加对比度
v_abs = 300
norm = TwoSlopeNorm(vcenter=0, vmin=-v_abs, vmax=v_abs)

im = ax_diff.imshow(mat_diff, norm=norm, cmap=CMAP_DIFF)
#ax_diff.set_title(f"Difference: {title}")

ax_diff.set_xticks([])
ax_diff.set_yticks([])

#for spine in ax_diff.spines.values():
#    spine.set_linewidth(0)

#ax_diff.axis('off')

# 添加颜色条
cbar = plt.colorbar(im, ax=ax_diff, fraction=0.046, pad=0.04, shrink=0.6)
cbar.set_label('$\Delta$ Interaction Strength')
'''
# --- 4. 绘制下方轨道 (与原脚本一致) ---
ctcf_left_plot = sorted(list(set(CTCF_left_positions // 10))) 
ctcf_right_plot = sorted(list(set(CTCF_right_positions // 10)))
loading_plot = sorted(list(set(np.nonzero(site_types)[0] // 10)))

ax_track.scatter(ctcf_left_plot, [1]*len(ctcf_left_plot), s=15, marker='>', color='green', label='CTCF R')
ax_track.scatter(ctcf_right_plot, [0]*len(ctcf_right_plot), s=15, marker='<', color='blue', label='CTCF L')
ax_track.scatter(loading_plot, [2]*len(loading_plot), s=15, c='red', marker='|', label='Loading')

ax_track.set_yticks([0, 1, 2])
ax_track.set_yticklabels(['CTCF R', 'CTCF L', 'Loading'], fontsize=8)
#ax_track.set_xlim(0, 400) #mat_diff.shape[1])
#ax_track.set_xlim(0, mat_diff.shape[1])
ax_track.set_xlim(start_bin, end_bin)
ax_track.set_ylim(-1, 3)
ax_track.grid(axis='x', linestyle='--', alpha=0.3)

# 调整轨道位置
pos = ax_diff.get_position()
ax_track.set_position([pos.x0, pos.y0 - 0.08, pos.width, 0.06])
'''

# 保存
save_name = f'{mut_dir}-{wild_dir}.diff.zoom2.svg'
plt.savefig(save_name, bbox_inches='tight')
print(f"Difference plot saved as {save_name}")
