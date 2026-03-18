import joblib
import numpy as np
import matplotlib.pyplot as plt
import sys
import cooler
import glob
import os
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from plot_loops import plot_loops

# --- 1. 配置与参数加载 ---
paramfile = "./data/all_ctcf_loading/cr139.pkl" 
loop_id, CHROM, START, END, RES, CTCF_left_positions, CTCF_right_positions, site_types = joblib.load(paramfile)

data_dir = sys.argv[1]
# 获取该 loop 所有的演化 cool 文件，并按时间排序
#cool_files = sorted(glob.glob(f'{data_dir}/cool_with_time/{loop_id}_upto_t*.cool'), 
cool_files = sorted(glob.glob(f'{data_dir}/{loop_id}_upto_t*.cool'), 
                    key=lambda x: int(x.split('_t')[-1].split('.')[0]))

if not cool_files:
    print(f"No evolution cool files found in {data_dir}")
    sys.exit()

# 设置绘图坐标
ctcf_left_plot = sorted(list(set((CTCF_left_positions // 10))))
ctcf_right_plot = sorted(list(set((CTCF_right_positions // 10))))
loading_plot = sorted(list(set(np.nonzero(site_types)[0] // 10))) # 根据你原脚本的 site_types 步长调整

# --- 2. 绘图布局 ---
n_timepoints = len(cool_files)
fig, axes = plt.subplots(2, n_timepoints, 
                         gridspec_kw={'height_ratios': [5, 2], 'wspace': 0.1})
#figsize=(4 * n_timepoints, 6), 

vmax, vmin = 0.0602, 0.004
CMAP = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#F70000'])

# --- 3. 循环绘制每个时间点 ---
for i, cool_path in enumerate(cool_files):
    # 提取时间标签
    time_label = os.path.basename(cool_path).split('_t')[-1].split('.')[0]
    
    # 获取热图和轨道图的 axis
    ax_sim = axes[0, i] if n_timepoints > 1 else axes[0]
    ax_track = axes[1, i] if n_timepoints > 1 else axes[1]
    
    # 读取数据
    c = cooler.Cooler(cool_path)
    mat = c.matrix(balance=False, sparse=False).fetch(c.chromnames[0]).astype(float)
    mat /= np.median(np.diag(mat, 2)) # 归一化
    
    # 绘制热图
    im = ax_sim.imshow(mat, norm=LogNorm(vmin, vmax), cmap=CMAP)
    loop_file = '/home/wangxiaotao/jianglinghan/project/yijun/data/K562.CTCF_loops.annotated.bedpe'
    #plot_loops(ax_sim, loop_file, CHROM, START, END, res=2000, color='cyan', marker_size=3, linewidths=0.2)
    ax_sim.set_title(f"Time: {time_label}", fontsize=10)
    
    # 移除热图刻度
    ax_sim.set_xticks([])
    ax_sim.set_yticks([])

    # 绘制下方轨道图
    ax_track.scatter(ctcf_left_plot, [1]*len(ctcf_left_plot), s=1, marker='>', color='green')
    ax_track.scatter(ctcf_right_plot, [0]*len(ctcf_right_plot), s=1, marker='<', color='blue')
    ax_track.scatter(loading_plot, [2]*len(loading_plot), s=1, c='red', marker='|')
    
    ax_track.set_xlim(0, mat.shape[1])
    ax_track.set_ylim(-1, 3)
    ax_track.set_xticks([])
    if i == 0:
        ax_track.set_yticks([0, 1, 2])
        ax_track.set_yticklabels(['R', 'L', 'M'], fontsize=3)
    else:
        ax_track.set_yticks([])
    pos = ax_sim.get_position()
    ax_track.set_position([pos.x0, pos.y0 - 0.06, pos.width, 0.04])

# 在最右侧添加颜色条
#cbar_ax = fig.add_axes([0.96, 0.5, 0.015, 0.1])
#fig.colorbar(im, cax=cbar_ax, label='Intensity')

# 保存
output_name = f'{data_dir}.{loop_id}_evolution.svg'
plt.savefig(output_name, bbox_inches='tight', dpi=300)
print(f"Evolution plot saved as {output_name}")
