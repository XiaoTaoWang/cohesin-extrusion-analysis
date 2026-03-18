import joblib
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.colors import TwoSlopeNorm
import matplotlib.cm as cm

# --- 1. 参数与路径配置 ---
# 假设用法: python agg_diff.py ./mut_dir ./wild_dir title
if len(sys.argv) < 3:
    print("Usage: python agg_diff.py <mut_dir> <wild_dir> [title]")
    sys.exit(1)

mut_dir = sys.argv[1]
wild_dir = sys.argv[2]
title = sys.argv[3] # if len(sys.argv) > 3 else "Mutant - Wild"
thr = sys.argv[4]

# 加载聚合后的数据
mat_mut = joblib.load(f'{mut_dir}/sim_rescale_agg.{thr}.pkl')
mat_wild = joblib.load(f'{wild_dir}/sim_rescale_agg.{thr}.pkl')

# 计算差异矩阵
mat_diff = mat_mut - mat_wild

# --- 3. 绘图设置 ---
# 修改为单图布局，删除 height_ratios
fig, ax_diff = plt.subplots(figsize=(6, 5))

# 设置差异图的 Colormap: 蓝-白-红 (bwr)
#CMAP_DIFF = cm.bwr 
CMAP_DIFF = 'RdYlBu_r'

# 动态计算 vmin/vmax，确保 0 对应白色
v_abs = 500
norm = TwoSlopeNorm(vcenter=0, vmin=-v_abs, vmax=v_abs)

# 绘制热图
im = ax_diff.imshow(mat_diff, norm=norm, cmap=CMAP_DIFF, interpolation='none')
ax_diff.set_title(f"Difference: {title}")

# 添加颜色条，保留原始比例
cbar = plt.colorbar(im, ax=ax_diff, fraction=0.046, pad=0.04, shrink=0.6)
cbar.set_label('$\Delta$ Interaction Strength')

# 装饰：去除刻度
ax_diff.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

# 保存文件
save_name = f'{mut_dir}-{wild_dir}.diff.agg.{thr}.svg'
plt.savefig(save_name, bbox_inches='tight', dpi=300)

print(f"Aggregated difference plot saved as {save_name}")
