import os
import sys
import glob
import joblib
import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from skimage import transform

data_dir = sys.argv[1]     
avg = joblib.load(f'{data_dir}/sim_rescale_agg.pkl')
# 绘图
plt.figure(figsize=(2.2, 2))
ax = plt.axes([0.15, 0.15, 0.7, 0.7])

# 模拟数据信号通常很强，vmin/vmax 需根据实际情况微调
#im = ax.imshow(avg, cmap='RdBu_r', norm=LogNorm(np.percentile(avg, 5), np.percentile(avg, 85))) #vmax=np.percentile(avg, 98)) #, norm=LogNorm(vmin=0.05, vmax=1.5))
#im = ax.imshow(avg, cmap='RdBu_r', norm=LogNorm(np.percentile(avg, 5), np.percentile(avg, 85))) #vmax=np.percentile(avg, 98)) #, norm=LogNorm(vmin=0.05, vmax=1.5))
#im = ax.imshow(avg, cmap=CMAP, norm=LogNorm(np.percentile(avg, 3), np.percentile(avg, 85))) #vmax=np.percentile(avg, 98)) #, norm=LogNorm(vmin=0.05, vmax=1.5))
#im = ax.imshow(avg, cmap=CMAP, norm=LogNorm(vmin=0.05, vmax=1.5))
vmax = np.percentile(avg, 90)
vmin = np.percentile(avg, 3)
print(vmax, vmin)
#vmax, vmin = 2000, 50
vmax, vmin = 3000, 50
#vmax, vmin = 0.205881383896166, 0.00446084646430251
# loading0 4964.141669663833 42.320541970676594
CMAP = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#F70000'])
#fall_colors = ["#FFFFFF", "#FFFAD0", "#FFD700", "#FF8C00", "#FF0000", "#000000"]
#CMAP = LinearSegmentedColormap.from_list("fall", fall_colors)
colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
CMAP = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
#im = ax.imshow(avg, cmap=CMAP, vmin=vmin, vmax=vmax)
im = ax.imshow(avg, cmap=CMAP, norm=LogNorm(vmin=vmin, vmax=vmax))

cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.3)

# 装饰
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# 画辅助线标识 S, M, E 的位置 (20, 40, 60 bin 处)
# ax.axvline(20, color='gray', lw=0.5, ls='--')
# ax.axvline(60, color='gray', lw=0.5, ls='--')
plt.title(f'{data_dir}')
plt.savefig(f'{data_dir}.agg.svg', dpi=300, bbox_inches='tight')
