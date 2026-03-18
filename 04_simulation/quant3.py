import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import joblib
import sys

indir, outdir = sys.argv[1:3]

CMAP = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#F70000'])

def get_window_mean(matrix, r, c, size=0):
    """提取以 (r, c) 为中心的 (2*size+1)x(2*size+1) 窗口均值"""
    r_start, r_end = max(0, r - size), min(matrix.shape[0], r + size + 1)
    c_start, c_end = max(0, c - size), min(matrix.shape[1], c + size + 1)
    window = matrix[r_start:r_end, c_start:c_end]
    return np.median(window) if window.size > 0 else 0

def distance_norm(M):
    """
    输入矩阵 M (1000x1000)，输出 O/E 矩阵 M_oe
    """
    n = M.shape[0]
    M = M.astype(float)
    expected = np.zeros_like(M)

    for s in range(n):
        # 提取第 s 条对角线
        diag = np.diag(M, k=s)
        if diag.size > 0:
            mean_val = np.nanmean(diag)
            if mean_val > 0:
                # 填充对角线
                idx = np.arange(n - s)
                expected[idx, idx + s] = mean_val
                expected[idx + s, idx] = mean_val

    # 计算 Observed / Expected
    # 避免除以 0
    M_oe = np.divide(M, expected, out=np.ones_like(M), where=expected != 0)
    return M_oe

def plot_and_quantify(avg_matrix, output_prefix):
    if 'loading' in output_prefix:
        title = output_prefix.replace('loading', 'Cohesin')
    elif 'wild' in output_prefix:
        title = 'Wild'
    elif 'ctcf' in output_prefix:
        title = output_prefix.replace('ctcf', 'CTCF')
    else:
        title = output_prefix
    title = title.replace('_', ' ')
    title = title.replace('pct', '%')

    # 标准化坐标 (基于 80x80 矩阵)
    S_bin, M_bin, E_bin = 20, 39, 59
    s = 1
    
    # 1. 计算 C-C 强度 (左锚点与右锚点的交互)
    val_cc = get_window_mean(avg_matrix, S_bin, E_bin, size=s)
    
    # 2. 计算 L-C 强度 (加载位点与两个锚点的交互均值)
    val_ls = get_window_mean(avg_matrix, S_bin, M_bin, size=s)
    val_le = get_window_mean(avg_matrix, M_bin, E_bin, size=s)
    val_lc = (val_ls + val_le) / 2.0
    
    # 3. 计算相对强度 (Ratio)
    ratio = val_cc / val_lc if val_lc > 0 else 0
    
    # --- 绘图 ---
    vmin, vmax = 0.02, 0.8
    fig, ax = plt.subplots(figsize=(2.5, 2.2))
    
    # 绘制热图
    #im = ax.imshow(avg_matrix, cmap='RdBu_r', norm=LogNorm(vmin=vmin, vmax=vmax))
    #im = ax.imshow(avg, cmap='RdBu_r', norm=LogNorm(np.percentile(avg, 5), np.percentile(avg, 100)))
    vmax = np.percentile(avg, 99.5)
    vmin = np.percentile(avg, 25)
    #vmax, vmin = 3000, 50
    #vmax, vmin = 2, 0.5
    #vmax, vmin = 1, 0.1
    vmax, vmin = 3, 0.5
    print(vmax, vmin)
    #im = ax.imshow(avg, cmap=CMAP, norm=LogNorm(np.percentile(avg, 5), np.percentile(avg, 100)))
    #im = ax.imshow(avg, cmap=CMAP, norm=LogNorm(vmin, vmax))
    colors = ['#FFFFFF', '#FFEDA0', '#FEB24C', '#F03B20', '#000000']
    #CMAP = LinearSegmentedColormap.from_list('custom_fall_black', colors, N=256)
    im = ax.imshow(avg, cmap=CMAP, norm=LogNorm(vmin, vmax))
    #im = ax.imshow(avg, cmap=CMAP, vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    # 标注采样点 (参考 plot_all8 风格)
    # L-C 点 (红色圆圈)
    ax.scatter([M_bin, S_bin], [S_bin, M_bin], s=30, edgecolors='blue', facecolors='none', label='L-C')
    ax.scatter([M_bin, E_bin], [E_bin, M_bin], s=30, edgecolors='blue', facecolors='none')
    
    # C-C 点 (洋红色方块)
    #ax.scatter([E_bin, S_bin], [S_bin, E_bin], s=30, edgecolors='magenta', facecolors='none', marker='s', label='A-A')
    ax.scatter([E_bin, S_bin], [S_bin, E_bin], s=30, edgecolors='cyan', facecolors='none', marker='s', label='A-A')
    
    # 设置标题显示量化结果
    title_str = f"{title}\nA-A: {val_cc:.3f} | L-A: {val_lc:.3f}\nRatio: {ratio:.3f}"
    ax.set_title(title_str, fontsize=8, fontweight='bold')
    
    ax.axis('off')
    plt.savefig(f'{output_prefix}_quant.new.svg', bbox_inches='tight', dpi=500)
    print(f"Quantification: {title_str}")
    
    return val_cc, val_lc, ratio

avg = joblib.load(f'{indir}/sim_rescale_agg.pkl')
avg = distance_norm(avg)
#plot_and_quantify(avg, 'sim_result')
plot_and_quantify(avg, outdir)
