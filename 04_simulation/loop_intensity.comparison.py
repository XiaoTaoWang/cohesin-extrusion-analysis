import joblib
import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

def compare_loop_intensity(group_names):
    #fig = plt.figure(figsize=(1.7, 1.2))
    fig = plt.figure(figsize=(3, 1))
    ax = fig.add_subplot(111)
   
    ls = [] 
    for group in group_names:
        # 自动根据 group 名字读取 pkl 文件
        filename = f'{group}.loop_intensity.pkl'
        try:
            times, intensity_sums = joblib.load(filename)
            #if 'wild' in group:
            if 'speed1' in group:
                color='#7FC97F'
            else:
                color='#FDC086'
            #plt.plot(times, intensity_sums, label=group, marker='o', markersize=4)
            def smooth_data(data, window_size=1):
                return np.convolve(data, np.ones(window_size)/window_size, mode='same')
            #intensity_sums = smooth_data(intensity_sums, window_size=3)
            #intensity_sums = savgol_filter(intensity_sums, window_length=5, polyorder=2)
            #intensity_sums = gaussian_filter1d(intensity_sums, sigma=2)
            #'''
            spline = make_interp_spline(times, intensity_sums, k=3)
            times = np.linspace(times[0], times[-1], 100)
            intensity_sums = spline(times)
            #'''
            if 'wild' in group:
                l, = ax.plot(times, intensity_sums, color=color, linewidth=1)
            else:
                l, = ax.plot(times, intensity_sums, color=color, linewidth=1, linestyle=':')
            print(times) 
            ls.append(l)
        except FileNotFoundError:
            print(f"错误: 找不到文件 {filename}")
    #plt.yscale('log')
    ax.set_xlabel('Simulation Time', fontsize=7)
    ax.set_ylabel('Loop Intensity', fontsize=7)
    #ax.legend(ls, ['One-sided', 'Two-sided'], fontsize=6) #, loc='lower right')
    ax.legend(ls, ['One-sided-speed2', 'One-sided-speed1'], fontsize=6)
    ax.xaxis.set_tick_params(width=0.8, labelsize=6)
    ax.yaxis.set_tick_params(width=0.8, labelsize=6)
    #xticks = list(range(5, 131, 25))
    #ax.set_xticks(xticks)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['left'].set_linewidth(0.8)
    '''
    plt.xlabel('Simulation Time (t)')
    plt.ylabel('Sum of Loop Intensities')
    plt.title('Loop Intensity Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    '''
    # 保存对比图
    output_name = "_vs_".join(group_names) + ".comparison.pdf"
    plt.savefig(output_name, bbox_inches='tight', dpi=500)
    plt.show()

# --- 调用示例 ---
# 只需要输入你想比较的 group 名字即可
#compare_loop_intensity(['wild_stability0.01', 'control_speed0.5'])
compare_loop_intensity(['wild_stability0.01', 'wild_speed1'])
