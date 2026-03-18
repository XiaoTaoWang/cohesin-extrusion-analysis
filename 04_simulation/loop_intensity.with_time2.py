import re
import sys
import glob
import joblib
import matplotlib.pyplot as plt
import pandas as pd
from loop_intensity import extract_loop_intensities
#from loop_intensity1 import extract_loop_intensities

# 假设之前的 extract_loop_intensities 函数已定义

def plot_loop_intensity_over_time(cool_files, bedpe_path, target_chrom, region_start, region_end):
    # 1. 提取时间点并去重
    # 使用正则表达式匹配文件名中 't' 后面数字的部分，例如 't23.2000' -> 23.2
    file_data = []
    for f in cool_files:
        time_match = re.search(r'_t([\d\.]+)\.', f)
        if time_match:
            time_val = int(float(time_match.group(1)))
            file_data.append({'time': time_val, 'path': f})
    print(file_data) 
    # 转换为 DataFrame 并根据时间点去重（保留每个时间点的一个路径即可）
    df_files = pd.DataFrame(file_data).drop_duplicates(subset='time').sort_values('time')
    
    times = []
    intensity_sums = []
    
    print(f"开始处理 {len(df_files)} 个唯一时间点...")
    
    # 2. 循环处理每个 cool 文件
    for _, row in df_files.iterrows():
        t = row['time']
        path = row['path']
        
        # 调用之前定义的函数提取该时间点的所有 Loop 强度
        # 注意：这里的 window_bins 可以根据需要调整，例如设为 2
        #df_res = extract_loop_intensities(path, bedpe_path, target_chrom, region_start, region_end, window_bins=2, bg_bins=5)
        df_res = extract_loop_intensities(path, bedpe_path, target_chrom, region_start, region_end, window_bins=2)#, bg_bins=10)
        
        # 计算该时间点所有 Loop 强度的总和
        total_intensity = df_res['intensity'].sum()
        
        times.append(t)
        intensity_sums.append(total_intensity)
        print(f"Time: {t}, Total Intensity Sum: {total_intensity:.4f}")
    joblib.dump([times, intensity_sums], f'{group}.loop_intensity.pkl')
    '''
    # 3. 绘制曲线
    plt.figure(figsize=(10, 6))
    plt.plot(times, intensity_sums, marker='o', linestyle='-', color='teal', linewidth=2)
    plt.xlabel('Simulation Time (t)', fontsize=12)
    plt.ylabel('Sum of Loop Intensities', fontsize=12)
    plt.title(f'Loop Intensity Evolution ({target_chrom}:{region_start}-{region_end})', fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.savefig(f'{group}.loop_intensity.svg', bbox_inches='tight')
    '''
# --- 参数设置 ---
#group = 'wild_stability0.01'
group = sys.argv[1]
thr = int(sys.argv[2])

cool_list = glob.glob(f'{group}/cool_with_time_u*/*cool')

def filter_by_t(files, n):
    """筛选 t < n 的文件"""
    pattern = r'_upto_t(\d+).2000\.cool$'
    result = []
    for f in files:
        match = re.search(pattern, f)
        if match:
            t_value = int(match.group(1))
            if t_value < n:
                result.append(f)
    return result

cool_list = filter_by_t(cool_list, thr)


# 调用绘图（请确保 region 参数与你之前的 fetch 区域对应）
plot_loop_intensity_over_time(
    cool_list, 
    bedpe_path='small_loop.bedpe',
    target_chrom='chr4', 
    region_start=36702000,
    region_end=38702000
)
'''
    region_start=37302000, 
    region_end=38102000
'''

