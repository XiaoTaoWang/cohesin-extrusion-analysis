import pandas as pd
import cooler
import numpy as np

def extract_loop_intensities(clr_path, bedpe_path, target_chrom, region_start, region_end, window_bins=0):
    """
    提取指定区域内所有 Loop 的强度值。
    强度定义为：Loop 中心 bin 及其周围 (window_bins) 范围内的像素总和。
    """
    # 1. 加载数据
    clr = cooler.Cooler(clr_path)
    res = clr.binsize
    
    # 读取并筛选 Loop (参考 plot_loops.py 逻辑)
    loops = pd.read_csv(bedpe_path, sep='\t', comment='#', header=None)
    loops = loops[loops[6] == 'convergent'].iloc[:, :6]
    loops.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    
    # 筛选区域内的 Loop
    mask = (loops['chrom1'] == target_chrom) & (loops['start1'] >= region_start) & (loops['end1'] <= region_end) & \
           (loops['chrom2'] == target_chrom) & (loops['start2'] >= region_start) & (loops['end2'] <= region_end)
    filtered_loops = loops[mask].copy()
    
    # 2. 获取该区域的稀疏矩阵
    # 使用 clr.matrix().fetch 获取整个区域以提高批量查询速度
    #mat = clr.matrix(balance=True).fetch(f"{target_chrom}:{region_start}-{region_end}")
    mat = clr.matrix(balance=False, sparse=False).fetch(clr.chromnames[0]).astype(float)
    mat /= np.median(np.diag(mat, 2))

    #mat = mat[300:700, 300:700]
    
    intensities = []
    
    for _, row in filtered_loops.iterrows():
        # 计算中心位置 (以 bin 为单位，相对于 region_start)
        y_center = ((row['start1'] + row['end1']) / 2 - region_start) // res
        x_center = ((row['start2'] + row['end2']) / 2 - region_start) // res
        
        # 定义提取窗口 (以中心为圆心或方块)
        y_start, y_end = int(y_center - window_bins), int(y_center + window_bins + 1)
        x_start, x_end = int(x_center - window_bins), int(x_center + window_bins + 1)
        
        # 边界检查
        if y_start < 0 or x_end > mat.shape[0]:
            intensities.append(np.nan)
            continue
            
        # 提取局部窗口并求和
        loop_window = mat[y_start:y_end, x_start:x_end]
        loop_intensity = np.nansum(loop_window)
        intensities.append(loop_intensity)
        
    filtered_loops['intensity'] = intensities
    return filtered_loops

if __name__ == '__main__':
    # --- 使用示例 ---
    clr_file = 'wild_stability0.01/cool_with_time_u999/cr139_upto_t999.2000.cool'
    loop_file = 'small_loop.bedpe'
    chrom = 'chr4'
    start, end = 37302000, 38102000

    df_intensities = extract_loop_intensities(clr_file, loop_file, chrom, start, end)
    print(df_intensities[['start1', 'start2', 'intensity']])
