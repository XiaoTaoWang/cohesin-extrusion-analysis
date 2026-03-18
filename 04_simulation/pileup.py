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

# 从命令行获取参数
if len(sys.argv) < 3:
    print("Usage: python plot_agg3.py <COOL_DIR> <OUT_DIR>")
    sys.exit(1)

COOL_DIR, out_dir = sys.argv[1:3]

# --- 1. 配置 ---
RES = 2000
# 严格 1-bin overlap 逻辑：
# 子块 40x40，偏移量 39。
# 块1：0-39；块2：39-78。总长 79。
# 为了凑够 80，我们将最终矩阵扩充为 80，并在边缘留白或补齐。
FIXED_SUB_SHAPE = (40, 40) 
FINAL_LENGTH = 80  
MIN_SIZE = 5       

TABLE_PATH = '/home/wangxiaotao/jianglinghan/software/Targeted_Cohesin_Loading/Analysis/sequential_barriers/sandbox/203_speedup/simu_all_loop/all_ctcf_loading/K562-loading-sites-within-convergent-loops.bed'
CMAP = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#F70000'])

# --- 2. 坐标解析 ---

def parse_loading_table(fil):
    groups = {}
    with open(fil, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4: continue
            try:
                start, end, label = int(parts[1]), int(parts[2]), parts[3]
                mid_pos = (start + end) // 2
                loop_id, type_info = label.split('_')
                if loop_id not in groups:
                    groups[loop_id] = {'M': []}
                if type_info == 'S': groups[loop_id]['S'] = mid_pos
                elif type_info == 'E': groups[loop_id]['E'] = mid_pos
                elif type_info.startswith('M'): groups[loop_id]['M'].append(mid_pos)
            except (ValueError, IndexError):
                continue
    return groups

def get_dynamic_diag_norm(mat_list, target_vmax=3000):
    """
    mat_list: 包含四个 sub-matrices (r_sub1...r_sub4) 的列表
    通过搜索找到均值最接近 target_vmax 的对角线 offset
    """
    # 临时拼成一个 81x81 的大矩阵来寻找全局对角线背景
    # 这里不需要处理 overlap，只是为了计算 offset 的均值
    temp_big = np.zeros((81, 81))
    offsets = [(0, 0), (0, 40), (40, 0), (40, 40)]
    for i, (r_off, c_off) in enumerate(offsets):
        temp_big[r_off:r_off+41, c_off:c_off+41] = mat_list[i]

    # 搜索对角线 k (从 2 开始，避开极高值)
    selected_mean = 1.0
    for k in range(2, 40):
        # 提取第 k 条对角线上的值
        diag_vals = np.diag(temp_big, k)
        current_mean = np.nanmean(diag_vals)
        
        if current_mean <= target_vmax:
            selected_mean = current_mean
            # print(f"Found suitable offset k={k}, mean={selected_mean}") # 调试用
            break
    else:
        # 如果没找到，退而求其次用边缘均值
        selected_mean = np.nanmean(temp_big[temp_big > 0])

    # 对传入的四个矩阵分别进行归一化
    return [m / selected_mean for m in mat_list]

# --- 3. 核心提取与重缩放逻辑 ---

def collect_rescaled_submatrices(cool_path, loop_data):
    anchor_mid = (loop_data['S'] + loop_data['E']) // 2
    sim_start_bp = anchor_mid - 1000000
    
    c = cooler.Cooler(cool_path)
    M_full = c.matrix(balance=False, sparse=False).fetch(c.chromnames[0])
    '''
    M_full = M_full.astype(np.float64)
    diag_med = np.median(np.diag(M_full, 2))
    if diag_med > 0: M_full /= diag_med
    '''
    loop_subs = []
    
    for m_pos in loop_data['M']:
        s = (loop_data['S'] - sim_start_bp) // RES
        m = (m_pos - sim_start_bp) // RES
        e = (loop_data['E'] - sim_start_bp) // RES
        
        #el, er = (m - s) // 2, (e - m) // 2
        el, er = (m - s), (e - m)
        
        if (m - s < MIN_SIZE) or (e - m < MIN_SIZE): continue
        
        if (s - el >= 0) and (e + er + 1 <= M_full.shape[0]):
            # 四象限切割
            sub1 = M_full[s-el:m, s-el:m] #左上
            sub2 = M_full[s-el:m, m:e+er+1] #右上
            sub3 = M_full[m:e+er+1, s-el:m] #左下
            sub4 = M_full[m:e+er+1, m:e+er+1] #右下
            
            # Resize 为 40x40
            r_sub1 = transform.resize(sub1, FIXED_SUB_SHAPE, preserve_range=True, mode='edge')
            r_sub2 = transform.resize(sub2, FIXED_SUB_SHAPE, preserve_range=True, mode='edge')
            r_sub3 = transform.resize(sub3, FIXED_SUB_SHAPE, preserve_range=True, mode='edge')
            r_sub4 = transform.resize(sub4, FIXED_SUB_SHAPE, preserve_range=True, mode='edge')
            
            # 创建 80x80 画布
            big = np.zeros((FINAL_LENGTH, FINAL_LENGTH))
            count = np.zeros((FINAL_LENGTH, FINAL_LENGTH))
            
            # 严格 1-bin overlap：偏移量为 39
            # 这样索引 39 (第40个bin) 是四个矩阵唯一的交汇点
            offsets = [
                (0, 0),   # Sub1: 0-39
                (0, 39),  # Sub2: 行0-39, 列39-78
                (39, 0),  # Sub3: 行39-78, 列0-39
                (39, 39)  # Sub4: 行39-78, 列39-78
            ]
            '''
            def norm_sub(mat):
                avg_val = np.nanmean(mat)
                return mat / avg_val if avg_val > 0 else mat

            r_sub1, r_sub2, r_sub3, r_sub4 = map(norm_sub, [r_sub1, r_sub2, r_sub3, r_sub4])
            
            '''
            #r_sub1 /= 1.4
            #r_sub4 /= 1.4

            subs = [r_sub1, r_sub2, r_sub3, r_sub4]
            #subs = get_dynamic_diag_norm(subs) 
            for (r_off, c_off), s_mat in zip(offsets, subs):
                # 将 40x40 填入
                big[r_off:r_off+40, c_off:c_off+40] += s_mat
                count[r_off:r_off+40, c_off:c_off+40] += 1
            
            # 仅在 count > 0 的地方除以计数（只有索引 39 相关的线和点 count 会 > 1）
            final_sub = np.divide(big, count, out=np.zeros_like(big), where=count!=0)
           
            # 注意：此时 final_sub 的第 79 行和第 79 列（Index 79）全是 0
            # 这是因为 39 + 40 = 79。我们将结果切片为 80x80，保持原样即可。
            loop_subs.append(final_sub)
            
    return loop_subs

# --- 4. 主流程 ---

def run_pipeline():
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print("Parsing table...")
    table_info = parse_loading_table(TABLE_PATH)
    
    all_matrices = []
    cool_files = glob.glob(os.path.join(COOL_DIR, '*.cool'))
    
    for cp in cool_files:
        sim_id = os.path.basename(cp).split('.')[0]
        if sim_id in table_info:
            print(f"Processing: {sim_id}")
            subs = collect_rescaled_submatrices(cp, table_info[sim_id])
            all_matrices.extend(subs)

    if not all_matrices:
        print("No matrices collected.")
        return

    avg = np.mean(all_matrices, axis=0)
    
    # 保存数据
    save_path = os.path.join(out_dir, 'sim_rescale_agg.pkl')
    joblib.dump(avg, save_path)

    # 绘图：切掉最后一列/行全零（可选），或直接画 80x80
    plt.figure(figsize=(2.2, 2))
    ax = plt.axes([0.15, 0.15, 0.7, 0.7])
    
    vmax = np.percentile(avg[avg>0], 90) # 排除空值计算分位数
    vmin = np.percentile(avg[avg>0], 3)
    
    # 绘图时只显示有数据的 79x79 区域以保证视觉对称，或者直接显示 80x80
    # 这里建议显示 avg[:79, :79] 更加美观
    display_mat = avg[:79, :79]
    im = ax.imshow(display_mat, cmap=CMAP, vmin=vmin, vmax=vmax, interpolation='none')
    
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    plt.savefig(os.path.join(out_dir, 'sim_rescale_agg.svg'), dpi=500, bbox_inches='tight')
    print("Done.")

if __name__ == "__main__":
    run_pipeline()
