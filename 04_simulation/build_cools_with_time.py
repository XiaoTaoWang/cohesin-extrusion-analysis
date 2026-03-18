import os
import glob
import numpy as np
import polychrom.hdf5_format
import polychrom.contactmaps
from polykit.analysis import contact_maps as cms

def build_cumulative_coolers(indir, outdir, interval=5):
    """
    interval: 每隔多少个 step 输出一个累积热图
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    sim_folders = sorted(glob.glob(f"{indir}/sims/folder*cr139"))
    
    for folder_path in sim_folders:
        sim_name = os.path.basename(folder_path).split('_')[-1]
        print(f"Processing evolution for {sim_name}...")
        
        try:
            URIs = polychrom.hdf5_format.list_URIs(folder_path)
            # 确保时间步严格升序
            #URIs_eq = sorted([u for u in URIs if int(u.split("::")[-1]) >= 0], 
            URIs_eq = sorted([u for u in URIs if int(u.split("::")[-1]) <= timeup], 
                            key=lambda x: int(x.split("::")[-1]))

            # --- 关键修改：起点相同，终点延伸 (累积切片) ---
            for end_idx in range(interval, len(URIs_eq) + interval, interval):
                # 确保不越界
                actual_end = min(end_idx, len(URIs_eq))
                # 始终从 0 开始切片：[0:50], [0:100], [0:150]...
                cumulative_uris = URIs_eq[0:actual_end]
                
                current_step = cumulative_uris[-1].split("::")[-1]
                
                # 生成接触图
                mrc = polychrom.contactmaps.monomerResolutionContactMapSubchains(
                    cumulative_uris,
                    np.arange(0, 10000, 1000), # monomer_per_replica
                    1000,
                    cutoff=2.3, n=10
                )

                # 命名区分时间点
                cool_uri = f'{outdir}/{sim_name}_upto_t{current_step}'
                cms.coolify(mrc, cool_uri, binsize=2000)
                print(f"Saved cumulative snapshot: {cool_uri} (steps: 0-{current_step})")
                
                if end_idx >= len(URIs_eq): break

        except Exception as e:
            print(f"Error at {sim_name}: {e}")

if __name__ == "__main__":
    import sys
    indir, outdir, timeup = sys.argv[1:4]
    timeup = int(timeup)
    #interval = timeup // 5 + 1
    interval = timeup // 10 + 1
    build_cumulative_coolers(indir, outdir, interval=interval)
