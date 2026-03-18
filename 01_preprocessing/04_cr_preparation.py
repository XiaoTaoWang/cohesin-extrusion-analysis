import pandas as pd

def process_loading_sites():
    loading_file = 'k562-chip-nipbl-overlap-rad21-rnapii-h3k27ac-h3k4me1.bed'
    loop_file = 'distance_Convergent_Loops_top10k.bedpe'
    
    out_all = 'K562-cohesin-specific-regions_all.bed'
    out_25_75 = 'K562-cohesin-specific-regions_filt_25-75percentile.bed'
    out_anchor5 = 'K562-cohesin-specific-regions_filt_anchor5percentile.bed'

    loading = pd.read_csv(loading_file, sep='\t', header=None, usecols=[0, 1, 2])
    loading.columns = ['chrom', 'start', 'end']
    loops = pd.read_csv(loop_file, sep='\t', header=None)
    
    res_all, res_25_75, res_anchor5 = [], [], []
    stats = {
        'all': {'loops': 0, 'sites': 0},
        '25_75': {'loops': 0, 'sites': 0},
        'anchor5': {'loops': 0, 'sites': 0}
    }

    for idx, row in loops.iterrows():
        chrom = row[0]
        if chrom != row[3]:
            continue
            
        left_start, left_end = int(row[1]), int(row[2])
        right_start, right_end = int(row[4]), int(row[5])
        loop_id = f"cr{idx + 1}"
        
        dist = right_start - left_end
        if dist <= 0:
            continue
            
        search_start = left_end - 4000
        search_end = right_start + 4000
        
        p5_dist = dist * 0.05
        p25 = left_end + (dist * 0.25)
        p75 = left_end + (dist * 0.75)
        
        mask = (loading['chrom'] == chrom) & (loading['end'] > search_start) & (loading['start'] < search_end)
        loop_loading = loading[mask]
        
        if loop_loading.empty:
            continue
            
        valid_all, valid_25_75, valid_anchor5 = [], [], []
        m_idx_all, m_idx_25_75, m_idx_anchor5 = 1, 1, 1
        
        for _, site in loop_loading.iterrows():
            site_mid = (site['start'] + site['end']) / 2.0
            
            # 1. All loading sites within -4kb to +4kb
            if search_start <= site_mid <= search_end:
                valid_all.append([site['chrom'], int(site['start']), int(site['end']), f"{loop_id}_M-{m_idx_all}"])
                m_idx_all += 1
                
            # 2. Loading sites within 25 - 75 percentile
            if p25 <= site_mid <= p75:
                valid_25_75.append([site['chrom'], int(site['start']), int(site['end']), f"{loop_id}_M-{m_idx_25_75}"])
                m_idx_25_75 += 1
                
            # 3. Loading sites near anchors (+/- 5 percentile)
            if (left_end - p5_dist <= site_mid <= left_end + p5_dist) or \
               (right_start - p5_dist <= site_mid <= right_start + p5_dist):
                valid_anchor5.append([site['chrom'], int(site['start']), int(site['end']), f"{loop_id}_M-{m_idx_anchor5}"])
                m_idx_anchor5 += 1

        anchors_out = [
            [chrom, left_start, left_end, f"{loop_id}_S"],
            [chrom, right_start, right_end, f"{loop_id}_E"]
        ]

        if valid_all:
            stats['all']['loops'] += 1
            stats['all']['sites'] += len(valid_all)
            res_all.extend(anchors_out + valid_all)
            
        if valid_25_75:
            stats['25_75']['loops'] += 1
            stats['25_75']['sites'] += len(valid_25_75)
            res_25_75.extend(anchors_out + valid_25_75)
            
        if valid_anchor5:
            stats['anchor5']['loops'] += 1
            stats['anchor5']['sites'] += len(valid_anchor5)
            res_anchor5.extend(anchors_out + valid_anchor5)

    pd.DataFrame(res_all).to_csv(out_all, sep='\t', header=False, index=False)
    pd.DataFrame(res_25_75).to_csv(out_25_75, sep='\t', header=False, index=False)
    pd.DataFrame(res_anchor5).to_csv(out_anchor5, sep='\t', header=False, index=False)

    print("-" * 80)
    print("Loops with at least 1 loading site between left anchor - 4kb and right anchor + 4 kb\n")
    print(f"         {stats['all']['loops']} loops encompassing {stats['all']['sites']} loading sites")
    print(f"         {out_all}")
    print("-" * 80)
    print("Loops with loading within 25 - 75 percentile of distance\n")
    print(f"         {stats['25_75']['loops']} loops encompassing {stats['25_75']['sites']} loading sites")
    print(f"         {out_25_75}")
    print("-" * 80)
    print("Loops with loading near anchors: left anchor (+/- 5 percentile) or right anchor (+/- 5 percentile) of distance\n")
    print(f"         {stats['anchor5']['loops']} loops encompassing {stats['anchor5']['sites']} loading sites near anchors")
    print(f"         {out_anchor5}")
    print("-" * 80)

if __name__ == "__main__":
    process_loading_sites()