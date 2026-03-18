
# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
import os
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import bisect

def parse_peaks(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            if not c in D:
                D[c] = []
            D[c].append((s, e))
    
    for c in D:
        D[c].sort()
    
    return D

def parse_ChIP(fil):

    peaks = []
    with open(fil, 'r') as source:
        for line in source:
            if not line.startswith('chr'):
                continue
            c, s, e = line.rstrip().split()[:3]
            peaks.append((c, int(s), int(e)))
    
    return peaks

def bisect_search(p, ref):

    cache = set()
    if not p[0] in ref:
        return cache
    
    List = ref[p[0]]
    idx = max(0, bisect.bisect(List, p[1:])-1)
    for q in List[idx:]:
        if q[1] <= p[1]:
            continue
        if q[0] >= p[2]:
            break
        cache.add((p[0],)+q)
    
    return cache

ctcf_peaks = parse_peaks(f'{DATA_ROOT}/ChIA-PET/WT_CTCF_merge/WT_CTCF_merge_peaks.narrowPeak')
cohesin_peaks = parse_peaks(f'{DATA_ROOT}/ChIA-PET/WT_RAD21_merge/WT_RAD21_merge_peaks.narrowPeak')
nipbl = parse_ChIP(f'{DATA_ROOT}/mly/K562_CUT_TAG/K562_CTCF_0h_NIPBL_CUT_TAG005/3_macs2/K562_CTCF_0h_NIPBL_CUT_TAG005_peaks.narrowPeak')
filtered = []
for c, s, e in nipbl:
    cache1 = bisect_search((c, s, e), ctcf_peaks)
    cache2 = bisect_search((c, s, e), cohesin_peaks)
    if len(cache2) and (not len(cache1)):
        filtered.append((c, s, e))

with open('K562.cohesin-loading-sites.new.bed', 'w') as out:
    for c, s, e in filtered:
        out.write('\t'.join([c, str(s), str(e)])+'\n')