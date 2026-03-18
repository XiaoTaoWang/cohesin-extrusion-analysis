
# ── Path configuration ─────────────────────────────────────────────
# Set these to match your local environment before running.
import os
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_ROOT    = os.environ.get("DATA_ROOT", os.path.dirname(PROJECT_ROOT))
GENOME_ROOT  = os.environ.get("GENOME_ROOT", os.path.join(os.path.dirname(DATA_ROOT), "genome"))
# ──────────────────────────────────────────────────────────────────

import bisect

def mapPeak(p, peakpool):
    
    # p: (chrom, start, end), loop anchor
    # peakpool: either ChIP-Seq peaks or CTCF motifs
    cache = set()
    if not p[0] in peakpool:
        return cache

    bychrom = peakpool[p[0]]
    idx = max(0, bisect.bisect(bychrom, p[1:])-1)
    for q in bychrom[idx:]:
        if q[1] <= p[1]:
            continue
        if q[0] >= p[2]:
            break
        cache.add(tuple(q))

    return cache

def load_loops(infil):

    collect = []
    with open(infil, 'r') as source:
        for line in source:
            c1, s1, e1, c2, s2, e2 = line.rstrip().split()
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            collect.append((c1, s1, e1, c2, s2, e2))
    
    return sorted(collect)

def parseMotif(fil):
    
    data = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if '_' in chrom:
                continue
            if chrom in data:
                data[chrom].append([int(parse[1]), int(parse[2]), parse[3]])
            else:
                data[chrom] = [[int(parse[1]), int(parse[2]), parse[3]]]
    for chrom in data:
        data[chrom].sort()
    
    return data

loops = load_loops(f'{PROJECT_ROOT}/mESC/loading_site_define/ENCFF550QMW_6col.bedpe')
motifs = parseMotif('mESC.CTCF-motifs.mm10.txt')

annot = []
for c1, s1, e1, c2, s2, e2 in loops:
    a1 = [c1, s1, e1]
    a2 = [c2, s2, e2]
    cache1 = set([m[-1] for m in mapPeak(a1, motifs)])
    cache2 = set([m[-1] for m in mapPeak(a2, motifs)])
    if '+' in cache1 and '-' in cache2:
        label = 'convergent'
    else:
        label = '.'
    
    annot.append([c1, s1, e1, c2, s2, e2, label])

with open('mESC_micro-c_mustache_loop.annotated.bedpe', 'w') as out:
    for loop in annot:
        out.write('\t'.join(list(map(str, loop)))+'\n')