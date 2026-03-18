import sys, bisect
import numpy as np
from palettable.colorbrewer.sequential import OrRd_5
import matplotlib.pyplot as plt
from collections import Counter

colors = OrRd_5.hex_colors

def parse_tads(fil):
    
    D = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom, start, end = parse[0], int(parse[1]), int(parse[2])
            if not chrom in D:
                D[chrom] = []
            D[chrom].append((start, end))
    
    for c in D:
        D[c].sort()
    
    return D

def check_in(p, List):

    cache = []
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        
        cache.append(q)
    
    return cache

def read_frags(string):

    frags = string.split(';')
    intervals = []
    for f in frags:
        s, e = f.split('-')
        chrom, s = s.split(':')
        s = int(s)
        e = int(e.split('(')[0])
        mid = (s + e) // 2
        intervals.append((chrom, mid, mid+1))
    
    return intervals, chrom

def work_core(string, tads):

    intervals, chrom = read_frags(string)
    hit = set()
    for chrom, s, e in intervals:
        cache = check_in((s, e), tads[chrom])
        hit.update(cache)
    
    return hit

def count_complexes(fil, tads, ntads=[(1,2),(2,3),(3,4),(4,5),(5,np.inf)]):

    D = {}
    total = 0
    with open(fil, 'r') as source:
        for i, line in enumerate(source):
            if i == 0:
                continue
            parse = line.rstrip().split()
            intervals, chrom = read_frags(parse[-1])
            if not chrom in tads:
                continue
            hit = []
            for chrom, s, e in intervals:
                cache = check_in((s, e), tads[chrom])
                hit.extend(cache)
            tad_freq = Counter(hit)
            n = len([tad for tad in tad_freq if tad_freq[tad]>=3]) # a complex must have 3 fragments in each TAD
            for li, ri in ntads:
                if li <= n < ri:
                    if not (li, ri) in D:
                        D[(li, ri)] = 0
                    D[(li, ri)] += 1
                    total += 1
                    break
    
    for k in D:
        D[k] = D[k] / total
    
    return D

tads = parse_tads(sys.argv[1])
D = count_complexes(sys.argv[2], tads)
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)
x = sorted(D)
print(x)
ratios = [D[i] for i in x][::-1]
print(ratios)
colors = colors[::-1]
patches, texts, autotexts = ax.pie(ratios, colors=colors, autopct='%1.1f%%',
                                startangle=90, labeldistance=1.03, pctdistance=0.6)
plt.savefig('complexes-TAD-span-new.svg', dpi=300, bbox_inches='tight')
plt.close()
