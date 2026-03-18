import numpy as np
from collections import Counter
import sys, bisect, joblib
import matplotlib.pyplot as plt

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

def check_in(p, List):

    hit = None
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        
        hit = q
    
    return hit

def assign_complexes_to_TADs(tads, fil):

    pool = {}
    with open(fil, 'r') as source:
        for i, line in enumerate(source):
            if i == 0:
                continue
            parse = line.rstrip().split()
            intervals, chrom = read_frags(parse[-1])
            if not chrom in tads:
                continue
            cache = {}
            for chrom, s, e in intervals:
                hit = check_in((s, e), tads[chrom])
                if hit is None:
                    continue
                key = (chrom, hit[0], hit[1])
                if not key in cache:
                    cache[key] = []
                cache[key].append(s)
            # the fragment midpoints should already be sorted
            for key in cache:
                if not key in pool:
                    pool[key] = []
                pool[key].append(cache[key])
    
    return pool

def count_frags(complex_by_tad):

    pool = np.zeros((len(complex_by_tad), 5))
    for i, key in enumerate(complex_by_tad):
        Ns = []
        for C in complex_by_tad[key]:
            if len(C) > 6:
                Ns.append(6)
            else:
                Ns.append(len(C))
        counts = Counter(Ns)
        for j in counts:
            pool[i, j-2] = counts[j] / len(Ns)
        
    avg = pool.mean(axis=0)
    std = pool.std(axis=0)
    stderr = std / np.sqrt(len(pool))

    return avg, stderr

def count_frags_new(complex_by_tad):

    pool = np.zeros(5)
    Ns = []
    for i, key in enumerate(complex_by_tad):
        for C in complex_by_tad[key]:
            if len(C) == 1:
                continue
            if len(C) > 6:
                Ns.append(6)
            else:
                Ns.append(len(C))
    counts = Counter(Ns)
    for j in counts:
        pool[j-2] = counts[j] / len(Ns)

    return pool, counts


tads = parse_tads(sys.argv[1])
complex_by_tad = assign_complexes_to_TADs(tads, sys.argv[2])
# average proportion of 2, 3, 4, 5, 6 + frags within a TAD
ratios, counts = count_frags_new(complex_by_tad)
print(counts)
#joblib.dump([ratios, stderr], 'average-proportion-TADs.pkl')

#ratios, stderr = joblib.load('average-proportion-TADs.pkl')
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111)
labels = ['2 frags', '3 frags', '4 frags', '5 frags', '6+ frags']
b = ax.bar(labels, ratios, width=0.6, color='#f7a04d')
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.xaxis.set_tick_params(width=1.5, labelsize=16)
ax.yaxis.set_tick_params(width=1.5, labelsize=16)
ax.set_ylabel('Average proportion', fontsize=18)
plt.savefig('frags-within-TADs.pdf', dpi = 300, bbox_inches='tight')
plt.close()
