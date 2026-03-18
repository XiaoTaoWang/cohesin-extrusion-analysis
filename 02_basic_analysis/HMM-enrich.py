import bisect, random, os, joblib
import numpy as np
from scipy.special import ndtr
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import LogNorm

random.seed(10)

def parseChromLens(chromfil):
    
    chromsizes = {}
    with open(chromfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0].lstrip('chr')
            if '_' in chrom:
                continue
            chromsizes[chrom] = int(parse[1])
            
    return chromsizes

def parseGapFile(gapfil):
    
    gaps = {}
    with open(gapfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[1].lstrip('chr')
            if '_' in chrom:
                continue
            if not chrom in gaps:
                gaps[chrom] = []
            sb = int(parse[2])
            eb = int(parse[3])
            gaps[chrom].append([sb, eb+1])

    for c in gaps:
        gaps[c].sort()
    
    return gaps

def generate_control_list(true_loci, chromfil, gapfil):

    gaps = parseGapFile(gapfil)
    chromsizes = parseChromLens(chromfil)

    random_list = []
    visited = set()
    for c, s, e in true_loci:
        interval = e - s
        check = True
        while check:
            rp = random.randint(0, chromsizes[c]-interval-1)
            tmp = [rp, rp+interval]
            cache = overlap_peaks(tmp, gaps[c])
            key = (c, rp, rp+interval)
            if (len(cache)==0) and (not key in visited):
                random_list.append(list(key))
                visited.add(key)
                check = False
    
    return random_list

def parsePeaks(fil):

    peaks = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0].lstrip('chr')
            if '_' in chrom:
                continue
            if not chrom in peaks:
                peaks[chrom] = []
            sb = int(parse[1])
            eb = int(parse[2])
            peaks[chrom].append([sb, eb])
    
    for c in peaks:
        peaks[c].sort()
    
    return peaks

def overlap_peaks(p, List):

    cache = set()
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        cache.add(tuple(q))
    
    return cache

def parse_hmm(fil, states_map):

    hmm = {}
    with open(fil, 'r') as source:
        for line in source:
            if line.startswith('#'):
                continue
            parse = line.rstrip().split('\t')
            chrom, s, e, label = parse[:4]
            chrom = chrom.lstrip('chr')
            s, e = int(s), int(e)
            label = states_map[label]
            if not label in hmm:
                hmm[label] = {}
            if not chrom in hmm[label]:
                hmm[label][chrom] = []
            hmm[label][chrom].append([s, e])
    
    for label in hmm:
        for chrom in hmm[label]:
            hmm[label][chrom].sort()
    
    return hmm

def check_in(p, List, mode='binary', min_len=15000):

    interval = p[1] - p[0]
    if interval < min_len:
        half = (min_len - interval) // 2
        tmp = [max(0, p[0]-half), p[1]+half]
        p = tmp
    interval = p[1] - p[0]

    cache = set()
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        minp, maxp = p[0], p[1]
        if q[0] > minp:
            minp = q[0]
        if q[1] < maxp:
            maxp = q[1]
        cache.add((minp, maxp))
    
    score = 0
    if mode=='binary':
        if len(cache):
            score = 1
    else:
        score = 0
        for s, e in cache:
            score += (e-s)
        score = score / interval

    
    return score

###### load and parse CTCF and cohesin peaks
fil1 = 'GM12878-CTCF-pooled_comp_sing_FDR_0.2_PASS_thresh70_merge3kb_peaks.bed'
fil2 = 'GM12878-cohesin-pooled_comp_sing_FDR_0.2_PASS_thresh400_merge3kb_peaks.bed'
ctcf_peaks = parsePeaks(fil1)
cohesin_peaks = parsePeaks(fil2)
ctcf_unique = set()
cohesin_unique = set()
common_peaks = set()
for c in ctcf_peaks:
    for p in ctcf_peaks[c]:
        cache = overlap_peaks(p, cohesin_peaks[c])
        if len(cache):
            for q in cache:
                common_peaks.add((c, q[0], q[1])) # overlap use the locations of cohesin peaks
        else:
            ctcf_unique.add((c, p[0], p[1]))

for c in cohesin_peaks:
    for p in cohesin_peaks[c]:
        k = (c, p[0], p[1])
        if not k in common_peaks:
            cohesin_unique.add(k)

peak_pool = [ctcf_unique, common_peaks, cohesin_unique]
peak_labels = ['CTCF Specific', 'Common Peaks', 'Cohesin Specific']

###### generate random loci
chromfil = 'hg38.chrom.sizes'
gapfil = 'hg38.gap.txt'
random_pool = []
print('generating random peaks ...')
for true_loci in peak_pool:
    tmp_pool = [true_loci]
    for i in range(100):
        random_loci = generate_control_list(true_loci, chromfil, gapfil)
        tmp_pool.append(random_loci)
    random_pool.append(tmp_pool)

###### load HMM states
states_map = {
    '1_Active_Promoter':'Active Promoter',
    '2_Weak_Promoter':'Weak Promoter',
    '3_Poised_Promoter':'Poised Promoter',
    '4_Strong_Enhancer':'Strong Enhancer',
    '5_Strong_Enhancer':'Strong Enhancer',
    '6_Weak_Enhancer':'Weak Enhancer',
    '7_Weak_Enhancer':'Weak Enhancer',
    '8_Insulator':'Insulator',
    '9_Txn_Transition':'Txn Transition',
    '10_Txn_Elongation':'Txn Elongation',
    '11_Weak_Txn':'Weak Txn',
    '12_Repressed':'Repressed',
    '13_Heterochrom/lo':'Heterochromatin',
    '14_Repetitive/CNV':'Repetitive/CNV',
    '15_Repetitive/CNV':'Repetitive/CNV'
}
states_code = [
    'Active Promoter',
    'Weak Promoter',
    'Poised Promoter',
    'Strong Enhancer',
    'Weak Enhancer',
    'Insulator',
    'Txn Transition',
    'Txn Elongation',
    'Weak Txn',
    'Repressed',
    'Heterochromatin',
    'Repetitive/CNV'
]


states = parse_hmm('BASIC_chromHMM_GM12878_hg38_Broad_UCSC.txt.tsv', states_map)

###### check hmm states within peaks
mode = 'proportion' # binary or proportion
min_len = 0 # minimum interval of a peak
score_pool = []
queue = list(zip(*random_pool)) # the 1st one is the true loci
print('calculate HMM occupation for each peak loci')
count = 0
for q in queue:
    print(count)
    arr = np.zeros((len(q), len(states_code)))
    for i, loci in enumerate(q):
        for j, s_n in enumerate(states_code):
            state = states[s_n]
            scores = []
            for c, s, e in loci:
                if not c in state:
                    scores.append(0)
                    continue
                scores.append(check_in([s, e], state[c], mode=mode, min_len=min_len))
            arr[i, j] = np.mean(scores)
    score_pool.append(arr)
    count += 1

score_pool = np.r_[score_pool]
joblib.dump(score_pool, 'intermediate.pkl')

###### calculate p-values by comparing with random peaks
true_scores = score_pool[0]
random_scores = score_pool[1:]
parr = np.ones_like(true_scores)
farr = np.ones_like(true_scores)
for i in range(true_scores.shape[0]):
    for j in range(true_scores.shape[1]):
        t_s = true_scores[i, j]
        r_s = random_scores[:, i, j]
        u = np.mean(r_s)
        std = np.std(r_s)
        z_s = (t_s - u) / std
        pvalue = 1 - ndtr(z_s)
        parr[i, j] = pvalue
        farr[i, j] = t_s / u

fig = plt.figure(figsize=(10, 2))
ax = fig.add_subplot(111)
data = pd.DataFrame(farr, columns=states_code, index=peak_labels)
cg = sns.heatmap(data, cmap='YlGnBu_r', annot=True, ax=ax, norm=LogNorm(vmin=farr.min(), vmax=farr.max()))
plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=12)
plt.setp(cg.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=12)
cg.set_clip_on(False)
plt.savefig('hmm-fc.pdf', dpi=300, bbox_inches='tight')
plt.close()
