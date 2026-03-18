import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from palettable.colorbrewer.sequential import GnBu_9
import sys, cooler
from collections import defaultdict, Counter

colors = GnBu_9.hex_colors[1:]

def parse_pairwise(fil, nfrags=[(2,3),(3,4),(4,5),(5,6),(6,11),(11,21),(21,51),(51,np.inf)]):

    D = {}
    with open(fil, 'r') as source:
        for i, line in enumerate(source):
            if i == 0:
                continue
            parse = line.rstrip().split()
            mids = read_frags(parse[-1])
            n = int(parse[3])
            dis = pairwise(mids)
            for li, ri in nfrags:
                if li <= n < ri:
                    if not (li, ri) in D:
                        D[(li, ri)] = []
                    D[(li, ri)].extend(dis)
                    break
    
    D = normalized_count(D)
    
    return D

def pairwise(mids, min_apart=50000):

    dis = []
    for i in range(len(mids)-1):
        for j in range(i+1, len(mids)):
            tmp = mids[j] - mids[i]
            if tmp < min_apart:
                continue
            dis.append(tmp)
    
    return dis

def read_frags(string):

    frags = string.split(';')
    mids = []
    for f in frags:
        s, e = f.split('-')
        s = int(s.split(':')[1])
        e = int(e.split('(')[0])
        mids.append((s+e)//2)
    
    return mids

def normalized_count(Dict, binsize=50000):

    norm = {}
    for k in Dict:
        extract = Dict[k]
        binned = [i//binsize*binsize for i in extract]
        tmp = Counter(binned)
        x = sorted(tmp)
        y = np.r_[[tmp[i] for i in x]]
        y = y / y.max()
        norm[k] = [x, y]
    
    return norm

def HiC_pairwise(clr):

    pool = defaultdict(float)
    res = clr.binsize
    for c in clr.chromnames:
        M = clr.matrix(balance=False, sparse=True).fetch(c)
        N = M.shape[0]
        for i in range(1, N-1):
            diag = M.diagonal(i)
            dis = i * res
            add = diag.sum()
            pool[dis] += add
            
    x = sorted(pool)
    y = np.r_[[pool[i] for i in x]]
    y = y / y.max()

    return x, y

# pairwise
D = parse_pairwise(sys.argv[1])
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
count = 0
for li, ri in sorted(D):
    if ri - li == 1:
        name = '{0} frags'.format(li)
    elif ri == np.inf:
        name = '{0}+ frags'.format(li-1)
    else:
        name = '{0}-{1} frags'.format(li, ri-1)
    x, y = D[(li, ri)]
    ax.plot(x, y, linewidth=2, label=name, color=colors[count])
    '''
    sns.distplot(D[(li,ri)], hist=False, kde=True, color=colors[count], ax=ax,
        kde_kws={'alpha':1, 'cumulative':True, 'linewidth':2}, label=name, bins=100,
        norm_hist=True)
    '''
    count += 1

clr = cooler.Cooler('GM19239-HindIII-allReps.mcool::resolutions/50000')
x, y = HiC_pairwise(clr)
ax.plot(x, y, linewidth=2, label='Hi-C', color='#990000', ls='--')
ax.legend(frameon=False, fontsize=10, title='ChIA-Drop')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.xaxis.set_tick_params(width=1.5, labelsize=14)
ax.yaxis.set_tick_params(width=1.5, labelsize=14)
ax.set_xscale('log')
#ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
ax.set_xlabel('Pairwise distance (bp)', fontsize=18)
ax.set_ylabel('Normalized Contact Frequency', fontsize=18)
plt.savefig('ChIADrop-power-law.svg', dpi = 300, bbox_inches='tight')
plt.close()
