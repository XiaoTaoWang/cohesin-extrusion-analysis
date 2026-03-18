import cooler, matplotlib, joblib
matplotlib.use('Agg')
from multiprocess import Pool
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from skimage import transform
from matplotlib.colors import LogNorm

cmap = LinearSegmentedColormap.from_list('interaction',
                                        ['#FFFFFF','#FFDFDF','#F70000'])

def parse_regions(fil, remove_chr=True):

    anchors = {}
    loadings = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e, label = line.rstrip().split()
            if remove_chr:
                c = c.lstrip('chr')
            s, e = int(s), int(e)
            p = (s + e) // 2
            key = label.split('_')[0]
            if not label.split('_')[1].startswith('M'):
                if not key in anchors:
                    anchors[key] = []
                anchors[key].append((c, p))
            else:
                if not key in loadings:
                    loadings[key] = []
                loadings[key].append((c, p))
    
    triples = {} # by chrom
    for key in anchors:
        S, E = anchors[key]
        c = S[0]
        for M in loadings[key]:
            if not c in triples:
                triples[c] = []
            triples[c].append((c, S[1], M[1], E[1]))
    
    return triples

def collect_submatrices(uri, triples, balance=False, minsize=10):
    
    lib = cooler.Cooler(uri)
    res = lib.binsize
    fixed_length = 80
    fixed_shape = (fixed_length//2, fixed_length//2)

    arrs = []
    filtered_triples = []
    for c in triples:
        if not c in lib.chromnames:
            continue

        M = lib.matrix(balance=balance, sparse=True).fetch(c).tocsr()
        for _, s, m, e in triples[c]:
            s = s // res
            m = m // res
            e = e // res
            el = (m - s) // 2
            er = (e - m) // 2
            if (m - s < minsize) or (e - m < minsize):
                continue

            if (s - el >= 0) and (e + er + 1 <= M.shape[0]):
                sub1 = M[s-el:m, s-el:m].toarray().astype(float)
                sub1[np.isnan(sub1)] = 0
                sub2 = M[s-el:m, m:e+er+1].toarray().astype(float)
                sub2[np.isnan(sub2)] = 0
                sub3 = M[m:e+er+1, s-el:m].toarray().astype(float)
                sub3[np.isnan(sub3)] = 0
                sub4 = M[m:e+er+1, m:e+er+1].toarray().astype(float)
                sub4[np.isnan(sub4)] = 0
                resized_sub1 = transform.resize(sub1, output_shape=fixed_shape)
                resized_sub2 = transform.resize(sub2, output_shape=fixed_shape)
                resized_sub3 = transform.resize(sub3, output_shape=fixed_shape)
                resized_sub4 = transform.resize(sub4, output_shape=fixed_shape)
                #resized_sub1 = resized_sub1 / resized_sub1.mean()
                #resized_sub2 = resized_sub2 / resized_sub2.mean()
                #resized_sub3 = resized_sub3 / resized_sub3.mean()
                #resized_sub4 = resized_sub4 / resized_sub4.mean()
                big = np.zeros((fixed_length, fixed_length))
                big[:fixed_shape[0], :fixed_shape[0]] = resized_sub1
                big[:fixed_shape[0], fixed_shape[0]:] = resized_sub2
                big[fixed_shape[0]:, :fixed_shape[0]] = resized_sub3
                big[fixed_shape[0]:, fixed_shape[0]:] = resized_sub4
                big = big / big.mean()

                arrs.append(big)
                filtered_triples.append((c, s, m, e))
    
    return arrs, filtered_triples

def pileup(arrs):

    pool = np.r_[arrs]
    avg = pool.mean(axis=0)

    return avg

uri = 'K562-MNase-R1-filtered.mcool::resolutions/2000'
outfig = 'test-wt.svg'

n_process = 6
remove_chr = False
balance = False

triples = parse_regions('K562-loading-sites-within-convergent-loops.bed',
                        remove_chr=remove_chr)

chroms = list(triples)
arrs, filtered_triples = collect_submatrices(uri, triples, balance=balance)
avg = pileup(arrs)

joblib.dump(avg, outfig.replace('.svg', '.pkl'))

#vmin = avg.min()
#vmax = avg.max()
vmin = 0.02
vmax = 0.8
size = (2.2, 2)
fig = plt.figure(figsize=size)
width = 0.7; Left = 0.1
HB = 0.1; HH = width * size[0] / size[1]
ax = fig.add_axes([Left, HB, width, HH])
#sc = ax.imshow(avg, cmap='seismic', vmax=vmax, vmin=vmin)
#sc = ax.imshow(avg, cmap='afmhot_r', norm=LogNorm(vmin=0.02, vmax=10))
sc = ax.imshow(avg, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
               labelbottom=False, labeltop=False, labelleft=False, labelright=False)
for spine in ['right', 'top', 'bottom', 'left']:
    ax.spines[spine].set_linewidth(0.9)
ax = fig.add_axes([Left+width+0.04, 0.72, 0.03, 0.15])
fig.colorbar(sc, cax=ax, ticks=[vmin, vmax], format='%.3g')
plt.savefig(outfig, dpi=500, bbox_inches='tight')
plt.close()