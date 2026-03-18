import cooler, matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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
    
    units = {} # store by chrom
    for key in anchors:
        S, E = sorted(anchors[key])
        c = S[0]
        if not c in units:
            units[c] = []
        if key in loadings:
            if len(loadings[key]) >= 2:
                units[c].append([[S, E], sorted(loadings[key])])
    
    return units

def collect_submatrices(uri, units, balance=False, w=10, min_dis=30):
    
    lib = cooler.Cooler(uri)
    res = lib.binsize

    l_to_l = []
    l_to_a = []
    a_to_a = []
    for c in units:
        if not c in lib.chromnames:
            continue

        M = lib.matrix(balance=balance, sparse=True).fetch(c).tocsr()
        for anchors, loadings in units[c]:
            S, E = anchors
            # anchor to anchor
            xi = S[1] // res
            yi = E[1] // res
            if yi - xi < min_dis:
                continue

            if (xi - w >= 0) and (yi + w + 1 <= M.shape[0]):
                sub = M[xi-w:xi+w+1, yi-w:yi+w+1].toarray().astype(float)
                sub[np.isnan(sub)] = 0
                if sub.max() == sub.min():
                    continue
                sub = sub / sub.mean()
                a_to_a.append(sub)
            
            # loading to loading
            for p1 in loadings:
                for p2 in loadings:
                    if p1 < p2:
                        xi = p1[1] // res
                        yi = p2[1] // res
                        if yi - xi < min_dis:
                            continue

                        if (xi - w >= 0) and (yi + w + 1 <= M.shape[0]):
                            sub = M[xi-w:xi+w+1, yi-w:yi+w+1].toarray().astype(float)
                            sub[np.isnan(sub)] = 0
                            if sub.max() == sub.min():
                                continue
                            sub = sub / sub.mean()
                            l_to_l.append(sub)
            
            # loading to anchor
            xi = S[1] // res
            for p in loadings:
                yi = p[1] // res
                if yi - xi < min_dis:
                    continue

                if (xi - w >= 0) and (yi + w + 1 <= M.shape[0]):
                    sub = M[xi-w:xi+w+1, yi-w:yi+w+1].toarray().astype(float)
                    sub[np.isnan(sub)] = 0
                    if sub.max() == sub.min():
                        continue
                    sub = sub / sub.mean()
                    l_to_a.append(sub)
            
            yi = E[1] // res
            for p in loadings:
                xi = p[1] // res
                if yi - xi < min_dis:
                    continue

                if (xi - w >= 0) and (yi + w + 1 <= M.shape[0]):
                    sub = M[xi-w:xi+w+1, yi-w:yi+w+1].toarray().astype(float)
                    sub[np.isnan(sub)] = 0
                    if sub.max() == sub.min():
                        continue
                    sub = sub / sub.mean()
                    l_to_a.append(sub)
           
    
    return a_to_a, l_to_a, l_to_l

def pileup(arrs):

    pool = np.r_[arrs]
    avg = pool.mean(axis=0)
    avg = avg / avg.mean()

    return avg

uri = 'ChIA-PET_hg38_HCT116-mAC-ctrl_none-likehic_LHH0146V_novaseq_pairs.5k.cool'
coord_fil = 'HCT116-cohesin-specific-regions_20230710_filt_25-75percentile_20231030.bed'
outpre = 'HCT116'
balance = False
outfig1 = '{0}.anchor-to-anchor.{1}.svg'.format(outpre, balance)
outfig2 = '{0}.loading-to-anchor.{1}.svg'.format(outpre, balance)
outfig3 = '{0}.loading-to-loading.{1}.svg'.format(outpre, balance)
ws = 10

units = parse_regions(coord_fil)
a_to_a, l_to_a, l_to_l = collect_submatrices(uri, units, w=ws, balance=balance)
queue = [
    [a_to_a, outfig1, 'anchor to anchor'],
    [l_to_a, outfig2, 'loading to anchor'],
    [l_to_l, outfig3, 'loading to loading']
]

for arrs, outfig, label in queue:
    avg = pileup(arrs)
    obs = avg[ws, ws]
    background = []
    for i in range(2*ws+1):
        for j in range(2*ws+1):
            if (i==ws) and (j==ws):
                continue
            background.append(avg[i,j])
    background = np.r_[background]
    z = (obs - background.mean()) / background.std()
    print('{0}: {1}'.format(label, z))

    vmin = 0.8
    vmax = 5
    size = (2.2, 2)
    fig = plt.figure(figsize=size)
    width = 0.7; Left = 0.1
    HB = 0.1; HH = width * size[0] / size[1]
    ax = fig.add_axes([Left, HB, width, HH])
    #sc = ax.imshow(avg, cmap='seismic', vmax=vmax, vmin=vmin)
    #sc = ax.imshow(avg, cmap=cmap, vmin=vmin, vmax=vmax)
    sc = ax.imshow(avg, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in ['right', 'top', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(0.9)
    ax = fig.add_axes([Left+width+0.04, 0.72, 0.03, 0.15])
    fig.colorbar(sc, cax=ax, ticks=[vmin, vmax], format='%.3g')
    plt.savefig(outfig, dpi=500, bbox_inches='tight')
    plt.close()
