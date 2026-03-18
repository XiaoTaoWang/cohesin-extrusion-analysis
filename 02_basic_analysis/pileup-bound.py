import matplotlib.pyplot as plt
import numpy as np
from tadlib.hitad.aligner import *
import cooler, itertools, sys
from neoloop.visualize.bigwig import plot_y_axis

class sigTrack(object):

    def __init__(self, fil, res=25000, col=3, skiprows=0):

        sig = {}
        with open(fil, 'r') as source:
            for i, line in enumerate(source):
                if i < skiprows:
                    continue
                p = line.rstrip().split()
                if not len(p):
                    continue
                chrom = p[0]
                if not chrom in sig:
                    sig[chrom] = []
                value = float(p[col])
                sig[chrom].append(value)
        
        for c in sig:
            sig[c] = np.r_[sig[c]]
        
        self.sig = sig
        self.res = res
    
    def get_interval(self, chrom, start, end):

        si = start // self.res
        ei = end // self.res

        return self.sig[chrom][si:ei]


def triangle_plot(M, vmin=0.4, vmax=1.2, cmap='RdYlBu_r', heatmap_pos=[0.1, 0.1, 0.8, 0.8],
    colorbar_pos=[0.08, 0.45, 0.02, 0.15]):

    fig = plt.figure(figsize=(7, 3.5))

    h_ax = fig.add_axes(heatmap_pos)
    c_ax = fig.add_axes(colorbar_pos)

    n = M.shape[0]
    # Create the rotation matrix
    t = np.array([[1,0.5], [-1,0.5]])
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)

    # Plot the Heatmap ...
    x = A[:,1].reshape(n+1, n+1)
    y = A[:,0].reshape(n+1, n+1)
    y[y<0] = -y[y<0]

    sc = h_ax.pcolormesh(x, y, np.flipud(M), vmin=vmin, vmax=vmax, cmap=cmap,
                        edgecolor='none', snap=True, linewidth=.001)
    # colorbar
    cbar = fig.colorbar(sc, cax=c_ax, ticks=[vmin, vmax], format='%.3g')
    c_ax.tick_params(labelsize=9)
    # Hide the bottom part
    xmin = A[:,1].min()
    xmax = A[:,1].max()
    ymin = A[:,0].min()
    ymax = 0
    h_ax.fill([xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax], 'w', ec='none')
    h_ax.axis('off')

    return fig

def distance_expected(uri, tad, max_apart=5000000, correct='weight'):

    lib = cooler.Cooler(uri)
    Len = max_apart // lib.binsize
    res = lib.binsize
    cache = {i:np.array([]) for i in range(1, Len+1)}

    bychroms = tad.bychroms
    for c in bychroms.keys():
        M = lib.matrix(balance=correct).fetch(c)
        M[np.isnan(M)] = 0
        # mask gap regions
        mask = M.sum(axis=0) == 0
        n = mask.size
        idx = np.arange(n) # bin index
        tads_c = bychroms[c]
        diag_idx = np.arange(1, Len+1) # diagonal index
        for i in diag_idx:
            # observed values
            ov = M.diagonal(i)
            # build an annotation array, 1,2,... for TAD index, -1 for gap regions, and 0 for background
            iv = np.zeros(n-i, dtype=np.int16)
            row = idx[:-i]
            col = idx[i:]
            for j, t in enumerate(tads_c): # TAD index
                tmp = (row >= t[0]//res) & (col < t[1]//res)
                iv[tmp] = j + 1
            invalid = mask[:-i] | mask[i:]
            iv[invalid] = -1

            # only intra-TAD interactions
            ov = ov[iv>0]
            cache[i] = np.r_[cache[i], ov]
    
    avg = {} # expected values
    for i in cache:
        if len(cache[i]) > 20:
            avg[i] = np.mean(cache[i])
    
    return avg

def pileup(uri, DI_track, IS_track, bounds, avg, correct='weight', halfr=1000000):

    Lib = cooler.Cooler(uri)
    chromsizes = Lib.chromsizes
    res = Lib.binsize

    refavg = avg[max(avg)]
    arrs = []
    DIs = []
    ISs = []
    for b in bounds:
        c = b[0]
        cen = b[1]
        left = cen - halfr
        if left < 0:
            continue
        right = cen + halfr
        if right > chromsizes[c]:
            continue
        M = Lib.matrix(balance=correct).fetch((c, left, right))
        M[np.isnan(M)] = 0
        idx = np.arange(len(M))
        M[idx, idx] = 1 # mask the diagonal
        for i in range(1, len(M)):
            if i in avg:
                exp = avg[i]
            else:
                exp = refavg
            M[idx[:-i], idx[i:]] = M[idx[:-i], idx[i:]] / exp
            M[idx[i:], idx[:-i]] = M[idx[i:], idx[:-i]] / exp
        
        arrs.append(M)

        di = DI_track.get_interval(c, left, right)
        di[np.isnan(di)] = 0
        DIs.append(di)
        is_ = IS_track.get_interval(c, left, right)
        is_[np.isnan(is_)] = 0
        ISs.append(is_)
    
    pool = np.r_[arrs]
    M = pool.mean(axis=0)
    ISs = np.r_[ISs]
    DIs = np.r_[DIs]
    DI = DIs.mean(axis=0)
    IS = ISs.mean(axis=0)

    return M, DI, IS

def clear_frame(ax):

    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
        labelbottom=False, labeltop=False, labelleft=False, labelright=False)

# command params: tad, uri, DI, IS, resolution, outpre

tad = readHierDomain(sys.argv[1])
res = int(sys.argv[5])
uri = sys.argv[2]
di_fil, di_col = sys.argv[3].split('::')
is_fil, is_col = sys.argv[4].split('::')
di_col = int(di_col)
is_col = int(is_col)

dset = DomainSet('pseudo', tad, res=res, hier=False)
bset = BoundSet('pseudo', tad, res=res)
DI_track = sigTrack(di_fil, res=res, col=di_col, skiprows=0)
IS_track = sigTrack(is_fil, res=res, col=is_col, skiprows=1)

avg = distance_expected(uri, dset, max_apart=3000000, correct='weight')
M, DI, IS = pileup(uri, DI_track, IS_track, bset.Bounds, avg, correct='weight', halfr=600000)


outpre = sys.argv[6]
outhm = outpre + '.hm.svg'
outDI = outpre + '.DI.svg'
outIS = outpre + '.IS.svg'
fig = triangle_plot(M, vmin=0.5, vmax=1.2)
fig.savefig(outhm, dpi=200, bbox_inches='tight')
plt.close()

# plot aggregated DI track
fig = plt.figure(figsize=(7, 0.8))
di_ax = fig.add_subplot(111)
pos_mask = DI >= 0
neg_mask = DI < 0
di_ax.fill_between(np.arange(DI.size), DI, where=pos_mask, color='#FB9A99',
                    edgecolor='none')
di_ax.fill_between(np.arange(DI.size), DI, where=neg_mask, color='#A6CEE3',
                    edgecolor='none')
ymin, ymax = di_ax.get_ylim()
ax_pos = di_ax.get_position().bounds
y_ax = fig.add_axes([ax_pos[0]-0, ax_pos[1],
                    0.01, ax_pos[3]])
plot_y_axis(y_ax, ymin, ymax, size=9)
clear_frame(di_ax)
clear_frame(y_ax)
plt.savefig(outDI, bbox_inches='tight')
plt.close()

# plot aggregated IS track
fig = plt.figure(figsize=(7, 0.8))
is_ax = fig.add_subplot(111)
is_ax.plot(IS, c='#666666', lw=2.5)
ymin, ymax = is_ax.get_ylim()
ax_pos = is_ax.get_position().bounds
y_ax = fig.add_axes([ax_pos[0]-0, ax_pos[1],
                    0.01, ax_pos[3]])
plot_y_axis(y_ax, ymin, ymax, size=9)
clear_frame(is_ax)
clear_frame(y_ax)
plt.savefig(outIS, bbox_inches='tight')
plt.close()