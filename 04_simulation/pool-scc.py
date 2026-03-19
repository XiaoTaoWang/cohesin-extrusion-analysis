from typing import Union
from scipy.sparse import SparseEfficiencyWarning
from scipy.interpolate import make_interp_spline
from scipy.interpolate import UnivariateSpline
import warnings, cooler, glob, os, matplotlib, joblib
import numpy as np
import scipy.sparse as sp
import scipy.stats as ss
import matplotlib.pyplot as plt

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def smooth(signal: sp.coo_matrix, h: int) -> sp.coo_matrix:
    """
    Smooth (blur) input sparse Hi-C sparse matrix using a uniform kernel of
    width 2*h+1.

    Parameters
    ----------
    signal : scipy.sparse.coo_matrix
        Input intrachromosomal Hi-C matrix to smooth (upper triangle matrix).
    h : int
        Half width of the kernel (neighbourhood size).

    Returns
    -------
    scipy.sparse.coo_matrix :
        The smoothed matrix.
    """

    ker_m = 2 * h + 1
    kernel = np.ones((ker_m, ker_m)) / (ker_m**2)
    sig_m, sig_n = signal.shape

    # Sanity checks
    if sp.issparse(kernel):
        raise ValueError("cannot handle kernel in sparse format")
    if not sp.issparse(signal):
        raise ValueError("cannot handle signal in dense format")

    constant_kernel = kernel[0, 0]

    # Simplified convolution for the special case where kernel is constant:
    l_subkernel_sp = sp.diags(
        np.ones(ker_m), np.arange(ker_m), shape=(sig_m - ker_m + 1, sig_m), format="csr"
    )
    r_subkernel_sp = sp.diags(
        np.ones(ker_m),
        -np.arange(ker_m),
        shape=(sig_n, sig_n - ker_m + 1),
        format="csr",
    )
    out = (l_subkernel_sp @ signal) @ r_subkernel_sp
    out *= constant_kernel
    # Resize matrix: increment rows and cols by half kernel and set shape to input
    # matrix, effectively adding margins.
    out = out.tocoo()
    rows, cols = out.row + h, out.col + h
    out = sp.coo_matrix((out.data, (rows, cols)), shape=(sig_m, sig_n), dtype=np.float64)

    return out


def vstrans(d1: "np.ndarray[float]", d2: "np.ndarray[float]") -> "np.ndarray[float]":
    """
    Variance stabilizing transformation to normalize read counts before
    computing stratum correlation. This normalizes counts so that different
    strata share similar dynamic ranges.

    Parameters
    ----------
    d1 : numpy.ndarray of floats
        Diagonal of the first matrix.
    d2 : numpy.ndarray of floats
        Diagonal of the second matrix.

    Returns
    -------
    r2k : numpy.ndarray of floats
        Array of weights to use to normalize counts.
    """
    # Get ranks of counts in diagonal
    ranks_1 = np.argsort(d1) + 1
    ranks_2 = np.argsort(d2) + 1
    # Scale ranks betweeen 0 and 1
    nranks_1 = ranks_1 / max(ranks_1)
    nranks_2 = ranks_2 / max(ranks_2)
    tot_ranks = len(ranks_1)
    r2k = np.sqrt(np.var(nranks_1 / tot_ranks) * np.var(nranks_2 / tot_ranks))
    return r2k

def diag_trim(
    mat: Union[sp.dia_matrix, np.ndarray], n: int
) -> Union[sp.dia_matrix, np.ndarray]:
    """
    Trim an upper triangle sparse matrix so that only the first n diagonals are
    kept.

    Parameters
    ----------

    mat : scipy.sparse.dia_matrix or numpy.ndarray
        The sparse matrix to be trimmed
    n : int
        The number of diagonals from the center to keep (0-based).

    Returns
    -------
    scipy.sparse.dia_matrix or numpy.ndarray:
        The diagonally trimmed upper triangle matrix with only the first n
        diagonal.
    """
    if not sp.issparse(mat):
        trimmed = mat.copy()
        n_diags = trimmed.shape[0]
        for diag in range(n, n_diags):
            set_mat_diag(trimmed, diag, 0)
        return trimmed

    if mat.format != "dia":
        raise ValueError("input type must be scipy.sparse.dia_matrix")
    # Create a new matrix from the diagonals below max dist (faster than removing them)
    keep_offsets = np.flatnonzero((mat.offsets <= n) & (mat.offsets >= 0))

    trimmed = sp.dia_matrix(
        (mat.data[keep_offsets], mat.offsets[keep_offsets]), shape=mat.shape
    )

    return trimmed


def set_mat_diag(mat: np.ndarray, diag: int = 0, val: int = 0) -> float:
    """
    Set the nth diagonal of a symmetric 2D numpy array to a fixed value.
    Operates in place.

    Parameters
    ----------
    mat : numpy.ndarray
        Symmetric 2D array of floats.
    diag : int
        0-based index of the diagonal to modify. Use negative values for the
        lower half.
    val : float
        Value to use for filling the diagonal
    """
    rows = mat.shape[0]
    step = rows + 1
    start = diag
    end = rows**2 - diag * rows
    mat.flat[start:end:step] = val

def get_scc(
    mat1: 'scipy.sparse.csr_matrix', mat2: 'scipy.sparse.csr_matrix', max_bins: int
) -> float:
    """
    Compute the stratum-adjusted correlation coefficient (SCC) between two
    Hi-C matrices up to max_dist. A Pearson correlation coefficient is computed
    for each diagonal in the range of 0 to max_dist and a weighted sum of those
    coefficients is returned.
    Parameters
    ----------
    mat1 : scipy.sparse.csr_matrix
        First matrix to compare.
    mat2 : scipy.sparse.csr_matrix
        Second matrix to compare.
    max_bins : int
        Maximum distance at which to consider, in bins.
    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    corr_diag = np.zeros(len(range(max_bins)))
    weight_diag = corr_diag.copy()
    for d in range(max_bins):
        d1 = mat1.diagonal(d)
        d2 = mat2.diagonal(d)
        # Silence NaN warnings: this happens for empty diagonals and will
        # not be used in the end.
        with warnings.catch_warnings():
            # Warning does not exist in older scipy versions (<=1.2)
            try:
                warnings.filterwarnings(
                    "ignore", category=ss.PearsonRConstantInputWarning
                )
            except AttributeError:
                pass
            # Compute raw pearson coeff for this diag
            corr_diag[d] = ss.pearsonr(d1, d2)[0]
        # Compute weight for this diag
        r2k = vstrans(d1, d2)
        weight_diag[d] = len(d1) * r2k
    # Normalize weights
    weight_diag /= sum(weight_diag)

    # Weighted sum of coefficients to get SCCs
    scc = np.sum(corr_diag * weight_diag)

    return scc

def array2sparse(arr):

    from scipy.sparse import coo_matrix

    x, y = np.nonzero(arr)
    value = arr[x, y]
    shape = arr.shape

    M = coo_matrix((value, (x, y)), shape=shape)

    return M

def calculate_scc(ref, uri, max_dist=300000, res=2000, h=3):
    

    max_bins = max_dist // res
    clr = cooler.Cooler(uri)
    query = clr.matrix(balance=False).fetch('chr_sim')

    ref = array2sparse(ref)
    query = array2sparse(query)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
        ref = diag_trim(ref.todia(), max_bins + h).tocoo()
        query = diag_trim(query.todia(), max_bins + h).tocoo()
    smooth_1 = smooth(ref, h)
    smooth_2 = smooth(query, h)

    scc = get_scc(smooth_1, smooth_2, max_bins=max_bins)
    
    return scc

def parse_folder(path):

    uris = glob.glob(os.path.join(path, '*.cool'))
    D = {}
    for uri in uris:
        name = os.path.split(uri)[-1].split('.')[0]
        t = int(name.split('_')[-1][1:])
        D[t] = uri
    
    return D

if __name__ == "__main__":

    '''
    ref = np.load('experi.npy')
    queue1 = parse_folder('one-sided-stability0')
    scc_dict1 = {}
    for t, uri in queue1.items():
        scc_dict1[t] = calculate_scc(ref, uri)
    
    time_ = sorted(scc_dict1.keys())
    one_sided = [scc_dict1[t] for t in time_]

    queue2 = parse_folder('two-sided-speed0.5')
    scc_dict2 = {}
    for t, uri in queue2.items():
        scc_dict2[t] = calculate_scc(ref, uri)
    
    time_ = sorted(scc_dict2.keys())
    two_sided = [scc_dict2[t] for t in time_]
    joblib.dump((time_, one_sided, two_sided), 'scc.300kb.h3.pkl')
    '''
    time_, one_sided, two_sided = joblib.load('scc.1mb.h3.pkl')
    fig = plt.figure(figsize=(3.5, 1))
    ax = fig.add_subplot(111)
    x_new = np.linspace(time_[0], time_[-1], 100)
    spline = make_interp_spline(time_, one_sided, k=3)
    l1, = ax.plot(x_new, spline(x_new), color='#7FC97F', linewidth=1)
    spline = make_interp_spline(time_, two_sided, k=3)
    l2, = ax.plot(x_new, spline(x_new), color='#FDC086', linewidth=1, linestyle=':')
    ax.set_xlabel('Simulation Time', fontsize=7)
    ax.set_ylabel('SCC', fontsize=7)
    #ax.legend([l1, l2], ['One-sided', 'Two-sided'], fontsize=6)
    ax.xaxis.set_tick_params(width=0.8, labelsize=6)
    ax.yaxis.set_tick_params(width=0.8, labelsize=6)
    #xticks = list(range(5, 131, 25))
    #ax.set_xticks(xticks)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['left'].set_linewidth(0.8)
    plt.savefig('scc.1mb.h3.pdf', dpi=500, bbox_inches='tight')
    plt.close()