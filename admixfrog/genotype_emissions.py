import numpy as np
from scipy.optimize import minimize
import pdb
from .read_emissions2 import p_snps_given_gt
from numba import njit
from math import exp, log
from .log import log_
from .utils import scale_mat, scale_mat3d
from .gllmode_emissions import update_Ftau_gllmode, _p_gt_homo, update_geno_emissions


@njit
def snp2bin(e_out, e_in, ix):
    for i, row in enumerate(ix):
        e_out[row] *= e_in[i]


@njit
def snp2bin2(e_out, e_in, ix, weight):
    """version with ld weight"""
    if np.all(weight == 1.0):
        for i, row in enumerate(ix):
            e_out[row] *= e_in[i]
    else:
        for i, row in enumerate(ix):
            e_out[row] += e_in[i] * weight[i] - weight[i]


def update_emissions(E, SNP, P, IX, bad_bin_cutoff=1e-250):
    """main function to calculate emission probabilities for each bin
    P(O | Z) = 1/S  \\ sum_G P(O, G | Z)  

    """
    #import pdb; pdb.set_trace()
    n_homo_states = P.alpha.shape[1]

    snp_emissions = np.sum(SNP, 2)

    E[:] = 1  # reset
    snp2bin(E, snp_emissions, IX.SNP2BIN)
    # snp2bin2(E, snp_emissions, IX.SNP2BIN, IX.snp_weight)
    log_.debug("mean emission %s" % np.mean(E))

    bad_bins = np.sum(E, 1) < bad_bin_cutoff
    if sum(bad_bins) > 0:
        log_.warning("%s underflow bins: %s", sum(bad_bins), np.where(bad_bins)[0])
    E[bad_bins] = bad_bin_cutoff / E.shape[1]

    log_scaling = scale_mat(E)
    return log_scaling


def update_post_geno(PG, SNP, Z, IX):
    """
    calculate P(G ,Z| O), the probability of genotype given 
    observations and hidden states Z. 

    P(G, Z| O) = P(Z|O) P(G | Z, O)
               = P(Z|O) P(O|G) P(G|Z) / sum_g P(G=g | Z)
               = P(Z|O) P(O, G|Z) / sum_g P(G=g | Z)
    Z.

    PG[n_snp x n_geno]: P( G| O'), prob of genotype given observations,
        parameters from previous iteration
    SNP[n_snp x n_states x n_geno]: P(O, G | Z) 
    Z[n_bin x n_states]: P(Z | O')

    """
    PG[:] = Z[IX.SNP2BIN, :, np.newaxis] * SNP  # P(Z|O) P(O, G | Z)
    PG /= np.sum(SNP, 2)[:, :, np.newaxis]
    PG[np.isnan(PG)] = 0.0
    PG = np.minimum(np.maximum(PG, 0), 1)  # rounding error

    try:
        assert np.all(PG >= 0)
        assert np.all(PG <= 1)
        assert np.allclose(np.sum(PG, (1, 2)), 1)
    except AssertionError:
        breakpoint()

    return PG


def update_snp_prob(
    SNP, P, IX, cont, error, F, tau, est_inbreeding=False, gt_mode=False
):
    """
    calculate P(O, G |Z) = P(O | G) P(G | Z)

    """
    #import pdb; pdb.set_trace()
    cflat = np.array([cont[lib] for lib in P.lib])
    eflat = np.array([error[lib] for lib in P.lib])

    # get P(G | Z)
    # save in the same array as SNP - size is the same, and
    # we do not need to allocate more memory
    update_geno_emissions(
        SNP, P, IX, F, tau, n_states=SNP.shape[1], est_inbreeding=est_inbreeding
    )

    assert np.allclose(np.sum(SNP, 2), 1)

    # get P(O | G)
    ll_snp = p_snps_given_gt(P, cflat, eflat, IX, gt_mode)

    SNP *= ll_snp[:, np.newaxis, :]
    log_scaling = scale_mat3d(SNP)

    return log_scaling


def update_Ftau(F, tau, PG, P, IX, est_options, gt_mode=False):
    return update_Ftau_gllmode(F, tau, PG, P, IX, est_options=est_options)
