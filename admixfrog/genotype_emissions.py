import numpy as np
from scipy.optimize import minimize
import pdb
from .read_emissions2 import p_snps_given_gt
from numba import njit
from math import exp, log
from .log import log_
from .utils import scale_mat, scale_mat3d
from .gtmode_emissions import update_Ftau_gtmode
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
    P(O | Z) = 1/S  \sum_G P(O, G | Z)  

    """
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
        pdb.set_trace()

    # PG *= IX.snp_weight[:, np.newaxis, np.newaxis]

    return PG


def update_snp_prob(
    SNP, P, IX, cont, error, F, tau, est_inbreeding=False, gt_mode=False
):
    """
    calculate P(O, G |Z) = P(O | G) P(G | Z)

    """
    cflat = np.array([cont[lib] for lib in P.lib])

    # get P(G | Z)
    # save in the same array as SNP - size is the same, and
    # we do not need to allocate more memory
    update_geno_emissions(
        SNP, P, IX, F, tau, n_states=SNP.shape[1], est_inbreeding=est_inbreeding
    )

    # get P(O | G)
    ll_snp = p_snps_given_gt(P, cflat, error, IX, gt_mode)

    SNP *= ll_snp[:, np.newaxis, :]
    log_scaling = scale_mat3d(SNP)

    return log_scaling


def update_F(F, tau, PG, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            x = np.log(_p_gt_homo(s, P, t[0], tau=exp(tau[s])) + 1e-10) * PG[:, s, :]
            if np.isnan(np.sum(x)):
                pdb.set_trace()
            return -np.sum(x)

        prev = f([F[s]])
        OO = minimize(
            f,
            [F[s]],
            bounds=[(0, 1)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        log_.info("[%s] \tF: [%.4f->%.4f]:\t%.4f" % (s, F[s], OO.x[0], prev - OO.fun))
        delta += abs(F[s] - OO.x[0])
        F[s] = OO.x[0]

    return delta


def update_Ftau(F, tau, PG, P, IX, gt_mode=False):
    if gt_mode:
        return update_Ftau_gtmode(F, tau, Z=PG, P=P, IX=IX)
    else:
        return update_Ftau_gllmode(F, tau, PG, P, IX)


def update_tau(F, tau, PG, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            x = np.log(_p_gt_homo(s, P, F=F[s], tau=exp(t[0])) + 1e-10) * PG[:, s, :]
            if np.isnan(np.sum(x)):
                pdb.set_trace()
            return -np.sum(x)

        prev = f([tau[s]])
        OO = minimize(
            f,
            [tau[s]],
            bounds=[(0, 10)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        log__ = "[%s] \tF: [%.4f->%.4f]\t:" % (s, F[s], F[s])
        log__ += "T: [%.4f->%.4f]:\t%.4f" % (tau[s], OO.x[0], prev - OO.fun)
        log_.info(log__)
        delta += abs(tau[s] - OO.x[0])
        tau[s] = OO.x[0]

    return delta
