import numpy as np
from scipy.optimize import minimize
import pdb
from .read_emissions2 import p_snps_given_gt
from numba import njit
from math import exp, log
from ..utils.log import log_
from ..utils.utils import scale_mat, scale_mat3d
from .gllmode_emissions import _p_gt_homo, update_geno_emissions


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


def update_emissions(E, SNP, P, bad_bin_cutoff=1e-250, scale_probs=True):
    """main function to calculate emission probabilities for each bin
    P(O | Z) = 1/S  \\ sum_G P(O, G | Z)

    """
    snp_emissions = np.sum(SNP, 2)

    E[:] = 1  # reset
    snp2bin(E, snp_emissions, P.SNP2BIN)
    log_.debug("mean emission %s" % np.mean(E))

    bad_bins = np.sum(E, 1) < bad_bin_cutoff
    if sum(bad_bins) > 0:
        log_.warning("%s underflow bins: %s", sum(bad_bins), np.where(bad_bins)[0])
    E[bad_bins] = bad_bin_cutoff / E.shape[1]

    log_scaling = scale_mat(E) if scale_probs else 0
    return log_scaling


def update_post_geno(X, P):
    """
    calculate P(G, Z | O), the probability of genotype given
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
    PG, SNP, Z = X.PG, X.SNP, X.Z
    PG[:] = Z[P.SNP2BIN, :, np.newaxis] * SNP  # P(Z|O) P(O, G | Z)
    PG /= np.sum(SNP, 2)[:, :, np.newaxis]
    PG[np.isnan(PG)] = 0.0
    np.clip(PG, 0, 1, out=PG)

    assert np.all(PG >= 0)
    assert np.all(PG <= 1)
    assert np.allclose(np.sum(PG, (1, 2)), 1)

    return PG


def update_snp_prob(
    SNP,
    P,
    pars,
    est_inbreeding=False,
    gt_mode=False,
    scale_probs=True,
):
    """
    calculate P(O, G |Z) = P(O | G) P(G | Z)

    """
    cflat = pars.cont[P.OBS2RG]
    eflat = pars.error[P.OBS2RG]

    # get P(G | Z)
    # save in the same array as SNP - size is the same, and
    # we do not need to allocate more memory
    update_geno_emissions(
        SNP, P, pars.F, pars.tau, n_states=SNP.shape[1], est_inbreeding=est_inbreeding
    )

    assert np.allclose(np.sum(SNP[P.diploid_snps], 2), 1)

    # get P(O | G)
    ll_snp = p_snps_given_gt(P, cflat, eflat, gt_mode)

    SNP *= ll_snp[:, np.newaxis, :]
    log_scaling = scale_mat3d(SNP) if scale_probs else 0

    return log_scaling


def update_Ftau(F, tau, PG, P, est_options):
    delta = 0.0
    for i, s in enumerate(P.states.homo_ids):

        def f(args):
            args = list(args)
            F = args.pop(0) if est_options.est_F else F[i]
            tau_ = exp(args.pop(0)) if est_options.est_tau else exp(tau[i])
            x = np.log(_p_gt_homo(s, P, F, tau_) + 1e-10) * PG[P.diploid_snps, i, :]
            if np.isnan(np.sum(x)):
                raise ValueError("nan in likelihood")
            return -np.sum(x)

        init, bounds = [], []
        if est_options.est_F:
            init.append(F[i])
            bounds.append((0, 1))
        if est_options.est_tau:
            init.append(tau[i])
            bounds.append((-10, 20))

        prev = f(init)
        OO = minimize(f, init, bounds=bounds, method="L-BFGS-B")
        opt = OO.x.tolist()

        old_F, old_tau = F[i], tau[i]
        if est_options.est_F:
            F[i] = opt.pop(0)

        if est_options.est_tau:
            tau[i] = opt.pop(0)

        log__ = "[%s] \tF: [%.4f->%.4f]\t:" % (s, old_F, F[i])
        log__ += "T: [%.4f->%.4f]:\t%.4f" % (old_tau, tau[i], prev - OO.fun)
        log_.info(log__)
        delta += abs(F[i] - old_F) + abs(tau[i] - old_tau)
    return delta
