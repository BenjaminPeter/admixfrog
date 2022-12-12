import numpy as np
from scipy.optimize import minimize
from numba import njit
from math import lgamma, exp


@njit
def binom_pmf(O, N, p):
    res = np.power(p, O) * np.power(1.0 - p, N - O)
    for i, (o, n) in enumerate(zip(O, N)):
        if o > 0 and n > 1:
            res[i] *= exp(lgamma(n + 1) - lgamma(o + 1) - lgamma(n - o + 1))
    return res


@njit
def p_reads_given_gt_gllmode(O, N, Pcont, c, error, n_obs):
    """calculates P(O | G); probabilty of anc/derived reads given genotype
    per read group

    """
    read_emissions = np.ones((n_obs, 3), float)
    for g in range(3):
        p = c * Pcont + (1 - c) * g / 2
        p = p * (1 - error) + (1 - p) * error
        read_emissions[:, g] = binom_pmf(O, N, p)

    return read_emissions


def p_reads_given_gt(*args, gt_mode=False, **kwargs):
    if gt_mode:
        return p_reads_given_gt_gtmode(*args, **kwargs)
    else:
        return p_reads_given_gt_gllmode(*args, **kwargs)


@njit
def p_reads_given_gt_gtmode(O, N, Pcont, c, error, n_obs):
    """calculates P(O | G); probabilty of anc/derived genotype given input genotype"""
    n_gt = 3
    read_emissions = np.ones((n_obs, n_gt))
    for g in range(n_gt):
        read_emissions[O / N == g / 2, g] = 1 - 2 * error[O / N == g / 2]
        read_emissions[O / N != g / 2, g] = error[O / N != g / 2]

    return read_emissions


@njit
def read2snp_emissions(read_emissions, n_snps, ix):
    n_gt = 3
    snp_emissions = np.ones((n_snps, n_gt))
    for i, row in enumerate(ix):
        snp_emissions[row] *= read_emissions[i]
    return snp_emissions


def p_snps_given_gt(P, c, error, gt_mode=False):
    """calculates probabilty of observed read data given genotype"""
    read_emissions = p_reads_given_gt(
        P.O, P.N, P.psi, c, error, P.n_obs, gt_mode=gt_mode
    )
    return read2snp_emissions(read_emissions, P.n_snps, P.OBS2SNP)
