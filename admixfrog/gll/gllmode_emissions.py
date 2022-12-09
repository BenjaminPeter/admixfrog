from scipy.optimize import minimize
from ..utils.log import log_
import numpy as np
import pdb
from ..utils.distributions import gt_homo_dist
from numba import njit
from math import exp, log


def _p_gt_homo(s, P, F=0, tau=1.0, res=None):
    """Pr(G | Z) for homozygous hidden states

    the basic version of the homozygous genotype emission. Assumes
    infinite reference population size
    """
    n_snps = P.alpha.shape[0]
    gt = np.ones((n_snps, 3)) if res is None else res
    gt_homo_dist(a=P.alpha[:, s], b=P.beta[:, s], F=F, tau=tau, n_snps=n_snps, res=gt)

    assert np.allclose(np.sum(gt, 1), 1)
    return np.minimum(np.maximum(gt, 0), 1)  # rounding error


@njit
def _p_gt_het(a1, b1, a2, b2, res=None):
    """Pr(G | Z) for heterozygous hidden states

    the basic version of the heterozygous genotype emission. Assumes
    infinite reference population size
    """
    n_snps = len(a1)
    gt = np.empty((n_snps, 3)) if res is None else res
    D = (a1 + b1) * (a2 + b2)

    gt[:, 0] = b1 * b2 / D
    gt[:, 2] = a1 * a2 / D
    gt[:, 1] = 1 - gt[:, 0] - gt[:, 2]
    return gt


@njit
def _p_gt_hap(a1, b1, res=None):
    """Pr(G | Z) for haploid hidden states

    the basic version of the haploid genotype emission. Assumes
    infinite reference population size
    """
    n_snps = len(a1)
    gt = np.empty((n_snps, 2)) if res is None else res

    gt[:, 0] = b1 / (a1 + b1)
    gt[:, 1] = 0.0
    gt[:, 2] = a1 / (a1 + b1)
    return gt


def update_geno_emissions_diploid(GT, P, IX, F, tau, n_states, est_inbreeding):
    """P(G | Z) for diploid SNP.
    build table giving the probabilities of P(G | Z)
    """

    n_snps = P.alpha.shape[0]

    GT[:] = 0.0
    # P(G | Z)
    for i, s in enumerate(P.S.homo_ids):
        _p_gt_homo(s=s, P=P, F=F[i], tau=exp(tau[i]), res=GT[:, i, :])

    for i, (s1, s2) in enumerate(P.S.het_ids):
        _p_gt_het(
            P.alpha[:, s1],
            P.beta[:, s1],
            P.alpha[:, s2],
            P.beta[:, s2],
            res=GT[:, i + P.S.n_homo, :],
        )

    if est_inbreeding:
        n_non_inbred = P.S.n_homo + P.S.n_het
        for i, s in enumerate(P.S.homo_ids):
            _p_gt_hap(P.alpha[:, s], P.beta[:, s], res=GT[:, i + n_non_inbred, :])

    assert np.allclose(np.sum(GT, 2), 1)

    return GT


def update_geno_emissions_haploid(GT, P):
    """P(G | Z) for each SNP for haploid SNP (e.g. X-chromosome) only"""
    n_snps = P.alpha.shape[0]

    GT[:] = 0.0
    # P(G | Z)
    for i in range(P.S.n_hap):
        _p_gt_hap(P.alpha_hap[:, i], P.beta_hap[:, i], res=GT[:, i, :])

    assert np.allclose(np.sum(GT[:, : P.S.n_hap], 2), 1)
    return GT


def update_geno_emissions(GT, P, IX, *args, **kwargs):
    update_geno_emissions_diploid(GT[IX.diploid_snps], P, IX, *args, **kwargs)
    update_geno_emissions_haploid(GT[IX.haploid_snps], P)
