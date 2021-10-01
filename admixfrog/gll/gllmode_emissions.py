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


def update_Ftau_gllmode(F, tau, PG, P, IX, est_options):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(args):
            args = list(args)
            F = args.pop(0) if est_options["est_F"] else F[s]
            tau = exp(args.pop(0)) if est_options["est_tau"] else exp(tau[s])
            x = np.log(_p_gt_homo(s, P, F, tau) + 1e-10) * PG[IX.diploid_snps, s, :]
            if np.isnan(np.sum(x)):
                raise ValueError("nan in likelihood")
            return -np.sum(x)

        init, bounds = [], []
        if est_options["est_F"]:
            init.append(F[s])
            bounds.append((0, 1))
        if est_options["est_tau"]:
            init.append(tau[s])
            bounds.append((-10, 20))

        prev = f(init)
        OO = minimize(f, init, bounds=bounds, method="L-BFGS-B")
        opt = OO.x.tolist()

        old_F, old_tau = F[s], tau[s]
        if est_options["est_F"]:
            F[s] = opt.pop(0)

        if est_options["est_tau"]:
            tau[s] = opt.pop(0)

        log__ = "[%s] \tF: [%.4f->%.4f]\t:" % (s, old_F, F[s])
        log__ += "T: [%.4f->%.4f]:\t%.4f" % (old_tau, tau[s], prev - OO.fun)
        log_.info(log__)
        delta += abs(F[s] - old_F) + abs(tau[s] - old_tau)
    return delta


def update_geno_emissions_diploid(GT, P, IX, F, tau, n_states, est_inbreeding):
    """P(G | Z) for diploid SNP. 
    build table giving the probabilities of P(G | Z)
    """
    n_snps, n_homo_states = P.alpha.shape

    GT[:] = 0.0
    # P(G | Z)
    for s in range(n_homo_states):
        _p_gt_homo(s=s, P=P, F=F[s], tau=exp(tau[s]), res=GT[:, s, :])

    for s1 in range(n_homo_states):
        for s2 in range(s1 + 1, n_homo_states):
            s += 1
            _p_gt_het(
                P.alpha[:, s1],
                P.beta[:, s1],
                P.alpha[:, s2],
                P.beta[:, s2],
                res=GT[:, s, :],
            )

    if est_inbreeding:
        for i in range(n_homo_states):
            _p_gt_hap(P.alpha[:, i], P.beta[:, i], res=GT[:, s + i + 1, :])

    try:
        assert np.allclose(np.sum(GT, 2), 1)
    except AssertionError:
        pdb.set_trace()

    return GT


def update_geno_emissions_haploid(GT, P):
    """P(G | Z) for each SNP for haploid SNP (e.g. X-chromosome) only
    """
    n_snps, n_homo = P.alpha.shape

    GT[:] = 0.0
    # P(G | Z)
    for i in range(n_homo):
        _p_gt_hap(P.alpha_hap[:, i], P.beta_hap[:, i], res=GT[:, i, :])

    assert np.allclose(np.sum(GT[:, :n_homo], 2), 1)
    return GT


def update_geno_emissions(GT, P, IX, *args, **kwargs):
    update_geno_emissions_diploid(GT[IX.diploid_snps], P, IX, *args, **kwargs)
    update_geno_emissions_haploid(GT[IX.haploid_snps], P)


def update_geno_emissions_diploid_sfs(GT, P, IX, F, tau, n_states, est_inbreeding):
    """P(G | tau, F) for diploid SNP. 
    compared to admixfrog model, we do not have a latent state Z
    build table giving the probabilities of P(G | Z)
    """
    n_snps, n_homo_states = P.alpha.shape

    GT[:] = 0.0
    # P(G | Z)
    for s in range(n_homo_states):
        _p_gt_homo_sfs(s=s, P=P, F=F[s], tau=exp(tau[s]), res=GT[:, s, :])

    for s1 in range(n_homo_states):
        for s2 in range(s1 + 1, n_homo_states):
            s += 1
            _p_gt_het(
                P.alpha[:, s1],
                P.beta[:, s1],
                P.alpha[:, s2],
                P.beta[:, s2],
                res=GT[:, s, :],
            )

    if est_inbreeding:
        for i in range(n_homo_states):
            _p_gt_hap(P.alpha[:, i], P.beta[:, i], res=GT[:, s + i + 1, :])

    try:
        assert np.allclose(np.sum(GT, 2), 1)
    except AssertionError:
        pdb.set_trace()

    return GT
