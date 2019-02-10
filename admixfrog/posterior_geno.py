import numpy as np
from scipy.stats import binom
from scipy.optimize import minimize, minimize_scalar
from math import exp, log

try:
    from .distributions import dbetabinom, dbetabinom2
except (ModuleNotFoundError, ImportError):
    from distributions import dbetabinom, dbetabinom2

from numba import njit


def p_reads_given_gt(P, c, error):
    read_emissions = np.ones((P.O.shape[0], 3))
    for g in range(3):
        p = c * P.P_cont + (1 - c) * g / 2
        p = p * (1 - error) + (1 - p) * error
        read_emissions[:, g] = binom.pmf(P.O, P.N, p)

    return read_emissions


def _p_gt_homo_py(s, P, tau):
    n_snps = P.alpha.shape[0]
    gt = np.ones((n_snps, 3)) + 21
    dbetabinom2(P.alpha[:, s] * tau, P.beta[:, s] * tau, n_snps, gt)
    assert np.allclose(np.sum(gt, 1), 1)
    return gt


def _p_gt_het_py(s1, s2, P):
    n_snps = P.alpha.shape[0]
    gt = np.ones((n_snps, 3))
    a1, b1 = P.alpha[:, s1], P.beta[:, s1]
    a2, b2 = P.alpha[:, s2], P.beta[:, s2]
    D = (a1 + b1) * (a2 + b2)

    gt[:, 0] = b1 * b2 / D
    gt[:, 2] = a1 * a2 / D
    gt[:, 1] = 1 - gt[:, 0] - gt[:, 2]
    assert np.allclose(np.sum(gt, 1), 1)
    return gt


@njit
def _p_gt_het_py2(a1, b1, a2, b2):
    n_snps = len(a1)
    gt = np.empty((n_snps, 3))
    D = (a1 + b1) * (a2 + b2)

    gt[:, 0] = b1 * b2 / D
    gt[:, 2] = a1 * a2 / D
    gt[:, 1] = 1 - gt[:, 0] - gt[:, 2]
    return gt


def post_geno_py(P, cont, tau, IX, error):
    """
    calculate the probability of genotype given 
    observation and hidden state
    """
    n_homo = len(tau)
    n_het = int(n_homo * (n_homo - 1) / 2)
    n_snps = P.alpha.shape[0]
    n_states = n_homo + n_het
    cflat = np.array([cont[lib] for lib in P.lib])

    pg = np.zeros((n_snps, n_states, 3))

    # P(O | G)
    ll_read = p_reads_given_gt(P, cflat, error)
    ll_snp = np.ones((IX.n_snps, 3))
    for snp, obs in zip(IX.OBS2SNP, ll_read):
        ll_snp[snp] *= obs
    ll_snp[IX.HAPSNP, 1] = 0.0

    # P(G | Z, homo)
    for s in range(n_homo):
        prior_geno = _p_gt_homo_py(s, P, tau[s])
        prior_geno[IX.HAPSNP, 0] = P.beta[IX.HAPSNP, s] / (
            P.beta[IX.HAPSNP, s] + P.beta[IX.HAPSNP, s]
        )
        prior_geno[IX.HAPSNP, 1] = 0
        prior_geno[IX.HAPSNP, 2] = P.alpha[IX.HAPSNP, s] / (
            P.beta[IX.HAPSNP, s] + P.beta[IX.HAPSNP, s]
        )
        pg[:, s] = ll_snp * prior_geno
    for s1 in range(n_homo):
        for s2 in range(s1 + 1, n_homo):
            s += 1
            # prior_geno = _p_gt_het_py(s1, s2, P)
            prior_geno = _p_gt_het_py2(
                P.alpha[:, s1], P.beta[:, s1], P.alpha[:, s2], P.beta[:, s2]
            )
            pg[:, s] = ll_snp * prior_geno
            pg[IX.HAPSNP, s] = 0

    pg = pg / np.sum(pg, 2)[:, :, np.newaxis]
    pg[np.isnan(pg)] = 1.0 / 3  # n_states
    return np.minimum(np.maximum(pg, 0), 1), ll_snp, prior_geno  # rounding error


def update_tau(tau, Z, pg, P, IX):
    n_states = len(tau)
    delta = 0.
    for s in range(n_states):

        def f(t):
            x = (
                np.log(_p_gt_homo_py(s, P, t[0]))
                * Z[IX.SNP2BIN, s][:, np.newaxis]
                * pg[:, s, :]
            )
            x[IX.HAPSNP] = 0.0
            return -np.sum(x)

        prev = f([tau[s]])
        OO = minimize(
            f,
            [tau[s]],
            bounds=[(1e-4, 1e2)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        # OO =  minimize(f, [log(tau[s])])
        print("[%s] \ttau: [%.4f->%.4f]:\t%.4f" % (s, tau[s], OO.x[0], prev - OO.fun))
        delta += abs(tau[s] - OO.x[0])
        tau[s] = OO.x[0]
        # OO =  minimize_scalar(f,tol = 1e-2)
        # tau[s] = exp(OO.x)

    return delta
