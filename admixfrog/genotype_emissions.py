import numpy as np
from scipy.optimize import minimize, minimize_scalar
import pdb
from .distributions import gt_homo_dist
from .read_emissions2 import p_snps_given_gt
from numba import njit
from scipy.special import betainc



def _p_gt_homo(s, P, F, res=None):
    """Pr(G | Z) for homozygous hidden states

    the basic version of the homozygous genotype emission. Assumes
    infinite reference population size
    """
    n_snps = P.alpha.shape[0]
    gt = np.ones((n_snps, 3)) if res is None else res
    gt_homo_dist(a=P.alpha[:, s], b=P.beta[:, s], F=F, n_snps=n_snps, res=gt)
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


def e_tbeta(N, alpha, beta, M=1, tau=1.):
    """calculate how much the beta distribution needs to be truncated to 
    take the finite reference population size into account

    N : effective size
    M : M-th moment
    alpha, beta: number of derived/ancestral observations
    tau: fst-like population subdivision parameter
    """
    return (betainc(alpha*tau + M, beta*tau, 1 - (1 / 2 / N)) -
            betainc(alpha*tau + M, beta*tau, (1 / 2 / N))) / \
        (betainc(alpha*tau, beta*tau, 1 - (1 / 2 / N)) -
         betainc(alpha*tau, beta*tau, (1 / 2 / N))) * \
        alpha*tau / (alpha*tau + beta*tau)



def _p_gt_het_finite(a1, b1, a2, b2, N1, N2, tau1=1, tau2=1, res=None):
    """Pr(G | Z, V) for heterozygous hidden states

    emissions including a parameter Ne that measures how large the pop is...
    and a parameter tau measuring the population structure
    """
    v1, v2 = e_tbeta(N1, a1, b1, tau=tau1), e_tbeta(N2, a2, b2, tau=tau2)
    w1, w2 = 1 - v1, 1 - v2

    n_snps = len(a1)
    gt = np.empty((n_snps, 3)) if res is None else res
    gt[:, 0] = w1 * w2
    gt[:, 2] = v1 * v2
    gt[:, 1] = 1 - gt[:, 0] - gt[:, 2]
    return gt



def post_geno_py(P, cont, F, IX, error):
    """
    calculate the probability of genotype given 
    observation and hidden state
    """
    n_homo = len(F)
    n_het = int(n_homo * (n_homo - 1) / 2)
    n_snps = P.alpha.shape[0]
    n_states = n_homo + n_het
    cflat = np.array([cont[lib] for lib in P.lib])

    pg = np.zeros((n_snps, n_states, 3))

    # P(O | G)
    ll_snp = p_snps_given_gt(P, cflat, error, n_snps, IX.OBS2SNP)
    ll_snp[IX.HAPSNP, 1] = 0.0

    # P(G | Z, homo)
    for s in range(n_homo):
        _p_gt_homo(s, P, F[s], res=pg[:, s])
        pg[IX.HAPSNP, s, 0] = P.beta[IX.HAPSNP, s] / (
            P.alpha[IX.HAPSNP, s] + P.beta[IX.HAPSNP, s]
        )
        pg[IX.HAPSNP, s, 1] = 0
        pg[IX.HAPSNP, s, 2] = P.alpha[IX.HAPSNP, s] / (
            P.alpha[IX.HAPSNP, s] + P.beta[IX.HAPSNP, s]
        )
        pg[:, s] *= ll_snp

    for s1 in range(n_homo):
        for s2 in range(s1 + 1, n_homo):
            s += 1
            _p_gt_het(
                P.alpha[:, s1], P.beta[:, s1], P.alpha[:, s2], P.beta[:, s2],
                res = pg[:, s]
            )
            pg[:, s] *= ll_snp
            pg[IX.HAPSNP, s] = 0

    pg = pg / np.sum(pg, 2)[:, :, np.newaxis]
    pg[np.isnan(pg)] = 1.0 / 3  # n_states
    return np.minimum(np.maximum(pg, 0), 1)  # rounding error


def update_F(F, Z, pg, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            x = (
                np.log(_p_gt_homo(s, P, t[0]) + 1e-10)  # avoid underflow
                * Z[IX.SNP2BIN, s][:, np.newaxis]
                * pg[:, s, :]
            )
            if np.isnan(np.sum(x)):
                pdb.set_trace()
            x[IX.HAPSNP] = 0.0
            return -np.sum(x)

        prev = f([F[s]])
        OO = minimize(
            f,
            [F[s]],
            bounds=[(0, 1)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        print("[%s] \tF: [%.4f->%.4f]:\t%.4f" % (s, F[s], OO.x[0], prev - OO.fun))
        delta += abs(F[s] - OO.x[0])
        F[s] = OO.x[0]

    return delta


def geno_emissions(P, IX, F):
    """P(G | Z) for each SNP
    build table giving the probabilities of 

    """
    n_snps = P.alpha.shape[0]
    n_homo_states = P.alpha.shape[1]
    n_het_states = int(n_homo_states * (n_homo_states - 1) / 2)
    n_states = n_homo_states + n_het_states

    GT = np.zeros((n_snps, n_states, 3)) + 300
    # P(GT | Z)
    for s in range(n_homo_states):
        for g in range(3):
            _p_gt_homo(s=s, P=P, F=F[s], res=GT[:, s, :])

    s = n_homo_states
    for s1 in range(n_homo_states):
        for s2 in range(s1 + 1, n_homo_states):
            _p_gt_het(P.alpha[:, s1],
                      P.beta[:, s1],
                      P.alpha[:, s2],
                      P.beta[:, s2],
                      res=GT[:, s])
            s += 1
    assert np.allclose(np.sum(GT, 2), 1)

    GT[IX.HAPSNP, :, 1] = 0.0  # no het emissions
    GT[IX.HAPSNP, n_homo_states:] = 0.0  # no het hidden state
    for s in range(n_homo_states):
        a, b = P.alpha[IX.HAPSNP, s], P.beta[IX.HAPSNP, s]
        GT[IX.HAPSNP, s, 0] = b / (a + b)  # if a +b > 0 else 0
        GT[IX.HAPSNP, s, 2] = a / (a + b)  # if a +b > 0 else 0
    return GT

