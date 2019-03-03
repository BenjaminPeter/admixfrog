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

    UNTESTED / UNUSED
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

    UNTESTED / UNUSED
    """
    v1, v2 = e_tbeta(N1, a1, b1, tau=tau1), e_tbeta(N2, a2, b2, tau=tau2)
    w1, w2 = 1 - v1, 1 - v2

    n_snps = len(a1)
    gt = np.empty((n_snps, 3)) if res is None else res
    gt[:, 0] = w1 * w2
    gt[:, 2] = v1 * v2
    gt[:, 1] = 1 - gt[:, 0] - gt[:, 2]
    return gt


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
    Z[n_snp x n_states]: P(Z | O')
    """
    PG[:] = Z[IX.SNP2BIN, :, np.newaxis] * SNP
    PG /= np.sum(SNP, 2)[:, :, np.newaxis]
    PG[np.isnan(PG)] = 0 #1.0 / 3. / Z.shape[1]

    assert np.all(PG >= 0)
    assert np.all(PG <= 1)
    assert np.allclose(np.sum(PG, (1, 2)), 1)

    return np.minimum(np.maximum(PG, 0), 1)  # rounding error


def update_F(F, Z, PG, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            x = (
                np.log(_p_gt_homo(s, P, t[0]) + 1e-10)  # avoid underflow
                * PG[:, s, :]
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


def get_snp_prob(SNP, P, IX, cont, error, F):
    """
    calculate P(O, G |Z) = P(O | G) P(G | Z)
    """
    n_homo = len(F)
    n_snps = P.alpha.shape[0]
    cflat = np.array([cont[lib] for lib in P.lib])

    # get P(O | G)
    # save in the same array as SNP - size is the same, and
    # we do not need to allocate more memory
    update_geno_emissions(SNP, P, IX, F, n_states=SNP.shape[1])

    # get P(G | Z)
    ll_snp = p_snps_given_gt(P, cflat, error, n_snps, IX)

    SNP *= ll_snp[:, np.newaxis, :]

    # haploidize x chrom
    s = n_homo
    for s1 in range(n_homo):
        for s2 in range(s1 + 1, n_homo):
            SNP[IX.HAPSNP, s] = 0
            s += 1


def update_geno_emissions(GT, P, IX, F, n_states):
    """P(G | Z) for each SNP
    build table giving the probabilities of P(G | Z)
    """
    n_snps, n_homo_states = P.alpha.shape

    #GT = np.zeros((n_snps, n_states, 3)) + 300
    # P(G | Z)
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

    GT[IX.HAPSNP, :, 1] = 0.0            # no het emissions
    GT[IX.HAPSNP, n_homo_states:] = 0.0  # no het hidden state
    for s in range(n_homo_states):
        a, b = P.alpha[IX.HAPSNP, s], P.beta[IX.HAPSNP, s]
        GT[IX.HAPSNP, s, 0] = b / (a + b)  # if a +b > 0 else 0
        GT[IX.HAPSNP, s, 2] = a / (a + b)  # if a +b > 0 else 0
    return GT

