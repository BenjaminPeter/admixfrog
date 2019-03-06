import numpy as np
from scipy.optimize import minimize
import pdb
from .distributions import gt_homo_dist
from .read_emissions2 import p_snps_given_gt
from numba import njit
from scipy.special import betainc
import logging


@njit(fastmath=True)
def snp2bin(e_out, e_in, ix):
    for i, row in enumerate(ix):
        e_out[row] *= e_in[i]


def scale_mat(M):
    """scale a matrix of probabilities such that it's highest value is one

    modifies M and returns log(scaling)
    """
    scaling = np.max(M, 1)[:, np.newaxis]
    M /= scaling
    assert np.allclose(np.max(M, 1), 1)
    log_scaling = np.sum(np.log(scaling))
    return log_scaling


def update_emissions(E, SNP, P, IX, est_inbreeding=False, bad_bin_cutoff=1e-150):
    """main function to calculate emission probabilities

    """
    n_homo_states = P.alpha.shape[1]

    # P(O|Z) = sum_G P(O, G | Z)
    snp_emissions = np.sum(SNP, 2)
    log_scaling = scale_mat(snp_emissions)

    E[:] = 1  # reset
    snp2bin(E, snp_emissions, IX.SNP2BIN)
    if not est_inbreeding:
        E[IX.HAPBIN, n_homo_states:] = 0

    bad_bins = np.sum(E, 1) < bad_bin_cutoff
    if sum(bad_bins) > 0:
        logging.warning("bad bins %s", sum(bad_bins))
    E[bad_bins] = bad_bin_cutoff / E.shape[1]

    log_scaling += scale_mat(E)
    return log_scaling


def _p_gt_homo(s, P, F=0, tau=1.0, res=None):
    """Pr(G | Z) for homozygous hidden states

    the basic version of the homozygous genotype emission. Assumes
    infinite reference population size
    """
    n_snps = P.alpha.shape[0]
    gt = np.ones((n_snps, 3)) if res is None else res
    gt_homo_dist(a=P.alpha[:, s], b=P.beta[:, s], F=F, tau=tau, n_snps=n_snps, res=gt)
    try:
        assert np.allclose(np.sum(gt, 1), 1)
    except AssertionError:
        pdb.set_trace()
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
    PG[np.isnan(PG)] = 0.0  # 1.0 / 3. / Z.shape[1]
    PG = np.minimum(np.maximum(PG, 0), 1)  # rounding error
    try:
        assert np.all(PG >= 0)
        assert np.all(PG <= 1)
        assert np.allclose(np.sum(PG, (1, 2)), 1)
    except AssertionError:
        pdb.set_trace()

    return PG


def update_F(F, tau, PG, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            x = np.log(_p_gt_homo(s, P, t[0], tau=tau[s]) + 1e-10) * PG[:, s, :]
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
        logging.info(
            "[%s] \tF: [%.4f->%.4f]:\t%.4f" % (s, F[s], OO.x[0], prev - OO.fun)
        )
        delta += abs(F[s] - OO.x[0])
        F[s] = OO.x[0]

    return delta


def update_Ftau(F, tau, PG, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            F, tau = t
            x = np.log(_p_gt_homo(s, P, F, tau) + 1e-10) * PG[:, s, :]
            if np.isnan(np.sum(x)):
                pdb.set_trace()
            x[IX.HAPSNP] = 0.0
            return -np.sum(x)

        prev = f([F[s], tau[s]])
        OO = minimize(
            f,
            [F[s], tau[s]],
            bounds=[(0, 1), (0, 1)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        logging.info("[%s] \tF: [%.4f->%.4f]:" % (s, F[s], OO.x[0]), end="\t")
        logging.info("T: [%.4f->%.4f]:\t%.4f" % (tau[s], OO.x[1], prev - OO.fun))
        delta += abs(F[s] - OO.x[0]) + abs(tau[s] - OO.x[1])
        F[s], tau[s] = OO.x

    return delta


def update_tau(F, tau, PG, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            x = np.log(_p_gt_homo(s, P, F=F[s], tau=t[0]) + 1e-10) * PG[:, s, :]
            if np.isnan(np.sum(x)):
                pdb.set_trace()
            x[IX.HAPSNP] = 0.0
            return -np.sum(x)

        prev = f([tau[s]])
        OO = minimize(
            f,
            [tau[s]],
            bounds=[(0, 1)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        logging.info("[%s] \tF: [%.4f->%.4f]:" % (s, F[s], F[s]), end="\t")
        logging.info("T: [%.4f->%.4f]:\t%.4f" % (tau[s], OO.x[0], prev - OO.fun))
        delta += abs(tau[s] - OO.x[0])
        tau[s] = OO.x[0]

    return delta


def update_snp_prob(SNP, P, IX, cont, error, F, tau, est_inbreeding=False):
    """
    calculate P(O, G |Z) = P(O | G) P(G | Z)
    """
    n_snps = P.alpha.shape[0]
    cflat = np.array([cont[lib] for lib in P.lib])

    # get P(G | Z)
    # save in the same array as SNP - size is the same, and
    # we do not need to allocate more memory
    update_geno_emissions(
        SNP, P, IX, F, tau, n_states=SNP.shape[1], est_inbreeding=est_inbreeding
    )

    # get P(O | G)
    ll_snp = p_snps_given_gt(P, cflat, error, n_snps, IX)
    # log_scaling = scale_mat(ll_snp)

    SNP *= ll_snp[:, np.newaxis, :]
    return 0  # log_scaling


def update_geno_emissions(GT, P, IX, F, tau, n_states, est_inbreeding):
    """P(G | Z) for each SNP
    build table giving the probabilities of P(G | Z)
    """
    n_snps, n_homo_states = P.alpha.shape

    GT[:] = 0.0
    # P(G | Z)
    for s in range(n_homo_states):
        for g in range(3):
            _p_gt_homo(s=s, P=P, F=F[s], tau=tau[s], res=GT[:, s, :])

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

    if not est_inbreeding:
        GT[IX.HAPSNP, :, 1] = 0.0  # no het emissions
        GT[IX.HAPSNP, n_homo_states:] = 0.0  # no het hidden state
        for s in range(n_homo_states):
            a, b = P.alpha[IX.HAPSNP, s], P.beta[IX.HAPSNP, s]
            GT[IX.HAPSNP, s, 0] = b / (a + b)  # if a +b > 0 else 0
            GT[IX.HAPSNP, s, 2] = a / (a + b)  # if a +b > 0 else 0

    return GT


def e_tbeta(N, alpha, beta, M=1, tau=1.0):
    """calculate how much the beta distribution needs to be truncated to
    take the finite reference population size into account

    N : effective size
    M : M-th moment
    alpha, beta: number of derived/ancestral observations
    tau: fst-like population subdivision parameter

    UNTESTED / UNUSED
    """
    return (
        (
            betainc(alpha * tau + M, beta * tau, 1 - (1 / 2 / N))
            - betainc(alpha * tau + M, beta * tau, (1 / 2 / N))
        )
        / (
            betainc(alpha * tau, beta * tau, 1 - (1 / 2 / N))
            - betainc(alpha * tau, beta * tau, (1 / 2 / N))
        )
        * alpha
        * tau
        / (alpha * tau + beta * tau)
    )


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
