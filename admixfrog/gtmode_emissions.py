from scipy.optimize import minimize
import numpy as np
from math import exp, log
from .log import log_
from numba import njit
from .distributions import gt_homo_gtmode
from .utils import scale_mat

def _p_gt_homo_gtmode(s, P, F=0, tau=1.0, res=None):
    """Pr(G | Z) for homozygous hidden states

    this version assumes that genotypes are known, the number of alt alleles
    is give in P.O
    """
    n_snps = P.alpha.shape[0]
    gt = np.ones(n_snps) if res is None else res

    gt_homo_gtmode(o=P.O, n = P.N, 
                   a=P.alpha[:, s], b=P.beta[:, s], F=F, tau=tau, n_snps=n_snps, res=gt)
    return np.minimum(np.maximum(gt, 0), 1)  # rounding error

@njit
def _p_gt_het_gtmode(o, n, a1, b1, a2, b2, res=None):
    """Pr(G | Z) for heterozygous hidden states

    assumes genotypes are known, o is the number of derived alleles
    """
    n_snps = len(a1)
    gt = np.empty(n_snps) if res is None else res

    for i in range(n_snps):
        if n[i] == 0:
            gt[i] = 1
        elif n[i] == 2:
            if o[i] == 0:
                gt[i] = b1[i] * b2[i] / (a1[i] + b1[i]) * (a2[i] + b2[i])
            elif o[i] == 2:
                gt[i] = a1[i] * a2[i] / (a1[i] + b1[i]) * (a2[i] + b2[i])
            elif o[i] == 1:
                gt[i] = (a1[i] * b2[i] + a2[i] * b1[i]) / (a1[i] + b1[i]) * (a2[i] + b2[i])
            else:
                raise ValueError("bad gt")
                gt[i] = -1
        else:
            raise ValueError("bad gt")
            gt[i] = -1

    return gt

def update_Ftau_gtmode(F, tau, Z, P, IX):
    n_states = len(F)
    delta = 0.0
    for s in range(n_states):

        def f(t):
            F, tau = t
            x = np.log(_p_gt_homo_gtmode(s, P, F, exp(tau)) + 1e-10) * Z[IX.SNP2BIN, s]
            if np.isnan(np.sum(x)):
                pdb.set_trace()
            #x[IX.HAPSNP] = 0.0
            return -np.sum(x)

        prev = f([F[s], tau[s]])
        OO = minimize(
            f,
            [F[s], tau[s]],
            bounds=[(0, 1), (-10, 10)],
            method="L-BFGS-B",
            options=dict([("gtol", 1e-2)]),
        )
        log__ = "[%s] \tF: [%.4f->%.4f]\t:" % (s, F[s], OO.x[0])
        log__ += "T: [%.4f->%.4f]:\t%.4f" % (tau[s], OO.x[1], prev - OO.fun)
        log_.info(log__)
        delta += abs(F[s] - OO.x[0]) + abs(tau[s] - OO.x[1])
        F[s], tau[s] = OO.x

    return delta


def update_geno_emissions_gt(GT, P, IX, F, tau, est_inbreeding):
    """P(G | Z) for each SNP
    build table giving the probabilities of P(G | Z)

    in gt mode, P.alpha contins the observed number of derived alleles
    """
    n_snps, n_homo_states = P.alpha.shape

    GT[:] = 0.0
    # P(G | Z)
    for s in range(n_homo_states):
        _p_gt_homo_gtmode(s=s, P=P, F=F[s], tau=exp(tau[s]), res=GT[:, s])

    for s1 in range(n_homo_states):
        for s2 in range(s1 + 1, n_homo_states):
            s += 1
            _p_gt_het_gtmode(
                P.O,
                P.N,
                P.alpha[:, s1],
                P.beta[:, s1],
                P.alpha[:, s2],
                P.beta[:, s2],
                res=GT[:, s],
            )

    if est_inbreeding:
        raise NotImplemented()
        for i in range(n_homo_states):
            _p_gt_hap_gtmode(P.alpha[:, i], P.beta[:, i], res=GT[:, s + i + 1, :])

    log_scaling = scale_mat(GT)

    return log_scaling
