import numpy as np
from scipy.optimize import minimize, minimize_scalar
from numba import njit
from math import lgamma, log, exp

#@njit(fastmath=True)
def get_po_given_c2(c, e, O, N, P_cont, Z, pg, rg2obs, obs2bin, obs2snp):
    ll = np.zeros(1)
    BIN = obs2bin[rg2obs]
    SNP = obs2snp[rg2obs]
    for g in range(3):
        p = c * P_cont[rg2obs] + (1. - c) * g / 2.
        p = p * (1 - e) + (1 - p) * e
        p = O[rg2obs] * np.log(p) + (N[rg2obs] - O[rg2obs]) * np.log(1 - p)
        ll += np.sum(Z[BIN] * pg[SNP, :, g] * p[:, np.newaxis])
    return ll

@njit(fastmath=True)
def get_po_given_c(c, e, O, N, P_cont, Z, pg, rg2obs, obs2bin, obs2snp):
    n_obs = len(rg2obs)
    n_states = pg.shape[1]
    ll = np.zeros(1)
    for i in range(n_obs):
        obs = rg2obs[i]
        bin_ = obs2bin[obs]
        #BIN = obs2bin[rg2obs[i]]
        #SNP = obs2snp[rg2obs[i]]
        snp = obs2snp[obs]
        for g in range(3):
            p = c * P_cont[obs] + (1. - c) * g / 2.
            p = p * (1 - e) + (1 - p) * e
            p = O[obs] * np.log(p) + (N[obs] - O[obs]) * np.log(1 - p)
            for s in range(n_states):
                ll += Z[bin_, s] * pg[snp, s, g] * p
    return ll

def update_contamination(cont, error, P, Z, pg, IX, libs):
    """
    update emissions by maximizing contamination parameter

    cont: dict of contamination rates (by library)
    Z : Pr(Z | O)
    pg: Pr(G | Z, O)



    """
    n_libs = len(libs)
    delta = 0.
    for i in range(n_libs):
        lib = libs[i]
        f_ = IX.RG2OBS[lib]
        #  assert all(lib == P.lib[f_])

        def get_po_given_c_all(cc):
            prob = get_po_given_c2(c=cc,
                                 e=error,
                                 O=P.O,
                                 N=P.N,
                                 P_cont=P.P_cont,
                                 Z=Z,
                                 pg=pg,
                                 rg2obs = IX.RG2OBS[lib],
                                 obs2snp = IX.OBS2SNP,
                                 obs2bin = IX.OBS2BIN)
            return -prob


        p0 = get_po_given_c_all(cont[lib])

        #OO =  minimize_scalar(get_po_given_c_all, bounds=(0., 1), method="Bounded")
        #print("[%s/%s]minimizing \tc: [%.4f->%.4f]:\t%.4f" % (lib, len(f_),
        #                                                           cont[lib], OO.x, p0-OO.fun))
        #cont[lib] = OO.x
        OO =  minimize(get_po_given_c_all, [cont[lib]], bounds=[(0., 1-1e-10)])
        print("[%s/%s]minimizing \tc: [%.4f->%.4f]:\t%.4f" % (lib, len(f_),
                                                                   cont[lib], OO.x[0], p0-OO.fun))
        delta += abs(cont[lib] - OO.x[0])
        cont[lib] = OO.x[0]

    return delta

@njit(fastmath=True)
def binom_pmf(O, N, p):
    res = np.power(p, O) * np.power(1.-p, N-O)
    for i, (o, n) in enumerate(zip(O, N)):
        if o >0 and n>1:
            res[i] *= exp(lgamma(n+1) - lgamma(o+1) - lgamma(n-o+1))
    return res


@njit
def p_reads_given_gt(O, N, Pcont, c, error, n_obs):
    """calculates probabilty of anc/derived reads given genotype
    """
    read_emissions = np.ones((n_obs, 3))
    for g in range(3):
        p = c * Pcont + (1 - c) * g / 2
        p = p * (1 - error) + (1 - p) * error
        # = binom.pmf(P.O, P.N, p)
        read_emissions[:, g] = binom_pmf(O, N, p)

    return read_emissions


@njit(fastmath=True)
def read2snp_emissions(read_emissions, n_snps, ix):
    snp_emissions = np.ones((n_snps, 3))
    for i, row in enumerate(ix):
        snp_emissions[row] *= read_emissions[i]
    return snp_emissions


def p_snps_given_gt(P, c, error, n_snps, ix):
    """calculates probabilty of anc/derived reads given genotype
    """
    n_obs = P.O.shape[0]
    read_emissions = p_reads_given_gt(P.O, P.N, P.P_cont, c, error, n_obs)
    return read2snp_emissions(read_emissions, n_snps, ix)
