import pandas as pd
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln
from scipy.stats import binom
from numba import njit
import pdb


@njit(fastmath=True)
def get_po_given_c(c, e, O, N, P_cont, Z, pg, rg2obs, obs2bin, obs2snp):
    n_obs = len(rg2obs)
    n_states = pg.shape[1]
    ll = 0.
    for i in range(n_obs):
        obs = rg2obs[i]
        bin_ = obs2bin[obs]
        snp = obs2snp[obs]
        for g in range(3):
            p = c * P_cont[obs] + (1. - c) * g / 2.
            p = p * (1 - e) + (1 - p) * e
            p = O[obs] * np.log(p) + (N[obs] - O[obs]) * np.log(1 - p)
            for s in range(n_states):
                print(Z.shape, pg.shape)
                print(Z[bin_, s])
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
        assert all(lib == P.lib[f_])

        def get_po_given_c_all(cc):
            prob = get_po_given_c(c=cc,
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

def p_reads_given_gt(P, c, error):
    """calculates probabilty of anc/derived reads given genotype
    """
    read_emissions = np.ones((P.O.shape[0], 3))
    for g in range(3):
        p = c * P.P_cont + (1 - c) * g / 2
        p = p * (1 - error) + (1 - p) * error
        read_emissions[:, g] = binom.pmf(P.O, P.N, p)

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
    read_emissions = p_reads_given_gt(P, c, error)
    return read2snp_emissions(read_emissions, n_snps, ix)
