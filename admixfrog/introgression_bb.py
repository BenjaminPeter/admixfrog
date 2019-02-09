import numpy as np
import pandas as pd
import sys
from collections import namedtuple, defaultdict
from scipy.stats import binom
from scipy.optimize import minimize
from .distributions import dbetabinom

#np.set_printoptions(suppress=True, precision=4)


def p_reads_given_gt(P, n_obs, c, error):
    read_emissions = np.ones((n_obs, 3)) 
    for g in range(3):
        p = c * P.P_cont + (1 - c) * g / 2
        p = p * (1-error) + (1-p) * error
        read_emissions[:,g] = binom.pmf(P.O, P.N, p)

    return read_emissions

def get_po_given_c_py(
    c, 
    e,
    O,
    N,
    P_cont,
    G,
    pg,
    IX
    ) :

    ll = 0.
    n_states = pg.shape[1]
    for i in range(IX.n_obs):
        chrom, bin_ = IX.OBS2BIN[i]
        snp = IX.OBS2SNP[i]
        for s in range(n_states):
            for g in range(3):
                p = c * P_cont[i] + (1.-c) * g / 2.
                p = p * (1-e) + (1-p) * e
                ll += G[chrom][bin_, s] * pg[snp, s, g] * (  O[i] * log(p) + (N[i] - O[i]) * log(1-p))
    return ll

def get_emissions_bb_py(P,
                        IX,
                        cont,
                        tau,
                        error, bad_bin_cutoff=1e-100
):
    n_homo_states = P.alpha.shape[1]
    n_het_states = int(n_homo_states * (n_homo_states - 1) / 2)
    n_states = n_homo_states + n_het_states 
    n_obs = P.O.shape[0]
    n_snps = P.alpha.shape[0]
    c = np.array([cont[l] for l in P.lib])

    read_emissions = p_reads_given_gt(P, n_obs, c, error)
    #P(SNP | GT)
    gt_emissions = np.ones((n_snps, 3)) 
    for i, row in enumerate(IX.OBS2SNP):
        gt_emissions[row] *= read_emissions[i]


    GT = np.zeros((n_snps, n_states, 3)) + 300
    #P(GT | Z)
    for s in range(n_homo_states):
        for g in range(3):
            GT[:,s, g] = dbetabinom(np.zeros(n_snps, int) + g,
                            np.zeros(n_snps, int) + 2,
                            tau[s] * (P.alpha[:, s] ),
                            tau[s] * (P.beta[:, s] ))

    s = n_homo_states
    fa, fb = (P.alpha ).T, (P.beta ).T
    pp = fa / (fa + fb)
    qq = 1 - pp
    for s1 in range(n_homo_states):
        for s2 in range(s1+1, n_homo_states):
            for g in range(3):
                GT[:, s, 0] = qq[s1] * qq[s2]
                GT[:, s, 2] = pp[s1] * pp[s2]
                GT[:, s, 1] = 1 - GT[:, s, 0] - GT[:, s, 2]
            s += 1
    assert np.allclose( np.sum(GT, 2), 1)

    snp_emissions = np.sum(GT * gt_emissions[:, np.newaxis, :],2)   

    bin_emissions = [np.ones((IX.bin_sizes[i], n_states), dtype=np.float64) for i in range(IX.n_chroms)]
    for (chrom, bin_), snp in zip(IX.SNP2BIN, snp_emissions):
        print(snp)
        bin_emissions[chrom][bin_] *= snp

    return bin_emissions


    
