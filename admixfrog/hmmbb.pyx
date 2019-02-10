#!python
#cython: language_level=3
#cython: infer_types=True

import pandas as pd
import numpy as np
cimport numpy as np
cimport scipy.special.cython_special as scs
cimport cython
from libc.math cimport pow, log, exp
from libc.stdio cimport printf
from scipy.optimize import minimize, minimize_scalar
from .distributions cimport *


def get_emissions_bb_cy(P,
                        IX,
                        cont,
                        tau,
                        error=2e-2, bad_bin_cutoff=1e-100,
                        has_hap = False,
):
    n_homo = P.alpha.shape[1]
    n_het = int(n_homo * (n_homo - 1) / 2)
    n_states = n_homo + n_het

    c = np.array([cont[r] for r in IX.OBS2RG])

    #read_emissions = np.ones((n_obs, 3))
    gt_emissions = np.ones((IX.n_snps, 3))
    _p_reads_given_gt(gt_emissions, P.O, P.N, P.P_cont, IX.n_obs, c, error,
                      snp_id=IX.OBS2SNP)


    GT = np.empty((IX.n_snps, n_states, 3))
    #P(GT | Z)
    for s in range(n_homo):
           _p_gt_homo(gt=GT[:, s, :], 
                      alpha=P.alpha[:, s] ,
                      beta=P.beta[:, s],
                      tau=tau[s],
                      n_snps = IX.n_snps)

    fa, fb = (P.alpha ).T, (P.beta ).T
    pp = fa / (fa + fb)
    qq = 1 - pp
    for s1 in range(n_homo):
        for s2 in range(s1+1, n_homo):
           s += 1
           _p_gt_het(GT[:, s, :], 
                     alpha1 = fa[s1],
                     beta1 = fb[s1],
                     alpha2 = fa[s2],
                     beta2 = fb[s2],
                     n_snps=IX.n_snps)

    snp_emissions = np.sum(GT * gt_emissions[:, np.newaxis, :],2)   

    bin_emissions = [np.ones((IX.bin_sizes[i], n_states)) for i in range(IX.n_chroms)]
    for (chrom, bin_), snp in zip(IX.SNP2CHROMBIN, snp_emissions):
        #print(chrom, bin_, snp)
        bin_emissions[chrom][bin_] *= snp

        if np.sum(bin_emissions[chrom][bin_]) < bad_bin_cutoff:
            bin_emissions[chrom][bin_] += bad_bin_cutoff / n_states

    return bin_emissions


#cdef void _p_read_given_snp(double [:,:] gt, double [:,:] read, long [:] snp_id,
#                            long n_snps):
#    cdef long i, g
#    for i in range(n_snps):
#        for g in range(3):
#            gt[ snp_id[i], g] *= read[i, g]
#    return

#def p_reads_given_gt_cy(P, n_obs, c, error, has_hap = False):
#    read_emissions = np.ones((n_obs, 5)) if has_hap else np.ones((n_obs, 3))
#    _p_reads_given_gt(read_emissions, P.O, P.N, P.P_cont, n_obs, c, error)
#    return read_emissions

    
#cdef void _p_read_given_snp(double [:,:] gt, double [:,:] read, long [:] snp_id,
#                            long n_snps):
#    cdef long i, g
#    for i in range(n_snps):
#        for g in range(3):
#            gt[ snp_id[i], g] *= read[i, g]
#    return

#def p_reads_given_gt_cy(P, n_obs, c, error, has_hap = False):
#    read_emissions = np.ones((n_obs, 5)) if has_hap else np.ones((n_obs, 3))
#    _p_reads_given_gt(read_emissions, P.O, P.N, P.P_cont, n_obs, c, error)
#    return read_emissions

    
