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
from scipy.special import betaln
from .distributions cimport *

ctypedef np.int_t INT_T
ctypedef np.float64_t DOUBLE_T

def posterior_geno_homo(state, P, snp_id, error, c, tau):
    n_snps, n_obs = P.alpha.shape[0], P.O.shape[0]
    ll_geno = np.ones((n_snps, 3))
    _p_reads_given_gt(ll_geno=ll_geno,
                      O=P.O, 
                      N=P.N,
                      p_cont=P.P_cont, 
                      n_obs = n_obs, 
                      c=c, 
                      error=error,
                      snp_id = snp_id)
    prior_geno = np.empty((n_snps, 3))
    _p_gt_homo(prior_geno,
               alpha=P.alpha[:, state] ,
               beta=P.beta[:, state],
               tau=tau,
               n_snps=n_snps)
    G =  ll_geno * prior_geno
    G = G / np.sum(G, 1)[:, np.newaxis]
    G[np.isnan(G)] = 1. / G.shape[1]
    return np.minimum(np.maximum(G,0), 1)#, ll_geno, prior_geno #rounding error
def posterior_geno_het(s1, s2, P, snp_id, error, c):
    n_snps, n_obs = P.alpha.shape[0], P.O.shape[0]
    ll_geno = np.ones((n_snps, 3))
    _p_reads_given_gt(ll_geno=ll_geno,
                      O=P.O, 
                      N=P.N,
                      p_cont=P.P_cont, 
                      n_obs = n_obs, 
                      c=c, 
                      error=error,
                      snp_id = snp_id)
    prior_geno = np.ones((n_snps, 3))
    _p_gt_het(prior_geno,
               alpha1=P.alpha[:, s1],
               beta1=P.beta[:, s1],
               alpha2=P.alpha[:, s2],
               beta2=P.beta[:, s2],
               n_snps=n_snps)
    G =  ll_geno * prior_geno
    G = G / np.sum(G, 1)[:, np.newaxis]
    G[np.isnan(G)] = 1. / G.shape[1]
    return np.minimum(np.maximum(G,0), 1) #rounding error


def pll_F(tau,  gamma, postg, a, b, n_snps):
    cdef double LL = 0.
    cdef int i, j, g
    for i in range(n_snps):
        for g in range(3):
            a_, b_ = a[i] * tau, b[i] * tau
            LL += betaln(a_ + g, b_ + 2 - g) - betaln(a_, b_) + log(2) * (g == 1) #binom(2,g) if g==1
    return LL
#@cython.boundscheck(False) # turn off bounds-checking for entire function
#@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double ll_F(double tau, 
                 double[:] gamma, 
                 double[:, :] postg, 
                 double[:] a, 
                 double[:] b,
                 int n_snps) nogil:
    cdef double LL = 0.
    cdef int i, g
    for i in range(n_snps):
        for g in range(3):
            LL += _log_p_gt_homo_single(g, a[i], b[i], tau) * gamma[i] * postg[i, g]
    return LL


def update_F(tau, gamma, pg, P, IX, error, cont):

    c = np.array([cont[l] for l in IX.OBS2RG])
    n_states = len(tau)
    n_snps = P.alpha.shape[0]
    G = np.array([gamma[c][b+1] for c, b in IX.SNP2BIN])


    for s in range(n_states):
        def f(tau):
            x= ll_F(tau,
                    gamma=G[:, s],
                    postg=pg[:,s, :],
                  a=P.alpha[:, s],
                  b=P.beta[:, s],
                  n_snps=n_snps)
            return -x
        OO =  minimize(f, [tau[s]], bounds=[(1e-4, 1.-1e-4)])
        print("[%s]minimizing \tF: [%.4f->%.4f]:\t%.4f" % (s, tau[s], OO.x[0], OO.fun))
        tau[s] = OO.x[0]

def post_geno(P, cont, tau, IX, error):
    n_homo = len(tau)
    n_het = int(n_homo * (n_homo - 1) / 2)
    n_snps = P.alpha.shape[0]
    n_states = n_homo + n_het
    c = np.array([cont[l] for l in P.lib])

    pg = np.empty((n_snps, n_states, 3))
    for s in range(n_homo):
        pg[:, s] = posterior_geno_homo(s, P, IX.OBS2SNP,
                        error=error, c=c, tau=tau[s])
    s = n_homo
    for s1 in range(n_homo):
        for s2 in range(s1+1,n_homo):
            pg[:, s] = posterior_geno_het(s1, s2, P, IX.OBS2SNP,
                        error=error, c=c)
            s+=1
    return pg

def update_pars(tau, gamma, P, IX,
                error, cont, libs, split_ids):

    n_homo = len(tau)
    n_het = int(n_homo * (n_homo - 1) / 2)
    
    pg = post_geno(P=P, 
                   cont=cont,
                   tau=tau,
                   IX = IX,
                   error=error)
    assert np.all(pg >=0)
    assert np.all(pg <=1)
    assert np.allclose(np.sum(pg, 2), 1)

    #G = np.array([gamma[c][b+1] for c, b in bin_data[['chrom_id', 'loc_id']]])
    update_contamination(cont, error, P=P,
                         gamma=gamma, pg=pg,
                         IX=IX,
                         libs=libs)
    #update_F(tau, gamma, pg[:, :n_homo], P, bin_data, error, cont)



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double get_po_given_c(
    double c, 
    double e,
    long [:] O,
    long [:] N,
    double [:] P_cont,
    double [:, :] Z,
    double [:, :, :] pg,
    long [:] rg2obs,
    long [:] obs2bin,
    long [:] obs2snp,
    ) :
    cdef int i, s, g,
    cdef long n_obs, n_states, 
    cdef int obs, bin_, snp
    cdef double p,  ll = 0.

    n_obs = len(rg2obs) 
    n_states = pg.shape[1]
    for i in range(n_obs):
        obs = rg2obs[i]
        bin_ = obs2bin[obs]
        snp = obs2snp[obs]
        for g in range(3):
            p = c * P_cont[obs] + (1.-c) * g / 2.
            p = p * (1-e) + (1-p) * e
            p = O[obs] * log(p) + (N[obs] - O[obs]) * log(1-p)
            for s in range(n_states):
                ll += Z[bin_, s] * pg[snp, s, g] * p
    return ll

def update_contamination(cont, error, P, Z, pg, IX, libs):
    """
    update emissions by maximizing contamination parameter

    cont: dict of contamination rates (by library)
    gamma : Pr(Z | O)
    postg: Pr(G | Z, O)



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

