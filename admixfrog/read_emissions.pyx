#!python
#cython: language_level=3
#cython: infer_types=True

import pandas as pd
import numpy as np
cimport scipy.special.cython_special as scs
cimport cython
from libc.math cimport pow, log, exp
from libc.stdio cimport printf
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln
from scipy.stats import binom
from .log import log_

q=0



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double get_po_given_c(
    double c, 
    double e,
    char [:] O,
    char [:] N,
    double [:] P_cont,
    double [:, :, :] PG,
    long [:] rg2obs,
    long [:] obs2snp
    ) :
    cdef int i, s, g,
    cdef long n_obs, n_states, 
    cdef int obs, snp
    cdef double p,  ll = 0.

    n_obs = len(rg2obs) 
    n_states = PG.shape[1]
    for i in range(n_obs):
        obs = rg2obs[i]
        snp = obs2snp[obs]
        for g in range(3):
            p = c * P_cont[obs] + (1.-c) * g / 2.
            p = p * (1-e) + (1-p) * e
            p = O[obs] * log(p) + (N[obs] - O[obs]) * log(1-p)
            for s in range(n_states):
                ll += PG[snp, s, g] * p
    return ll

def update_contamination(cont, error, P, PG, IX, libs):
    """
    update emissions by maximizing contamination parameter

    cont: dict of contamination rates (by library)
    PG: Pr(G, Z | O)



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
                                 PG=PG,
                                 rg2obs = IX.RG2OBS[lib],
                                 obs2snp = IX.OBS2SNP)
            return -prob


        p0 = get_po_given_c_all(cont[lib])

        #OO =  minimize_scalar(get_po_given_c_all, bounds=(0., 1), method="Bounded")
        #print("[%s/%s]minimizing \tc: [%.4f->%.4f]:\t%.4f" % (lib, len(f_),
        #                                                           cont[lib], OO.x, p0-OO.fun))
        #cont[lib] = OO.x
        OO =  minimize(get_po_given_c_all, [cont[lib]], bounds=[(0., 1-1e-10)])
        log_.info("[%s]\t%s\tc:[%.4f->%.4f]:\t%.4f" % (lib, len(f_),
                                                                   cont[lib], OO.x[0], p0-OO.fun))
        delta += abs(cont[lib] - OO.x[0])
        cont[lib] = OO.x[0]

    return delta
