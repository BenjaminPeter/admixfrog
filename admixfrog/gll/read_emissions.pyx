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
from ..utils.log import log_



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

def update_contamination(cont, error, P, PG, IX,
                         est_options):
    """
    update emissions by maximizing contamination parameter

    cont: dict of contamination rates (by read group)
    error: dict of contamination rates (by read group)
    P: data structure with reference data
    PG: Pr(G, Z | O)

    """
    delta = 0.

    for i in range(len(IX.rgs)):
        rg = IX.rgs[i]
        f_ = IX.RG2OBS[rg]
        assert all(rg == P.rg[f_])

        def get_po_given_c_all(args):
            args = list(args)
            C = args.pop(0) if est_options['est_contamination'] else cont[rg]
            E = args.pop(0) if est_options['est_error'] else error[rg]

            prob = get_po_given_c(c=C,
                                 e=E,
                                 O=P.O,
                                 N=P.N,
                                 P_cont=P.P_cont,
                                 PG=PG,
                                 rg2obs = IX.RG2OBS[rg],
                                 obs2snp = IX.OBS2SNP)
            return -prob

        init, bounds = [], []
        if est_options['est_contamination']:
            init.append(cont[rg])
            bounds.append((0, 1-1e-10))
        if est_options['est_error']:
            init.append(error[rg])
            bounds.append((0.00001, .1))

        prev = get_po_given_c_all(init)
        OO =  minimize(get_po_given_c_all, init, bounds=bounds, method="L-BFGS-B")
        opt = OO.x.tolist()

        old_c, old_e = cont[rg], error[rg]
        if est_options['est_contamination']:
            cont[rg] = opt.pop(0)

        if est_options['est_error']:
            error[rg] = opt.pop(0)

        log__ = "[%s|%s] \tc: [%.4f->%.4f]\t:" % (rg, len(f_), old_c, cont[rg])
        log__ += "e: [%.4f->%.4f]:\t%.4f" % (old_e, error[rg], prev - OO.fun)
        log_.info(log__)
        delta += abs(cont[rg] - old_c) + abs(error[rg] - old_e)

    return delta
