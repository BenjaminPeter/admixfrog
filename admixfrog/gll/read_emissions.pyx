#!python
#cython: language_level=3
#cython: infer_types=True

import numpy as np
cimport cython
from libc.math cimport pow, log, exp
from scipy.optimize import minimize
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
    ) :
    cdef int i, s, g,
    cdef long n_obs, n_states, 
    cdef double p,  ll = 0.

    n_obs = len(O) 
    n_states = PG.shape[1]
    for i in range(n_obs):
        for g in range(3):
            p = c * P_cont[i] + (1.-c) * g / 2.
            p = p * (1-e) + (1-p) * e
            p = O[i] * log(p) + (N[i] - O[i]) * log(1-p)
            for s in range(n_states):
                ll += PG[i, s, g] * p
    return ll

def update_contamination(cont, error, P, PG,
                         est_options):
    """
    update emissions by maximizing contamination parameter

    cont: array with contamination rates (for each read group)
    error: array with contamination rates (for each read group)
    P: data structure with reference data
    PG: Pr(G, Z | O)

    """
    delta = 0.

    for rg, i in P.rgs.items():
        f_ = P.OBS2RG == i
        log_.debug(np.sum(f_))

        def get_po_given_c_all(args):
            args = list(args)
            C = args.pop(0) if est_options['est_contamination'] else cont[i]
            E = args.pop(0) if est_options['est_error'] else error[i]

            prob = get_po_given_c(c=C,
                                 e=E,
                                 O=P.O[f_],
                                 N=P.N[f_],
                                 P_cont=P.psi[f_],
                                 PG=PG[P.OBS2SNP[f_]])
            return -prob

        init, bounds = [], []
        if est_options['est_contamination']:
            init.append(cont[i])
            bounds.append((0, 1-1e-10))
        if est_options['est_error']:
            init.append(error[i])
            bounds.append((0.00001, .1))

        prev = get_po_given_c_all(init)
        OO =  minimize(get_po_given_c_all, init, bounds=bounds, method="L-BFGS-B")
        opt = OO.x.tolist()

        old_c, old_e = cont[i], error[i]
        if est_options['est_contamination']:
            cont[i] = opt.pop(0)

        if est_options['est_error']:
            error[i] = opt.pop(0)

        log__ = "[%s|%s] \tc: [%.4f->%.4f]\t:" % (rg, np.sum(f_), old_c, cont[i])
        log__ += "e: [%.4f->%.4f]:\t%.4f" % (old_e, error[i], prev - OO.fun)
        log_.info(log__)
        delta += abs(cont[i] - old_c) + abs(error[i] - old_e)

    return delta
