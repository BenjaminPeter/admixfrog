#!python
#cython: language_level=3
#cython: infer_types=True

import pandas as pd
import numpy as np
cimport numpy as np
cimport scipy.special.cython_special as scs
cimport cython
from libc.math cimport pow, log, exp
from cython.parallel import prange
from libc.stdio cimport printf
from numpy.math cimport INFINITY
from scipy.optimize import minimize, minimize_scalar

ctypedef np.int_t INT_T
ctypedef np.float64_t DOUBLE_T

cdef double _pbinom_single(int k, int N, double p) nogil:
    return scs.binom(N, k) * pow(p,k) * pow(1.-p, N-k)


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef void _get_emissions(
    double[:, :] emissions,
    double [:] cont,
    long[:] O,
    long[:] N,
    double[:] p_cont,
    double[:, :] p_state,
    long[:] bin_id,
    int n_steps, int n_snps, int n_states,
    double error,
    ) nogil :
    """R * P(R = r | Z = k, o, c, n) = P(o | R, Z) P(R|Z c,n) / P(o | Z, c, n)

    return a vector of length n+1, giving the probabilty that r reads are contaminant
    """
    #cdef np.ndarray[DOUBLE_T, ndim=2] emissions = np.ones((n_steps, n_states))
    #cdef double[:, :] e = emissions
    cdef int s, i, row
    cdef double p
    #cdef np.ndarray[DOUBLE_T] p

    for s in range(n_states):
        #p = p_cont * cont + p_state[:, s] * (1.-cont)
        #p = p * (1. - error) + (1.-p) * error
        for i in range(n_snps):
            p = p_cont[i] * cont[i] + p_state[i, s] * (1.-cont[i])
            p = p * (1. - error) + (1.-p) * error
            row = bin_id[i]
            emissions[row, s] *= _pbinom_single(O[i], N[i], p)

    return 


def get_emissions_cy(cont, bins, bin_data, 
                     freqs, 
                     e=1e-2, 
                     bad_snp_cutoff=1e-10,
                     garbage_state=True):
    c = np.array([cont[l] for l in freqs.lib])
    bin_id = bin_data[:, 1]
    n_snps = len(freqs.O)
    n_steps = bins.shape[0]
    n_states = freqs.P.shape[1]

    emissions = np.ones((n_steps, n_states+garbage_state))
    
    _get_emissions(
                          emissions,
                            cont=c,
                          O=freqs.O,
                          N=freqs.N,
                          p_cont = freqs.P_cont,
                          p_state = freqs.P,
                          bin_id = bin_id,
                          n_snps = n_snps,
                          n_steps = n_steps,
                          n_states = n_states,
                          error = e)

    if garbage_state:
        emissions[:, n_states] = bad_snp_cutoff
    else:
        em_sum = np.sum(emissions, 1)
        bad_snps = em_sum < bad_snp_cutoff
        emissions[bad_snps, :] = bad_snp_cutoff / n_states

    chroms = np.unique(bins[:, 0])
    emissions = [emissions[bins[:,0] == chrom] for chrom in chroms]
    
    return emissions


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double get_po_given_c(
    double c, 
    double e,
    long [:] O,
    long [:] N,
    double [:] P_cont,
    double [:, :] P,
    double [:, :] G,
    long [:] ids,
    int n_states,
    ) nogil:
    """
    calculate Pr(O | Z, c, p), the probability of all observations given 
    the hidden states, allele frequencies in all potential donors and contamination rate

    - this is currenlty used to maximize the likelihood of c
    - the binomial coefficient is omitted
    - error is ignored

    cont : contamination estimate, by library
    freqs : allele frequencies, reads, by library/chromosome
    """
    cdef int i, s, n_snps, id_
    cdef double p, prob,  ll = 0.

    n_snps = len(ids) 
    for i in range(n_snps):
        prob = 0.
        id_ = ids[i]
        for s in range(n_states):
            p = c * P_cont[id_] + (1.-c) * P[id_, s]
            p = p * (1-e) + (1-p) * e
            prob += G[i,s] * pow(p, O[id_]) * pow(1-p, N[id_] - O[id_])

        ll += log(prob)
    return ll

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function

def split_freqs(libs, freqs):
    split_ids = dict()

    for lib in libs:
        split_ids[lib] = np.where(freqs.lib == lib)[0]

    return split_ids

def update_contamination_cy(cont, error, bin_data, 
                            freqs, gamma, 
                            #split_data,
                            split_ids,
                            libs,
                            garbage_state=True
                     ):
    """
    update emissions by maximizing contamination parameter

    cont: list of contamination rates (by library)
    gamma : Pr(Z | O, c_prev)
    snp_to_z : data structure connection snp to hidden state/ chrom coordinate
    data:

    we split up SNP 1. by library, 2. by sequence
    data is organized as data[lib_id][seq_id][state_id],
    i.e. the entry of data[lib_id][seq_id][state_id] is a Freq(O, N, p_cont, p_state)
        object


    """



    n_states = freqs.P.shape[1]

    n_libs = len(libs)
    for i in range(n_libs):
        lib = libs[i]
        f_ = split_ids[lib]
        G = np.array([gamma[i][j+1] for i, _, j in bin_data[f_]])
        assert all(lib == freqs.lib[f_])
        #print(lib, np.sum(G, 0), end = "\t")

        if garbage_state:
            G = G[:,:-1]

        def get_po_given_c_all(c):
            prob = get_po_given_c(c=c,
                                 e=error,
                                 O=freqs.O,
                                 N=freqs.N,
                                 P_cont=freqs.P_cont,
                                 P=freqs.P,
                                 G=G,
                                 ids = split_ids[lib],
                                 n_states = n_states)
            return -prob


        p0 = get_po_given_c_all(cont[lib])

        #OO =  minimize_scalar(get_po_given_c_all, bounds=(0., 1), method="Bounded")
        #print("[%s/%s]minimizing \tc: [%.4f->%.4f]:\t%.4f" % (lib, len(f_),
        #                                                           cont[lib], OO.x, p0-OO.fun))
        #cont[lib] = OO.x
        OO =  minimize(get_po_given_c_all, [cont[lib]], bounds=[(0., 1)])
        print("[%s/%s]minimizing \tc: [%.4f->%.4f]:\t%.4f" % (lib, len(f_),
                                                                   cont[lib], OO.x[0], p0-OO.fun))
        cont[lib] = OO.x[0]
    return dict(cont)

