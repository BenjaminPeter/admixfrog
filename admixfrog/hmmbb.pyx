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

ctypedef np.int_t INT_T
ctypedef np.float64_t DOUBLE_T

cdef double _dbeta_single(double p, double alpha, double beta) nogil:
    return pow(p, alpha-1) * pow(1-p, beta - 1) / scs.beta(alpha, beta) 

cdef double _dbeta_mu(double p, double mu, double M) nogil:
    return _dbeta_single(p, mu * M, (1 - mu) * M)

cdef double _dbetabinom_single(int k, int N, double alpha, double beta) nogil:
    cdef double x
    x = scs.betaln(alpha + k, beta + N - k) - scs.betaln(alpha, beta)
    return scs.binom(N, k) * exp(x)

cdef double _dbetabinom_mu(int k, int N, double mu, double M) nogil:
    return _dbetabinom_single(k, N, mu *M, (1-mu) * M)

cdef double _logdbeta_single(double p, double alpha, double beta) nogil:
    return log(p) * (alpha-1) + log(1-p) * (beta - 1) - scs.betaln(alpha, beta) 

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef void _get_emissions_beta_homo(
    double[:, :] emissions,
    double [:] cont,
    long[:] O,
    long[:] N,
    double[:] p_cont,
    double[:, :] alpha,
    double[:, :] beta,
    double[:] fst,
    long[:] snp_id,
    int n_pops, int n_snps, 
    double error,
    #double tau_state1, double tau_state2, #drift to reference
    #double tau_cont #drift to reference
    ) nogil :
    """model is now that each read is - 
        - contaminant w. p. c
        - S1 w.p. (1-c) /2
        - S2 w.p. (1-c) / 2
        
        P(O=1) = c p_c + (1-c)/2 Beta(a_1, M_1) + (1-c)/2 Beta(f_2, M_2)

        - c is shared between libraries
        - M is shared between refernce panels

    """
    cdef int s, i, g
    cdef int row = -1
    cdef double a, b, p
    cdef double G, D = 0.


    for s in range(n_pops):
        for i in range(n_snps):
            a, b = alpha[s, i] * fst[s], beta[s, i] * fst[s]
            for g in range(3):
                G = _dbetabinom_single(g, 2, a, b)
                p = cont[i] * p_cont[i] + (1.-cont[i]) / 2. * g
                p = p * (1. - error) + (1 - p) * error
                row = snp_id[i]
                emissions[row, s] += G  * _dbinom_single(O[i], N[i], p)

    return 

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef void _get_emissions_beta_het(
    double[:, :] emissions,
    double [:] cont,
    long[:] O,
    long[:] N,
    double[:] p_cont,
    double[:, :] alpha,
    double[:, :] beta,
    long[:] snp_id,
    int n_pops, int n_snps, 
    double error,
    #double tau_state1, double tau_state2, #drift to reference
    #double tau_cont #drift to reference
    ) nogil :
    """model is now that each read is - 
        - contaminant w. p. c
        - S1 w.p. (1-c) /2
        - S2 w.p. (1-c) / 2
        
        P(O=1) = c p_c + (1-c)/2 Beta(a_1, M_1) + (1-c)/2 Beta(f_2, M_2)

        - c is shared between libraries
        - M is shared between refernce panels

    """
    cdef int s, s1, s2, i, row
    cdef double p0, p1, p2 # probability of derivd allele
    cdef double G0, G1, G2 #Prob( GT)
    cdef double D

    s = n_pops
    for s1 in range(n_pops):
        for s2 in range(s1+1, n_pops):
            s += 1 # next column
            for i in range(n_snps):
                D = (alpha[s1, i] + beta[s1, i]) * (alpha[s2, i] + beta[s2, i])
                G0 = alpha[s1, i] * alpha[s2, i] / D
                G1 = (alpha[s1, i] * beta[s2, i] + alpha[s2, i] * beta[s1, i]) / D
                G2 = beta[s1, i] * beta[s2, i] / D

                p0 = cont[i] * p_cont[i]
                p1 = cont[i] * p_cont[i] + (1.-cont[i]) / 2.
                p2 = cont[i] * p_cont[i] + (1.-cont[i])

                p0 = p0 * (1. - error) + (1 - p0) * error
                p1 = p1 * (1. - error) + (1 - p1) * error
                p2 = p2 * (1. - error) + (1 - p2) * error

                D = G0  * _dbinom_single(O[i], N[i], p0)
                D += G1  * _dbinom_single(O[i], N[i], p1)
                D += G2  * _dbinom_single(O[i], N[i], p2)

                row = snp_id[i]
                emissions[row, s] *= D
    return 

def get_emissions_beta_cy(cont, bins, bin_data, 
                     freqs, #data object
                     fst,
                     snp_id,
                     e=1e-2, 
                     bad_snp_cutoff=1e-250,
                     garbage_state=False):
    """
    freqs needs:
    let n = n_snps, p = n_pops, l = n_libraries
        a, b for each pop [n x p]
        O, N for each snp [n x 1]
        cont for each library [l x 1]
        p_cont [n x 1]
    """
    c = np.array([cont[l] for l in freqs.lib])
    bin_id = bin_data[:, 1]

    n_snps = len(freqs.O)
    n_steps = bins.shape[0]
    n_pops = freqs.A.shape[1]
    n_states = int(n_pops + n_pops * (n_pops - 1) / 2)

    emissions = np.ones((n_steps, n_states+garbage_state, 3))
    
    _get_emissions_beta_homo(
                          emissions,
                            cont=c,
                          O=freqs.O,
                          N=freqs.N,
                          p_cont = freqs.P_cont,
                          alpha = freqs.A,
                          beta = freqs.B,
                          snp_id = snp_id,
                          n_snps = n_snps,
                          n_pops = n_pops,
                          fst = fst,
                          error = e)

    _get_emissions_beta_het(
        emissions=emissions,
        cont=c, #for each SNP/lib
        O=freqs.O,
        N=freqs.N,
        p_cont = freqs.P_cont,
        alpha = freqs.A,
        beta = freqs.B,
        snp_id = snp_id,
        n_snps = n_snps,
        n_pops = n_pops,
        error = e)

    if garbage_state:
        emissions[:, n_states, :] = bad_snp_cutoff
    else:
        em_sum = np.sum(emissions, 1)
        bad_snps = em_sum < bad_snp_cutoff
        emissions[bad_snps, :] = bad_snp_cutoff / n_states

    chroms = np.unique(bins[:, 0])
    emissions = [emissions[bins[:,0] == chrom] for chrom in chroms]
    
    return emissions

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double get_po_given_c_beta(
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
            prob += G[i,s] * (  O[id_] * log(p) + (N[id_] - O[id_]) * log(1-p))

        ll += prob
    return ll

def update_contamination_beta_cy(cont, error, bin_data, 
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
