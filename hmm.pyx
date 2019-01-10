#!python
#cython: language_level=3
# cython: infer_types=True
# distutils: language = c++

from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
cimport scipy.special.cython_special as scs
cimport cython
from libc.math cimport pow, log, exp
from libc.stdlib cimport malloc, free
from cython.parallel import prange
from libc.stdio cimport printf
from numpy.math cimport INFINITY

# Pr(O|R, Z, N)
#G_PO_RZ[ l, s, o, r] is Pr(O|R, Z) for observation o in library l, conditional on Z=s, R=r
cdef vector[vector[vector[vector[double]]]] G_PO_RZ

ctypedef np.int_t INT_T
ctypedef np.float64_t DOUBLE_T

cdef double _pbinom_single(int k, int N, double p) nogil:
    return scs.binom(N, k) * pow(p,k) * pow(1.-p, N-k)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double* _pbinom_all(
    int N, 
    double p) nogil:
    cdef double * ret_array = <double*>malloc( (N+1) * sizeof(double*))
    cdef double c
    cdef int k

    for k in range(N+1):
        ret_array[k] = _pbinom_single(k, N, p)
    return ret_array

# example R is 2, E is 5, O is 1 ==> K = 0,1    | O-K = 1,0
# example R is 2, E is 5, O is 4 ==> K = 0,1,2  | O-K = 4,3,2
# example R is 5, E is 2, O is 1 ==> K = 0,1    | O-K = 1,0
# example R is 5, E is 2, O is 4 ==> K = 2,3,4  | O-K = 2,1,0
# we know R out of N reads are contaminant, E are endogenous
# thus, if we observe O derived alleles we need to sum over all combination of
#endogenous / contaminant
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double _pbinom_convolution(int obs, int R, int E, double p_cont, double p_state) nogil:
    cdef double res = 0.
    cdef int k  # k is number of alt alleles from contaminants
    k_max = min(obs, R) 
    k_min = max(obs - E, 0) #q = obs - k is number of alt from endog, q_max is min(obs, E)
    for k in range(k_min, k_max+1):
        res += _pbinom_single(k, R, p_cont ) * _pbinom_single(obs-k, E, p_state)

    return res

cpdef double _pbinom_convolution2(int obs, int R, int E, double p_cont, double p_state) :
    cdef double res = 0.
    cdef int k  # k is number of alt alleles from contaminants
    k_max = min(obs, R) 
    k_min = max(obs - E, 0) #q = obs - k is number of alt from endog, q_max is min(obs, E)
    for k in range(k_min, k_max+1):
        res += _pbinom_single(k, R, p_cont ) * _pbinom_single(obs-k, E, p_state)

    return res

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(True)
cdef double _er_given_zoc_single_obs(
    double c, 
    int obs,
    int N,
    double p_cont,
    double p_state,
    vector[double] po_given_rz
    ) nogil:
    """sum R * P(R = r | Z = k, o, c, n) = P(o | R, Z) P(R|Z c,n) / P(o | Z, c, n)
    """
    cdef int r
    cdef double er_given_zoc = 0.
    cdef double psum = 0.
    cdef double po_given_rz_local = 0.
    # array of P(R=i | N, c)
    cdef double* pr_given_c = _pbinom_all(N, c)
    cdef double po_given_zc = _pbinom_single(obs, N, (1-c) * p_state + c * p_cont)
    for r in range(N+1): 
        po_given_rz_local = _pbinom_convolution(obs, r, N -r, p_cont, p_state)
        er_given_zoc += r * po_given_rz_local * pr_given_c[r]  / po_given_zc
        #er_given_zoc += r * po_given_rz[r] * pr_given_c[r]  / po_given_zc
        psum += po_given_rz[r] * pr_given_c[r]  / po_given_zc

    if er_given_zoc != er_given_zoc:
        er_given_zoc = 0#- INFINITY
        printf("\t%f\t%d\t%d\t%f\t%f\t%f\n", c, obs, N, p_cont, p_state, er_given_zoc)
    if psum - 1 > 1e-10:
        printf("ERROR:\t%f\n", psum)

    free(pr_given_c)
    return er_given_zoc


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef void init_global_po(obs, N, p_cont, p_state):
    cdef long n_states = p_state[0].shape[0]
    cdef long n_libs = len(obs)
    global G_PO_RZ
    cdef long l, s, i, r, n_obs, N_

    G_PO_RZ.resize(n_libs)
    for l in range(n_libs):
        G_PO_RZ[l].resize(n_states)
        for s in range(n_states):
            n_obs = obs[l].shape[0]
            G_PO_RZ[l][s].resize(n_obs)
            for i in range(n_obs):
                N_ = N[l][i]
                G_PO_RZ[l][s][i].resize(N_+1)
                for r in range(N_+1):
                    G_PO_RZ[l][s][i][r] = _pbinom_convolution(
                        obs[l][i], 
                        r,
                        N[l][i] - r,
                        p_cont[l][i],
                        p_state[l][s][i]
                    )
                #if i % 10000 == 0:
                #    print(l, s, i, G_PO_RZ[l][s][i])
    print("done init global")


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double* _er_given_zoc_all(
    double c,
    long[:] obs,
    long[:] N,
    double[:] p_cont,
    double[:] p_state,
    vector[vector[double]] po_given_rz
    ) :
    cdef int n_obs = obs.shape[0]
    cdef int i
    cdef double * er = <double*>malloc( n_obs * sizeof(double*))
    cdef double[:] x

    with nogil:
        for i in prange(n_obs, schedule="static", chunksize=10000):
            er[i] = _er_given_zoc_single_obs(c, 
                                             obs[i], 
                                             N[i], 
                                             p_cont[i],
                                             p_state[i], 
                                             po_given_rz[i])

    return er


cpdef np.ndarray[DOUBLE_T] er_given_zoc(
    double c,
    np.ndarray[INT_T] obs,
    np.ndarray[INT_T] N,
    np.ndarray[DOUBLE_T] p_cont,
    np.ndarray[DOUBLE_T] p_state,
    int lib, int state
    ) :
    """R * P(R = r | Z = k, o, c, n) = P(o | R, Z) P(R|Z c,n) / P(o | Z, c, n)

    return a vector of length n+1, giving the probabilty that r reads are contaminant
    """
    cdef double* er = _er_given_zoc_all(c, obs, N, p_cont, p_state, G_PO_RZ[lib][state])
    cdef long n_obs = obs.shape[0]
    cdef np.ndarray[DOUBLE_T] x = np.empty(n_obs, np.float64)
    for i in range(n_obs):
        x[i] = er[i]

    free(er)
    return x


