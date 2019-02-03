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

cdef double dbinom_single(int k, int N, double p) nogil:
    return scs.binom(N, k) * pow(p,k) * pow(1.-p, N-k)

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


def dbetabinom(long[:] k, long[:] N, double[:] a, double[:] b):
    res = np.empty(k.shape[0])
    for i in range(k.shape[0]):
        res[i] = _dbetabinom_single(k[i], N[i], a[i], b[i])
    return res

