#!python
# cython: language_level=3
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import pandas as pd
import numpy as np
cimport numpy as np
cimport scipy.special.cython_special as scs
cimport cython
from libc.math cimport pow, log, exp, lgamma
from libc.stdio cimport printf
from scipy.optimize import minimize, minimize_scalar

cdef double dbinom_single(int k, int N, double p) nogil:
    return scs.binom(N, k) * pow(p,k) * pow(1.-p, N-k)

cdef double lbeta(double a, double b) nogil:
    return lgamma(a) + lgamma(b) - lgamma(a + b);

cdef double _dbeta_single(double p, double alpha, double beta) nogil:
    return pow(p, alpha-1) * pow(1-p, beta - 1) / scs.beta(alpha, beta) 

cdef double _dbeta_mu(double p, double mu, double M) nogil:
    return _dbeta_single(p, mu * M, (1 - mu) * M)

cdef double _dbetabinom_single(int k, int N, double alpha, double beta) nogil:
    cdef double x
    x = lbeta(alpha + k, beta + N - k) - lbeta(alpha, beta)
    #x = scs.betaln(alpha + k, beta + N - k) - scs.betaln(alpha, beta)
    return scs.binom(N, k) * exp(x)

cdef double _dbetabinom_mu(int k, int N, double mu, double M) nogil:
    return _dbetabinom_single(k, N, mu *M, (1-mu) * M)

cdef double _logdbeta_single(double p, double alpha, double beta) nogil:
    return log(p) * (alpha-1) + log(1-p) * (beta - 1) - scs.betaln(alpha, beta) 

cdef double _p_gt_homo_single(int g, double a, double b, double tau) nogil:
    return _dbetabinom_single(g, 2, a * tau, b * tau)

cdef double _log_p_gt_homo_single(int g, double a, double b, double tau) nogil:
    cdef double a_, b_
    a_, b_ = a * tau, b * tau
    return lbeta(a_ + g, b_ + 2 - g) - lbeta(a_, b_) + log(2) * (g == 1) #binom(2,g) if g==1

cdef double _p_gt_het_single(int g, double a1, double b1, double a2, double b2) nogil:
    cdef double D = (a1 + b1) * (a2 + b2)
    if g == 0:
        return b1 * b2 / D
    if g == 2:
        return a1 * a2 / D
    return (a1 * b2 + a2 * b1) / D

cdef void _p_gt_homo(
    double[:, :] gt,
    double[:] alpha,
    double[:] beta,
    double tau,
    int n_snps
    ) nogil :
    cdef int i, g

    for i in range(n_snps):
        for g in range(3):
            gt[i, g] = _p_gt_homo_single(g, alpha[i], beta[i], tau)
    return 

cdef void _p_gt_het(
    double[:, :] gt,
    double[:] alpha1,
    double[:] beta1,
    double[:] alpha2,
    double[:] beta2,
    int n_snps
    ) nogil :

    cdef int i
    cdef double D

    for i in range(n_snps):
        D = (alpha1[i] + beta1[i]) * (alpha2[i] + beta2[i])
        gt[i, 0] = beta1[i] * beta2[i] / D
        gt[i, 2] = alpha1[i] * alpha2[i] / D
        gt[i, 1] = 1- gt[i,0] - gt[i, 2]
    return 

cdef void _p_reads_given_gt(double[:,:] ll_geno,
                      long [:] O,
                      long [:] N,
                      double[:] p_cont,
                      int n_obs,
                      double[:] c,
                      double error,
                      long[:] snp_id
                      ) nogil:
        
    cdef int g, s
    cdef double p

    for i in range(n_obs):
        for g in range(3):
            ll_geno[snp_id[i] ,g] *= _p_reads_given_gt_single(O[i],
                                                              N[i],
                                                              g,
                                                              p_cont[i],
                                                              c[i],
                                                              error)
    return

cdef double _p_reads_given_gt_single(long O, long N, int g, double p_cont, double c,
                                     double error) nogil:
    cdef double p
    p = c * p_cont + (1.-c) * g / 2.
    p = p * (1.-error) + (1.-p) * error
    return dbinom_single(O, N, p)

def P_reads_given_gt_single(O, N, g, p_cont, c, error):
    _p_reads_given_gt_single(O, N, g, p_cont, c, error)

def dbetabinom(long[:] k, long[:] N, double[:] a, double[:] b):
    res = np.empty(k.shape[0])
    for i in range(k.shape[0]):
        res[i] = _dbetabinom_single(k[i], N[i], a[i], b[i])
    return res

cdef double _log_p_gt_homo_single2(int g, double a, double b) nogil:
    return lbeta(a + g, b + 2 - g) - lbeta(a, b) + log(2) * (g == 1) #binom(2,g) if g==1
cdef double _p_gt_homo_single2(int g, double a, double b) nogil:
    return exp(_log_p_gt_homo_single2(g, a, b))
def dbetabinom2(double[:] a, double[:] b, long n_snps, double[:, :] res):
    cdef long i, g
    for i in range(n_snps):
        for g in range(3):
            res[i, g] = _p_gt_homo_single2(g, a[i], b[i])
