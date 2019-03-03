#cython: language_level=3
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

cimport scipy.special.cython_special as scs
cimport cython
from libc.math cimport pow, log, exp, lgamma

cdef double dbinom_single(int k, int N, double p) nogil:
    return scs.binom(N, k) * pow(p,k) * pow(1.-p, N-k)

def gt_homo_dist(double[:] a, double[:] b, double F, double tau, long n_snps, double[:, :] res):
    _gt_homo_dist_tau(a, b, F, tau, n_snps, res)
cdef void _gt_homo_dist(double[:] a, double[:] b, double F, long n_snps, double[:, :] res):
    cdef long i, g
    for i in range(n_snps):
        res[i, 0] =  (b[i]*b[i] + b[i] + a[i]*b[i]*F) / (a[i] + b[i]) / (a[i] + b[i] + 1)
        res[i, 2] =  (a[i]*a[i] + a[i] + a[i]*b[i]*F) / (a[i] + b[i]) / (a[i] + b[i] + 1)
        res[i, 1] =  (2.*a[i]*b[i]*(1-F)) / (a[i] + b[i]) / (a[i] + b[i] + 1)

cdef void _gt_homo_dist_tau(double[:] a, double[:] b, double F, double tau, long n_snps, double[:, :] res):
    cdef long i, g
    for i in range(n_snps):
        res[i, 0] =  (b[i]*b[i]*tau + b[i] + a[i]*b[i]*F * tau) / (a[i] + b[i]) / (tau *a[i] + tau *b[i] + 1)
        res[i, 2] =  (a[i]*a[i]*tau + a[i] + a[i]*b[i]*F * tau) / (a[i] + b[i]) / (tau *a[i] + tau *b[i] + 1)
        res[i, 1] =  1 - res[i, 0] - res[i, 2] #(2.*a[i]*b[i]*(1-F) * (1-tau)) / (a[i] + b[i]) / (a[i] + b[i] + 1)

