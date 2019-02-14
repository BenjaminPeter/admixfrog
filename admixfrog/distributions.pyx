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

def gt_homo_dist(double[:] a, double[:] b, double F, long n_snps, double[:, :] res):
    _gt_homo_dist(a, b, F, n_snps, res)
cdef void _gt_homo_dist(double[:] a, double[:] b, double F, long n_snps, double[:, :] res):
    cdef long i, g
    for i in range(n_snps):
        res[i, 0] =  (b[i]*b[i] + b[i] + a[i]*b[i]*F) / (a[i] + b[i]) / (a[i] + b[i] + 1)
        res[i, 2] =  (a[i]*a[i] + a[i] + a[i]*b[i]*F) / (a[i] + b[i]) / (a[i] + b[i] + 1)
        res[i, 1] =  (2.*a[i]*b[i]*(1-F)) / (a[i] + b[i]) / (a[i] + b[i] + 1)

