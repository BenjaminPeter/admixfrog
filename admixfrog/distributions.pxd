cdef double dbinom_single(int , int, double ) nogil
cdef double _dbeta_single(double p, double alpha, double beta) nogil
cdef double _dbeta_mu(double p, double mu, double M) nogil
cdef double _dbetabinom_single(int k, int N, double alpha, double beta) nogil
cdef double _dbetabinom_mu(int k, int N, double mu, double M) nogil
cdef double _logdbeta_single(double p, double alpha, double beta) nogil

