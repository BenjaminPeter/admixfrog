cdef double dbinom_single(int , int, double ) nogil
cdef double _dbeta_single(double p, double alpha, double beta) nogil
cdef double _dbeta_mu(double p, double mu, double M) nogil
cdef double _dbetabinom_single(int k, int N, double alpha, double beta) nogil
cdef double _dbetabinom_mu(int k, int N, double mu, double M) nogil
cdef double _logdbeta_single(double p, double alpha, double beta) nogil
cdef void _p_reads_given_gt(double[:,:] ll_geno, 
                            long[:] O, 
                            long[:] N, 
                            double[:] p_cont, 
                            int n_obs, 
                            double[:] c,
                            double error, 
                            long[:] snp_id) nogil
cdef void _p_gt_het(double[:,:] gt, double[:] alpha1, double[:] beta1, double[:] alpha2,
                    double[:] beta2, int n_snps) nogil

cdef void _p_gt_homo(double[:, :] gt, double[:] alpha, double[:] beta, double tau, int n_snps) nogil
cdef double _p_gt_homo_single(int g, double a, double b, double tau) nogil
cdef double _p_gt_het_single(int g, double a1, double b1, double a2, double b2) nogil
cdef double _log_p_gt_homo_single(int g, double a, double b, double tau) nogil
