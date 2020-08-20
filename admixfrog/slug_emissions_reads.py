import numpy as np
from numba import njit
from math import lgamma, exp, pow
from copy import deepcopy


def slug_p_gt_diploid(tau0, F0, SNP2SFS, FLIPPED=None, res=None):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l})

    Parameters
    ---------
    F0: np.ndarray[S x 1] float
        F-parameter measuring inbreeding, of length S
    tau0: np.ndarray[S x 1] float
        expected derived allele frequency, of length S
    SNP2SFS : np.ndarray[L x 1] int
        index array assigning each SNP to it's SFS entry
    FLIPPED : np.ndarray[L x 1] bool, optional
        index array asserting whether a SNP is flipped with respect to reference allele

    Returns
    -------
    res :  np.ndarray[L x 3]
        Pr(G_i = j  | F, tau)
    """

    if res is None:
        res = np.empty((FLIPPED.shape[0], 3))

    F = np.array(F0) if len(F0) == 1 else F0[SNP2SFS]
    tau = np.array(tau0) if len(tau0) == 1 else tau0[SNP2SFS]
    if FLIPPED is not None:
        tau[FLIPPED] = 1 - tau[FLIPPED]

    res[:, 0] = F * (1 - tau) + (1 - F) * (1 - tau) ** 2  # HOMO REF
    res[:, 2] = F * tau       + (1 - F) * tau * tau  # HOMO ALT
    res[:, 1] = (1 - F) * 2 * tau * (1 - tau)

    assert np.allclose(np.sum(res, 1), 1)
    np.clip(res, 0, 1, out=res)  # rounding

    return res


def slug_p_gt_haploid(*args, **kwargs):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l}) for haploid genotypes
    """
    res = slug_p_gt_diploid(F0=[1.], *args, **kwargs)
    return res 


def slug_fwd_p_g(data, pars):
    """calculate forward probabilities of genotypes"""
    fwd_g = slug_p_gt_diploid(pars.tau, pars.F, 
                              FLIPPED=data.FLIPPED, SNP2SFS=data.SNP2SFS)  # size [L x 3]
    if data.haploid_snps is not None:
        fwd_g[data.haploid_snps] = slug_p_gt_haploid(pars.tau, 
                                                     FLIPPED=data.FLIPPED[data.haploid_snps], 
                                                     SNP2SFS=data.SNP2SFS[data.haploid_snps])  
    return fwd_g


def slug_bwd_p_o_given_x(READS, e, b):
    """Caluclate Pr(O | X, e, b) 

    Parameters
    ----------
    READS : np.ndarray[R x 1] int 
    e : float
        error rate from reference to alt allele
    b : float
        error rate from alt to reference allele

    Returns
    -------
    res : np.ndarray [R x 1]
        Pr(O_i | X_i = j)
    """
    res = np.empty((READS.shape[0], 2))
    res[READS==1, 0] = e
    res[READS==1, 1] = 1 - b
    res[READS==0, 0] = 1 - e
    res[READS==0, 1] = b
    return res


@njit
def slug_bwd_p_one_o_given_g(bwd_x, fwd_a, fwd_c,  READ2SNP, READ2RG, n_reads):
    """calculate Pr(O_i | G ) 
    
    Parameters
    ----------
    bwd_x : np.ndarray[R x 2]
        Pr(O | X_i = j)
    fwd_a : np.ndarray[L x 1]
        Pr(X_i = 1 | C_r =1)
    fwd_c : np.ndarray[R x 1] 
        Pr(C_i)
    READ2SNP : np.ndarray[R x 1]
        index mapping each read to a SNP
    READ2RG : np.ndarray[R x 1] int
        index array mapping reads to corresponding read group
    n_reads : int
        number of reads R
    
    Returns
    -------
    Pr(O_i | G) : np.ndarray[R x 1]
        these are the naive probabilities of observations given genotypes of the 
        underlying SNP, not taking info from other Reads into account
    """

    fwd_a_mat = np.vstack((1-fwd_a, fwd_a)).T[READ2SNP]
    p_o_cont = np.sum(fwd_a_mat * bwd_x, 1) * fwd_c[READ2RG]

    # [i, j, k] Pr(O_i= k | G_i = j )
    bwd_g = np.zeros((n_reads, 3))
    bwd_g[:] = np.expand_dims(p_o_cont, 1) #cont is independent of genotype


    z = np.arange(3)/2
    M = np.vstack((1-z, z))
    bwd_g += (bwd_x @ M) * np.expand_dims((1-fwd_c[READ2RG]) , 1)

    return bwd_g

@njit
def slug_bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps):
    """
    calculate  ℙ(O_l | G_l) = ∏_r ∏_j ℙ(O_lrj | G_l)

    where ℙ(O_lrj = 1 | G_l) = (c_r)psi_l * (1-c_r) G_l / 2. 
    and   ℙ(O_lrj = 0 | G_l) = (c_r) (1 - psi_l) * (1-c_r) (1 - G_l / 2). 
    """
    bwd_g = np.ones((n_snps, 3))

    for read, snp in enumerate(READ2SNP):
        bwd_g[snp] *= bwd_g1[read]

    return bwd_g

def slug_fwd_p_x_cont(psi, READ2SNP):
    """slug_fwd_p_x_cont calculate Pr(X | C=1, psi)
    
    Parameters
    ----------
    psi : np.ndarray[L x 1]
        alt allele rate of contaminant
    READ2SNP : np.ndarray[R x 1]
        index mapping each read to a SNP
    
    Returns
    -------
    Pr(X_i = j | C=1, psi_l) : np.ndarray[R x 1]
    """
    res = np.empty((READ2SNP.shape[0], 2))

    res[:, 0] = 1 - psi[READ2SNP]
    res[:, 1] = psi[READ2SNP]
    return res


def message_fwd_p_x_nocont(fwd_g, bwd_g, bwd_g1, READ2SNP):
    """message passing version of Pr(X_i | C=0, G_l, O)

    calculate Pr(X_i | C=0, G_l, O) needed for updates;
    taking into account that other reads will change our beliefe about a particular read
    
    Parameters
    ----------
    fwd_g : np.ndarray [L x 1]
        Pr(G_l | Z)
    bwd_g : np.ndarray [L x 1]
        Pr(O_i1, O_i2, O_i3 | G_i) for all observations from this SNP
    bwd_g1 : np.ndarray [R x 1]
        Pr(O_i | G) for a single SNP
    READ2SNP : np.ndarray[R x 1]
        index mapping each read to a SNP
    
    Returns
    -------
    Pr(X_li | C=0, G_l, O_l1, O_l2, O_l...) : np.ndarray[R x 1]
    """
    v = np.arange(3)/2.; 
    M = np.vstack((1-v, v)).T

    res = ((fwd_g * bwd_g)[READ2SNP]  / (bwd_g1 + 1e-300) @ M)
    res = res / np.expand_dims(np.sum(res, 1), 1)

    return res

def slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, READ2RG):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    PC : array[n_obs]  of c_lr, the probability that the r'th read is
    contaminant
    PA : array[n_obs] the alt allele freq for observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = Σ_G Pr(X | G, C=0) Pr(G, C=0) + Σ_A Pr(X |A, C=1) Pr(A, C=1)
    """

    res = (1 - fwd_c[READ2RG])[:, np.newaxis] * fwd_x_nocont
    res += fwd_c[READ2RG][:, np.newaxis] * fwd_x_cont

    return res


def slug_post_g(bwd_g, fwd_g):
    """calculates marginalized posteriors for g, a and c

        post_g :   Pr(G_l | Z_l) x  (𝚷_rj Pr(O_lrj | G_l))    [L x 3]
    """

    post_g = bwd_g * fwd_g
    post_g /= np.expand_dims(np.sum(post_g, 1), 1)
    assert np.allclose(np.sum(post_g, 1), 1)

    return post_g

def slug_post_x(bwd_x, *args, **kwargs):
    fwd_x = slug_fwd_p_x(*args, **kwargs)
    post_x = fwd_x * bwd_x
    return post_x / np.sum(post_x, 1)[:, np.newaxis]


def full_posterior_genotypes(data, pars):
    fwd_g = slug_fwd_p_g(data, pars)
    bwd_x = slug_bwd_p_o_given_x(data.READS, e=pars.e, b = pars.b)
    bwd_g1 = slug_bwd_p_one_o_given_g(bwd_x, data.psi, pars.cont, 
                                       data.READ2SNP, data.READ2RG, data.n_reads)
    bwd_g = slug_bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    return bwd_g, slug_post_g(bwd_g, fwd_g)



#@njit
def slug_post_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, READ2RG):
    """ calculate Pr(C | O) = Pr(C, O) / Pr (O)
        with Pr(C, O)  = Σ_X  Pr(O | X) Pr(X | C) Pr(C)

        return array [i, j k] with Pr(C_i = k | O_i = j)
    """
    x_cont = np.sum(fwd_x_cont * bwd_x, 1) * fwd_c[READ2RG]
    x_nocont = np.sum(fwd_x_nocont * bwd_x,1) * (1 - fwd_c[READ2RG])
    
    post_c = x_cont / (x_cont +x_nocont)
    if np.any(np.isnan(post_c)):
        raise ValueError()

    return post_c


def calc_full_ll_reads(data, pars):
    fwd_g = slug_fwd_p_g(data, pars)
    bwd_x = slug_bwd_p_o_given_x(data.READS, e=pars.e, b=pars.b)  # Pr(O_lrj | G_lrj)
    bwd_g1 = slug_bwd_p_one_o_given_g(bwd_x, data.psi, pars.cont, data.READ2SNP, data.READ2RG, data.n_reads)
    bwd_g = slug_bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    ll2 = calc_ll(fwd_g, bwd_g)
    return ll2


def calc_ll(fwd_g, bwd_g):
    ll = np.sum(np.log(np.sum(fwd_g * bwd_g, 1)))
    return ll



