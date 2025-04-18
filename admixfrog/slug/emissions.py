import numpy as np
from numba import njit
from math import lgamma, exp, pow
from copy import deepcopy


@njit(cache=True)
def p_gt_diploid(tau0, F0, SNP2SFS, res=None):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l})

    Parameters
    ---------
    F0: np.ndarray[S x 1] float
        F-parameter measuring inbreeding, of length S
    tau0: np.ndarray[S x 1] float
        expected derived allele frequency, of length S
    SNP2SFS : np.ndarray[L x 1] int
        index array assigning each SNP to it's SFS entry

    Returns
    -------
    res :  np.ndarray[L x 3]
        Pr(G_i = j  | F, tau)
    """

    if res is None:
        res = np.empty((SNP2SFS.shape[0], 3))

    # F = np.array(F0) if len(F0) == 1 else F0[SNP2SFS]
    # tau = np.array(tau0) if len(tau0) == 1 else tau0[SNP2SFS]
    F = F0[SNP2SFS]
    tau = tau0[SNP2SFS]

    res[:, 0] = F * (1 - tau) + (1 - F) * (1 - tau) ** 2  # HOMO REF
    res[:, 2] = F * tau + (1 - F) * tau * tau  # HOMO ALT
    res[:, 1] = (1 - F) * 2 * tau * (1 - tau)

    # assert np.allclose(np.sum(res, 1), 1)
    # np.clip(res, 0, 1, out=res)  # rounding

    return res


@njit(cache=True)
def p_gt_haploid(tau, SNP2SFS, res=None):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l}) for haploid genotypes"""
    F0 = np.ones_like(tau)
    res = p_gt_diploid(tau, F0, SNP2SFS=SNP2SFS, res=res)
    return res


def fwd_p_g(data, pars):
    """calculate forward probabilities of genotypes"""
    fwd_g = p_gt_diploid(pars.tau, pars.F, SNP2SFS=data.SNP2SFS)  # size [L x 3]
    if data.haploid_snps is not None:
        fwd_g[data.haploid_snps] = p_gt_haploid(
            pars.tau,
            SNP2SFS=data.SNP2SFS[data.haploid_snps],
        )
    return fwd_g


@njit(cache=True)
def bwd_p_o_given_x(R, F, e, b):
    """Caluclate Pr(O | X, e, b)

    Parameters
    ----------
    R : READS; np.ndarray[R x 1] int
    F : Flipped, np.ndarray[R x 1] bool: whether read is flipped
    e : float
        error rate from reference to alt allele
    b : float
        error rate from alt to reference allele

    Returns
    -------
    res : np.ndarray [R x 2]
        Pr(O_i | X_i = j)
    """

    res = np.empty((R.shape[0], 2))
    res[R == 0 & ~F, 0] = 1 - e  # ref=anc, obs=ref, X=anc
    res[R == 1 & ~F, 0] = e  # ref=anc, obs=alt, X=anc, error away from ref
    res[R == 0 & ~F, 1] = b  # ref=anc, obs=ref, X=der, error towards ref
    res[R == 1 & ~F, 1] = 1 - b  # ref=anc, obs=alt, X=der

    res[R == 0 & F, 0] = b  # ref=der, obs=ref, X=anc, error towards ref
    res[R == 1 & F, 0] = 1 - b  # ref=der, obs=ref, X=anc
    res[R == 0 & F, 1] = 1 - e  # ref=der, obs=ref, X=anc
    res[R == 1 & F, 1] = e  # ref=der, obs=ref, X=anc, error away from ref
    return res


@njit(cache=True)
def bwd_p_one_o_given_g(bwd_x, fwd_a, fwd_c, READ2SNP, READ2RG, n_reads):
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

    # A is matrix with contamination panel allele freq
    # A[,0] is Pr(X_i=0 |psi_i)
    # A[,1] is Pr(X_i=1 |psi_i) = psi_i
    fwd_a_mat = np.vstack((1 - fwd_a, fwd_a)).T[READ2SNP]
    p_o_cont = np.sum(fwd_a_mat * bwd_x, 1) * fwd_c[READ2RG]

    # [i, j, k] Pr(O_i= k | G_i = j )
    bwd_g = np.zeros((n_reads, 3))
    bwd_g[:] = np.expand_dims(p_o_cont, 1)  # cont is independent of genotype

    z = np.arange(3) / 2
    M = np.vstack((1 - z, z))
    bwd_g += (bwd_x @ M) * np.expand_dims((1 - fwd_c[READ2RG]), 1)

    return bwd_g


@njit(cache=True)
def bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps):
    """
    calculate  ℙ(O_l | G_l) = ∏_r ∏_j ℙ(O_lrj | G_l)

    where ℙ(O_lrj = 1 | G_l) = (c_r)psi_l * (1-c_r) G_l / 2.
    and   ℙ(O_lrj = 0 | G_l) = (c_r) (1 - psi_l) * (1-c_r) (1 - G_l / 2).
    """
    bwd_g = np.ones((n_snps, 3))

    for read, snp in enumerate(READ2SNP):
        bwd_g[snp] *= bwd_g1[read]

    return bwd_g


@njit(cache=True)
def fwd_p_x_cont(psi, READ2SNP):
    """fwd_p_x_cont calculate Pr(X | C=1, psi)

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


@njit(cache=True)
def fwd_p_x_nocont(fwd_g, READ2SNP):
    """Pr(X_i=j | C=0, G_l)"""
    M = np.array([[1.0, 0.5, 0.0], [0.0, 0.5, 1.0]]).T  # P(X=x | G=j)
    return fwd_g[READ2SNP] @ M


@njit(cache=True)
def fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, READ2RG):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    PC : array[n_obs]  of c_lr, the probability that the r'th read is
    contaminant
    PA : array[n_obs] the alt allele freq for observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = Σ_G Pr(X | G, C=0) Pr(G, C=0) + Σ_A Pr(X |A, C=1) Pr(A, C=1)
    """

    res = np.empty_like(fwd_x_nocont)
    res[:, 0] = (1 - fwd_c[READ2RG]) * fwd_x_nocont[:, 0] + fwd_c[READ2RG] * fwd_x_cont[
        :, 0
    ]
    res[:, 1] = (1 - fwd_c[READ2RG]) * fwd_x_nocont[:, 1] + fwd_c[READ2RG] * fwd_x_cont[
        :, 1
    ]

    return res


def posterior_g(bwd_g, fwd_g):
    """calculates marginalized posteriors for g, a and c

    post_g :   Pr(G_l | Z_l) x  (𝚷_rj Pr(O_lrj | G_l))    [L x 3]
    """

    post_g = bwd_g * fwd_g
    post_g /= np.expand_dims(np.sum(post_g, 1), 1)
    assert np.allclose(np.sum(post_g, 1), 1)

    return post_g


def posterior_x(bwd_x, *args, **kwargs):
    fwd_x = fwd_p_x(*args, **kwargs)
    post_x = fwd_x * bwd_x
    return post_x / np.sum(post_x, 1)[:, np.newaxis]


def full_posterior_genotypes(data, pars):
    fwd_g = fwd_p_g(data, pars)
    bwd_x = bwd_p_o_given_x(data.READS, data.FLIPPED_READS, e=pars.e, b=pars.b)
    bwd_g1 = bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.READ2SNP, data.READ2RG, data.n_reads
    )
    bwd_g = bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    return bwd_g, posterior_g(bwd_g, fwd_g)


@njit(cache=True)
def posterior_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, READ2RG):
    """calculate Pr(C | O) = Pr(C, O) / Pr (O)
    with Pr(C, O)  = Σ_X  Pr(O | X) Pr(X | C) Pr(C)

    return array [i, j k] with Pr(C_i = k | O_i = j)
    """
    x_cont = np.sum(fwd_x_cont * bwd_x, 1) * fwd_c[READ2RG]
    x_nocont = np.sum(fwd_x_nocont * bwd_x, 1) * (1 - fwd_c[READ2RG])

    post_c = x_cont / (x_cont + x_nocont)
    if np.any(np.isnan(post_c)):
        raise ValueError()

    return post_c


def calc_full_ll(data, pars):
    fwd_g = fwd_p_g(data, pars)
    bwd_x = bwd_p_o_given_x(data.READS, data.FLIPPED_READS, pars.e, pars.b)
    bwd_g1 = bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.READ2SNP, data.READ2RG, data.n_reads
    )
    bwd_g = bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)
    ll2 = calc_ll(fwd_g, bwd_g)
    return ll2


def calc_ll(fwd_g, bwd_g):
    ll = np.sum(np.log(np.sum(fwd_g * bwd_g, 1)))
    return ll
