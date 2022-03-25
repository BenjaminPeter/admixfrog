import numpy as np
from numba import njit
from math import lgamma, exp, pow
from copy import deepcopy


@njit
def binom_pmf(ALT, REF, p):
    res = pow(p, ALT) * pow(1.0 - p, REF)
    # if ALT > 0 and REF > 1:
    #        res *= exp(lgamma(ALT + REF + 1) - lgamma(REF + 1) - lgamma(ALT + 1))
    return res


def slug_p_gt_diploid(tau, F, res=None):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l})
    l is the number of SNPs

    res: result, of size [l x 3]
    F: F-parameter measuring drift, either scalar or vector of length l
    tau: tau parameter for expected derived allele frequency
    """

    if res is None:
        res = np.empty((tau.shape[0], 3))

    res[:, 0] = F * (1 - tau) + (1 - F) * (1 - tau) * (1 - tau)  # HOMO REF
    res[:, 2] = F * tau + (1 - F) * tau * tau  # HOMO ALT
    res[:, 1] = (1 - F) * 2 * tau * (1 - tau)

    assert np.allclose(np.sum(res, 1), 1)
    np.clip(res, 0, 1, out=res)  # rounding

    return res


def slug_p_gt_haploid(*args, **kwargs):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l})
    l is the number of SNPs

    probabilities are the same as diploid case with 100% inbreeding

    res: result, of size [l x 3]
    tau: tau parameter for expected derived allele frequency
    """

    return slug_p_gt_diploid(F=1, *args, **kwargs)


@njit
def slug_fwd_p_o(fwd_x, REF, ALT, e, b):
    p = fwd_x * (1 - b) + (1 - fwd_x) * e
    fwd_o = np.empty_like(p)
    for i, (p_, a, r) in enumerate(zip(p, ALT, REF)):
        fwd_o[i] = binom_pmf(a, r, p_)

    return fwd_o


def slug_ll_o(slug_fwd_o):
    return np.sum(np.log(fwd_o))


def slug_fwd_p_x(PX_CONT, PX_NOCONT, PC, IX, res=None):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    PC : array[n_obs]  of c_lr, the probability that the r'th read is
    contaminant
    PA : array[n_obs] the alt allele freq for observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = Σ_G Pr(X | G, C=0) Pr(G, C=0) + Σ_A Pr(X |A, C=1) Pr(A, C=1)
    """
    if res is None:
        res = np.empty(IX.n_obs)

    res[:] = (1 - PC[IX.OBS2RG]) * PX_NOCONT[IX.OBS2SNP]
    res += PC[IX.OBS2RG] * PX_CONT[IX.OBS2SNP]

    return res


def reads_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, data):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    PC : array[n_obs]  of c_lr, the probability that the r'th read is
    contaminant
    PA : array[n_obs] the alt allele freq for observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = Σ_G Pr(X | G, C=0) Pr(G, C=0) + Σ_A Pr(X |A, C=1) Pr(A, C=1)
    """

    res = (1 - fwd_c[data.READ2RG]) * fwd_x_nocont[data.READ2SNP] + fwd_c[
        data.READ2RG
    ] * fwd_x_cont[data.READ2SNP]

    return res


def slug_fwd_p_x_nocont(PG, res=None):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    P(X | G, C = 0) = Σ_G Pr(X | G, C=0) Pr(G)
    """
    if res is None:
        res = np.empty(PG.shape[0])

    res[:] = np.sum(np.arange(3) / 2.0 * PG, 1)
    return res


def slug_fwd_p_x_cont(PA, res=None):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PA : array[n_obs] the alt allele freq for contamination observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = Σ_G Pr(X | G, C=0) Pr(G, C=0) + Σ_A Pr(X |A, C=1) Pr(A, C=1)
    """
    if res is None:
        res = np.empty(PA.shape[0])

    res[:] = PA
    return res


def slug_bwd_p_o_given_x(e, b, res=None):
    """Caluclate Pr(O | X, e, b)

    return [2 x 2]
    entry[i, j] is Pr(O_{lr} = j |X_{lr} = o)


    """
    if res is None:
        res = np.empty((2, 2))

    res[0] = 1 - e, e  # Pr(O_lrj = 0, 1 | X=0)
    res[1] = b, 1 - b  # Pr(O_lrj = 0, 1 | X=1)
    return res


def reads_bwd_p_o_given_x(READS, e, b):
    """Caluclate Pr(O | X, e, b)

    entry[i, j] is Pr(O_i |X_i = j)

    """
    res = np.empty((READS.shape[0], 2))

    res[READS == 0, 0] = 1 - e
    res[READS == 0, 1] = e
    res[READS == 1, 0] = b
    res[READS == 1, 1] = 1 - b
    return res


def slug_bwd_p_o_given_gac(bx):
    """calculate Pr(O_lrj | G_l, A_l, C_lr)
        = Pr(O|X=0) Pr(X=0 | G, A, C) + Pr(O|X=1)Pr(X=1, G, A, C)
        = (1-BX) Pr(X=0 | GAC) + B Pr(X=1 | GAC)

    return P[g, o] = Pr(O=o | G=g)

    """
    v = np.vstack(((1, 0.5, 0), (0, 0.5, 1))).T  # Pr(X | G)
    return v @ bx

    res = np.empty((3, 2))
    res[0] = bx[0] * 1.0 + bx[1] * 0.0
    res[1] = bx[0] * 0.5 + bx[1] * 0.5
    res[2] = (
        bx[0] * 0.0 + bx[1] * 1.0
    )  # Pr(O | G= 2) = Pr(O | X = 0) Pr(X=0 | G=2) + Pr(O | X=1) Pr(X=1 | G=2)

    return res


def slug_bwd_p_o_given_g_nocont(bx):
    return slug_bwd_p_o_given_gac(bx)


def slug_bwd_p_o_given_a_cont(bx):
    return slug_bwd_p_o_given_gac(bx)[[0, 2]]


def read_post_x(fwd_x, bwd_x):
    post_x = np.vstack(((1 - fwd_x) * bwd_x[:, 0], (fwd_x) * bwd_x[:, 1])).T
    return post_x / np.expand_dims(np.sum(post_x, 1), 1)


def slug_post_x(FX, BX):
    """calculate Pr(X | 0)

    given Pr(X | Z, psi) and Pr(O | X)


    calculate PX[i, j, k] = Pr(X_i = k | O_i = j)
    """
    post_x = np.empty((FX.shape[0], 2, 2))

    # Pr(X_i = 0 , O_i = 0) = Pr(X_i = 0) Pr(O_i = 0 | X_i = 0)
    post_x[:, 0, 0] = (1 - FX) * BX[0, 0]
    # Pr(X_i = 1 , O_i = 0) = Pr(X_i = 1) Pr(O_i = 0| X_i = 1)
    post_x[:, 0, 1] = FX * BX[1, 0]

    # Pr(X_i = 0 , O_i = 1) = Pr(X_i = 0) Pr(O_i = 1| X_i = 0)
    post_x[:, 1, 0] = (1 - FX) * BX[0, 1]
    # Pr(X_i = 1 , O_i = 1) = Pr(X_i = 1) Pr(O_i = 1| X_i = 1)
    post_x[:, 1, 1] = FX * BX[1, 1]

    # Pr(X | O)
    post_x = post_x / np.expand_dims(np.sum(post_x, 2), 2)
    assert np.allclose(np.sum(post_x, 2), 1)

    return post_x


@njit
def slug_bwd_p_o_given_a(bwd_gac, fwd_g, fwd_c, OBS2SNP, OBS2RG, n_obs):
    # [i, j, k] Pr(O_i= k | A_i = j )
    bwd_a = np.empty((n_obs, 2, 2))

    # contamination stuff, independent of g
    # Pr(A_i | O = 0, C=0) = Pr(C=0) Pr(G=2) Pr(O =0 | G =2, C=0)
    bwd_a[:, :, 0] = np.expand_dims(
        (1 - fwd_c[OBS2RG]) * fwd_g[OBS2SNP, 2] * bwd_gac[2, 0], 1
    )  # alt cont
    # Pr(A_i | O = 0, C=0) = Pr(C=0) Pr(G=1) Pr(O =0 | G =1, C=0)
    bwd_a[:, :, 0] = np.expand_dims(
        (1 - fwd_c[OBS2RG]) * fwd_g[OBS2SNP, 1] * bwd_gac[1, 0], 1
    )  # alt cont
    bwd_a[:, :, 0] = np.expand_dims(
        (1 - fwd_c[OBS2RG]) * fwd_g[OBS2SNP, 0] * bwd_gac[0, 0], 1
    )  # alt cont
    # Pr(A_i | O = 1, C=0) += Pr(C=1) Pr(G=0) Pr(O =1 | G =0, C=1)
    bwd_a[:, :, 1] = np.expand_dims(
        (1 - fwd_c[OBS2RG]) * fwd_g[OBS2SNP, 2] * bwd_gac[2, 1], 1
    )  # alt cont
    bwd_a[:, :, 1] = np.expand_dims(
        (1 - fwd_c[OBS2RG]) * fwd_g[OBS2SNP, 1] * bwd_gac[1, 1], 1
    )  # alt cont
    bwd_a[:, :, 1] = np.expand_dims(
        (1 - fwd_c[OBS2RG]) * fwd_g[OBS2SNP, 0] * bwd_gac[0, 1], 1
    )  # alt cont

    # Pr(O = 0 | A=0, C=1)
    bwd_a[:, 0, 0] += fwd_c[OBS2RG] * bwd_gac[0, 0]

    # Pr(O = 1 | A=0, C=1)
    bwd_a[:, 0, 1] += fwd_c[OBS2RG] * bwd_gac[0, 1]

    # Pr(O = 0 | A=1, C=1)
    bwd_a[:, 1, 0] += fwd_c[OBS2RG] * bwd_gac[2, 0]

    # Pr(O = 1 | A=1, C=1)
    bwd_a[:, 1, 1] += fwd_c[OBS2RG] * bwd_gac[2, 1]

    return bwd_a


@njit
def reads_bwd_p_one_o_given_g(bwd_x, fwd_a, fwd_c, READ2SNP, READ2RG, n_reads):
    """Pr(O | G ) = \sum Pr(O , C=0)  + Pr(O | C=1, G) Pr(C=1) as C=0 implies G is independent of O"""

    fwd_a_mat = np.vstack((1 - fwd_a, fwd_a)).T[READ2SNP]
    p_o_cont = np.sum(fwd_a_mat * bwd_x, 1) * fwd_c[READ2RG]

    # [i, j, k] Pr(O_i= k | G_i = j )
    bwd_g = np.zeros((n_reads, 3))
    bwd_g[:] = np.expand_dims(p_o_cont, 1)  # cont is independent of genotype

    z = np.arange(3) / 2
    M = np.vstack((z, 1 - z))
    bwd_g += (bwd_x @ M) * np.expand_dims((1 - fwd_c[READ2RG]), 1)

    return bwd_g


@njit
def slug_bwd_p_one_o_given_g(bwd_x, fwd_a, fwd_c, OBS2SNP, OBS2RG, n_obs):
    """Pr(O | G ) = \sum Pr(O , C=0)  + Pr(O | C=1, G) Pr(C=1) as C=0 implies G is independent of O"""

    bwd_g_nocont = slug_bwd_p_o_given_g_nocont(bwd_x)
    bwd_a_cont = slug_bwd_p_o_given_a_cont(bwd_x)
    fwd_a_mat = np.vstack((1 - fwd_a, fwd_a)).T[OBS2SNP]

    # [i, j, k] Pr(O_i= k | G_i = j )
    bwd_g = np.zeros((n_obs, 3, 2))

    p_cont = fwd_a_mat @ bwd_a_cont  # Pr(O_i = j  | C_i = 1) * Pr(C_i = 1)
    bwd_g[:] = np.expand_dims(p_cont * np.expand_dims(fwd_c[OBS2RG], 1), 1)
    bwd_g += np.expand_dims(bwd_g_nocont, 0) * np.expand_dims(
        (1 - fwd_c[OBS2RG]), (1, 2)
    )

    assert np.allclose(np.sum(bwd_g, 2), 1)

    return bwd_g


@njit
def reads_bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps):
    """
    calculate  ℙ(O_l | G_l) = ∏_r ∏_j ℙ(O_lrj | G_l)

    where ℙ(O_lrj = 1 | G_l) = (c_r)psi_l * (1-c_r) G_l / 2.
    and   ℙ(O_lrj = 0 | G_l) = (c_r) (1 - psi_l) * (1-c_r) (1 - G_l / 2).
    """
    bwd_g = np.ones((n_snps, 3))

    for read, snp in enumerate(READ2SNP):
        bwd_g[snp] *= bwd_g1[read]

    return bwd_g


# @njit
def slug_bwd_p_all_o_given_g(bwd_g1, READS, READSNP, n_snps):
    """
    calculate  ℙ(O_l | G_l) = ∏_r ∏_j ℙ(O_lrj | G_l)

    where ℙ(O_lrj = 1 | G_l) = (c_r)psi_l * (1-c_r) G_l / 2.
    and   ℙ(O_lrj = 0 | G_l) = (c_r) (1 - psi_l) * (1-c_r) (1 - G_l / 2).
    """
    bwd_g = np.ones((n_snps, 3))

    for obs, snp in enumerate(OBS2SNP):
        for g in range(3):
            # bwd_g[snp, g] *= binom_pmf(ALT[obs], REF[obs], bwd_g1[obs, g, 1])
            bwd_g[snp, g] *= bwd_g1[obs, g, 0] ** REF[obs]
            bwd_g[snp, g] *= bwd_g1[obs, g, 1] ** ALT[obs]
    return bwd_g


@njit
def slug_bwd_p_all_o_given_g_scaling(bwd_g1, REF, ALT, OBS2SNP, n_snps):
    """
    calculate  ℙ(O_l | G_l) = ∏_r ∏_j ℙ(O_lrj | G_l)

    where ℙ(O_lrj = 1 | G_l) = (c_r)psi_l * (1-c_r) G_l / 2.
    and   ℙ(O_lrj = 0 | G_l) = (c_r) (1 - psi_l) * (1-c_r) (1 - G_l / 2).
    """
    bwd_g = np.ones((n_snps, 3))
    scaling = 0.0

    for obs, snp in enumerate(OBS2SNP):
        for g in range(3):
            # bwd_g[snp, g] *= binom_pmf(ALT[obs], REF[obs], bwd_g1[obs, g, 1])
            bwd_g[snp, g] *= bwd_g1[obs, g, 0] ** REF[obs]
            bwd_g[snp, g] *= bwd_g1[obs, g, 1] ** ALT[obs]
        s = np.max(bwd_g[snp])
        bwd_g[snp] /= s
        scaling += np.log(s)
    return bwd_g, scaling


def full_posterior_genotypes(data, pars):
    fwd_g = slug_fwd_p_g(data, pars)
    bwd_x = reads_bwd_p_o_given_x(data.READS, e=pars.e, b=pars.b)  # Pr(O_lrj | G_lrj)
    bwd_g1 = reads_bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.READ2SNP, data.READ2RG, data.n_reads
    )
    bwd_g = reads_bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    return slug_post_g(bwd_g, fwd_g)


def reads_full_posterior_genotypes(data, pars):
    fwd_g = slug_fwd_p_g(data, pars)
    bwd_x = slug_bwd_p_o_given_x(pars.e, pars.b)
    bwd_g1 = slug_bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.OBS2SNP, data.OBS2RG, data.n_obs
    )
    bwd_g = slug_bwd_p_all_o_given_g(
        bwd_g1, data.REF, data.ALT, data.OBS2SNP, data.n_snps
    )

    return slug_post_g(bwd_g, fwd_g)


def full_posterior_cont(data, pars):
    fwd_g = slug_fwd_p_g(data, pars)
    bwd_x = slug_bwd_p_o_given_x(pars.e, pars.b)
    fwd_x_nocont = slug_fwd_p_x_cont(data.psi)
    fwd_x_cont = slug_fwd_p_x_nocont(fwd_g)
    fwd_c = pars.cont
    bwd_g1 = slug_bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.OBS2SNP, data.OBS2RG, data.n_obs
    )
    bwd_g = slug_bwd_p_all_o_given_g(
        bwd_g1, data.REF, data.ALT, data.OBS2SNP, data.n_snps
    )

    return slug_post_c(
        bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, data.OBS2SNP, data.OBS2RG, data.n_obs
    )


def slug_post_g(bwd_g, fwd_g):
    """calculates marginalized posteriors for g, a and c
    bwd_gac [2 x 3] array of Pr(O = i | G = j, c=0) which is the same
        as for c=1 with het removed
    fwd_g : forward probs of genotype [n_snps x 3]


    as events are disjoint

    return:
    post_g :   Pr(G_l | Z_l) 𝚷_rj Pr(O_lrj | G_l)    [L x 3]

    """

    post_g = bwd_g * fwd_g
    post_g /= np.expand_dims(np.sum(post_g, 1), 1)
    assert np.allclose(np.sum(post_g, 1), 1)

    return post_g


# @njit
def slug_post_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, OBS2SNP, OBS2RG, n_obs):
    """calculate Pr(C | O) = Pr(C, O) / Pr (O)
    with Pr(C, O)  = Σ_X  Pr(O | X) Pr(X | C) Pr(C)

    return array [i, j k] with Pr(C_i = k | O_i = j)
    """
    x_cont_mat = np.vstack((1 - fwd_x_cont[OBS2SNP], fwd_x_cont[OBS2SNP])).T @ bwd_x
    x_nocont_mat = (
        np.vstack((1 - fwd_x_nocont[OBS2SNP], fwd_x_nocont[OBS2SNP])).T @ bwd_x
    )
    x_cont_mat *= np.expand_dims(fwd_c[OBS2RG], 1)
    x_nocont_mat *= np.expand_dims(1 - fwd_c[OBS2RG], 1)
    post_c = np.stack((x_nocont_mat, x_cont_mat), axis=2)
    post_c /= np.expand_dims(np.sum(post_c, 2), 2)

    return post_c

    # calculate for single read [i x j x k] Pr(C_i = j , O = k)
    post_c = np.empty((n_obs, 2, 2))

    # post_c[:] = np.expand_dims(bwd_x, 0)
    # Pr(C_i =0 , O_i=0)
    # Pr(C, O)       = Pr(O=0 | X=0)   *  Pr(X=0 | C = 0)   * Pr(C = 0)
    post_c[:, 0, 0] = bwd_x[0, 0] * (1 - fwd_x_nocont[OBS2SNP]) * (1 - fwd_c[OBS2RG])
    post_c[:, 0, 0] += bwd_x[1, 0] * fwd_x_nocont[OBS2SNP] * (1 - fwd_c[OBS2RG])
    post_c[:, 1, 0] = bwd_x[0, 1] * (1 - fwd_x_nocont[OBS2SNP]) * (1 - fwd_c[OBS2RG])
    post_c[:, 1, 0] += bwd_x[1, 1] * fwd_x_nocont[OBS2SNP] * (1 - fwd_c[OBS2RG])
    post_c[:, 0, 1] = bwd_x[0, 0] * (1 - fwd_x_cont[OBS2SNP]) * fwd_c[OBS2RG]
    post_c[:, 0, 1] += bwd_x[1, 0] * fwd_x_cont[OBS2SNP] * fwd_c[OBS2RG]
    post_c[:, 1, 1] = bwd_x[0, 1] * (1 - fwd_x_cont[OBS2SNP]) * fwd_c[OBS2RG]
    post_c[:, 1, 1] += bwd_x[1, 1] * fwd_x_cont[OBS2SNP] * fwd_c[OBS2RG]

    # Pr(C | O)
    post_c /= np.expand_dims(np.sum(post_c, 2), 2)
    return post_c


def slug_fwd_p_g(data, pars):
    fwd_g = slug_p_gt_diploid(pars.tau, pars.F)[data.SNP2SFS]  # size [L x 3]
    if data.haploid_snps is not None:
        fwd_g[data.haploid_snps] = slug_p_gt_haploid(pars.tau)[
            data.SNP2SFS[data.haploid_snps]
        ]  # size [L x 3]

    return fwd_g


def calc_full_ll(data, pars):
    fwd_a = data.psi  # size [L x 1]
    fwd_c = pars.cont
    fwd_g = slug_fwd_p_g(data, pars)
    fwd_x_cont = slug_fwd_p_x_cont(fwd_a)  # size [L x 1]
    fwd_x_nocont = slug_fwd_p_x_nocont(fwd_g)  # size [L x 1]
    fwd_x = slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, data)
    # fwd_o = slug_fwd_p_o(fwd_x, data.REF, data.ALT, pars.e, pars.b)

    bwd_x = slug_bwd_p_o_given_x(pars.e, pars.b)  # Pr(O_lrj | G_lrj)
    bwd_g1 = slug_bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.OBS2SNP, data.OBS2RG, data.n_obs
    )
    bwd_g = slug_bwd_p_all_o_given_g(
        bwd_g1, data.REF, data.ALT, data.OBS2SNP, data.n_snps
    )

    ll1 = np.sum(np.log(fwd_o))
    ll2 = calc_ll(fwd_g, bwd_g)
    print(f"{ll1:.4f} | {ll2:.4f} ")
    return ll2


def calc_full_ll_reads(data, pars):
    fwd_g = slug_fwd_p_g(data, pars)
    bwd_x = reads_bwd_p_o_given_x(data.READS, e=pars.e, b=pars.b)  # Pr(O_lrj | G_lrj)
    bwd_g1 = reads_bwd_p_one_o_given_g(
        bwd_x, data.psi, pars.cont, data.READ2SNP, data.READ2RG, data.n_reads
    )
    bwd_g = reads_bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    ll2 = calc_ll(fwd_g, bwd_g)
    return ll2


def calc_ll(fwd_g, bwd_g):
    ll = np.sum(np.log(np.sum(fwd_g * bwd_g, 1)))
    return ll
