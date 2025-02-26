import numpy as np
from numba import njit


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
    # res[:, 1] = 1 - res[:, 0] - res[:, 2]

    assert np.allclose(np.sum(res, 1), 1)
    res[:] = np.minimum(np.maximum(res, 0), 1)  # rounding

    return res


def slug_p_gt_haploid(*args, **kwargs):
    """calculate Pr(G_l | F_{Z_l}, tau_{Z_l})
    l is the number of SNPs

    probabilities are the same as diploid case with 100% inbreeding

    res: result, of size [l x 3]
    tau: tau parameter for expected derived allele frequency
    """

    return slug_p_gt_diploid(F=1, *args, **kwargs)


def slug_p_c(cont, res=None):
    """calculate Pr(C = 1)"""
    if res is None:
        res = np.empty((cont.shape[0]))

    res[:] = cont
    return res


def slug_p_obs_given_x(O, N, e, b, res=None):
    if res is None:
        res = np.empty((O.shape[0]), 2)

    res[:, 0] = binom_pmf(O, N, e)
    res[:, 1] = binom_pmf(O, N, 1 - b)

    return res


def slug_p_obs1_given_x(data, e, b, res=None):
    """calculate Pr(O = 1 | X)

    the probability of the observed read O given the molecule state X
    should take into account technical issues like sequencing error,
    reference bias and deamination; but no covariates

    e: error rate  ref -> alt
    b: error rate  alt -> ref
    data: vector of reads X

    """

    if res is None:
        res = np.empty((data.shape[0]))

    res[data == 0] = e  # Pr(O=1 | X = 0)
    res[data == 1] = 1 - b  # Pr(O=1 | X = 1)
    return res


def slug_p_x_given_g_nocont(n_snps, res=None):
    """calculate Pr(X=1 | C=0, G)

    the probability that the molecule carries a derived allele given genotype,
    assuming no contamination
    """
    if res is None:
        res = np.empty(n_snps, 3)

    res[:] = G / 2.0
    return res


def slug_p_x_given_cag(C, A, G, res=None):
    """calculate Pr(X=1 | G, A, C)

    the probability that the molecule carries a derived allele given genotype,
    contamination and contamination rate
    """
    if res is None:
        res = np.empty((G.shape[0]))

    res[C == 0] = G[C == 0] / 2.0
    res[C == 1] = A[C == 1] / 2.0
    return res


def slug_p_a_given_psi(psi, res=None):
    """calculate Pr(A=1 | psi, delta=0)"""
    if res is None:
        res = np.empty((G.shape[0]))

    res[:] = psi
    return res


def test_slug_p_gt_diploid():
    tau0 = 0.4
    F = np.array([0, 0.1, 1, 1])
    tau = np.array([tau0])
    res = np.empty((4, 3))

    slug_p_gt_diploid(tau, F, res)

    pred0 = np.array([(1 - tau0) ** 2, 2 * tau0 * (1 - tau0), tau0**2])
    pred2 = np.array([(1 - tau0), 0, tau0])
    pred1 = F[1] * pred2 + (1 - F[1]) * pred0

    pred = np.vstack((pred0, pred1, pred2, pred2))

    assert np.allclose(res - pred, 0)

    return res, res - pred


def test_slug_p_gt_haploid():
    tau0 = np.arange(10) / 10.0
    res = np.empty((10, 3))

    slug_p_gt_haploid(tau=tau0, res=res)
    assert np.allclose(res[:, 1], 0)
    assert np.allclose(res[:, 2], tau0)
    assert np.allclose(res[:, 0], 1 - tau0)


def test_slug_fwd_p_x():
    pg = np.array([[0, 0.5, 0.5], [1, 0, 0]])
    pa = np.array([0.3, 1])
    pc = np.array([0.3, 1])

    class _IX:
        n_obs = 3
        N_rg = 2
        OBS2RG = [0, 0, 1]
        OBS2SNP = [0, 1, 1]

    IX = _IX()

    res = slug_fwd_p_x(pg, pa, pc, IX)

    assert res[0] == 0.3 * 0.3 + 0.75 * 0.7
    assert res[1] == 0.3
    assert res[2] == 1


def test_slug_post_x():
    pg = np.array([[0, 0.5, 0.5], [1, 0, 0]])
    pa = np.array([0.3, 1])
    pc = np.array([0.3, 1])

    class _IX:
        n_obs = 3
        N_rg = 2
        OBS2RG = [0, 0, 1]
        OBS2SNP = [0, 1, 1]

    IX = _IX()

    res = slug_fwd_p_x(pg, pa, pc, IX)


def slug_fwd_p_x(PX_CONT, PX_NOCONT, PC, IX, res=None):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    PC : array[n_obs]  of c_lr, the probability that the r'th read is
    contaminant
    PA : array[n_obs] the alt allele freq for observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = \sum_G Pr(X | G, C=0) Pr(G, C=0) + \sum_A Pr(X |A, C=1) Pr(A, C=1)
    """
    if res is None:
        res = np.empty(IX.n_obs)

    res[:] = (1 - PC[IX.OBS2RG]) * PX_NOCONT[IX.OBS2SNP]
    res += PC[IX.OBS2RG] * PX_CONT[IX.OBS2SNP]

    return res


def slug_fwd_p_x_nocont(PG, res=None):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PG : array[n_snps x 3]  of G_l, the genotype at locus l
    P(X | G, C = 0) = \sum_G Pr(X | G, C=0) Pr(G)
    """
    if res is None:
        res = np.empty(PG.shape[0])

    res[:] = np.sum(np.arange(3) / 2.0 * PG, 1)
    return res


def slug_fwd_p_x_cont(PA, res=None):
    """Forward probability Pr(X_{lrj} = 1 | G_l, c_lr, a_lr)

    PA : array[n_obs] the alt allele freq for contamination observation lrj

    P(X | G, C, A) = Pr(X | C=0, G) Pr(C=0) + Pr(X | C=1, A) Pr(C=1)
                   = \sum_G Pr(X | G, C=0) Pr(G, C=0) + \sum_A Pr(X |A, C=1) Pr(A, C=1)
    """
    if res is None:
        res = np.empty(PA.shape[0])

    res[:] = PA
    return res


def slug_bwd_p_o_given_x(O, N, e, b, res=None):
    """Caluclate Pr(O | N, e, b, X)

    return [R x 2]
    """
    if res is None:
        res = np.empty((np.sum(N), 2))

    r = 0
    for o, n in zip(O, N):
        for i in range(n - o):
            res[r] = 1 - e, b  # Pr(O=0 | X =0) , Pr(O=0 | X = 1)
            r += 1
        for i in range(o):
            res[r] = e, 1 - b  # Pr(O=1 | X =0) , Pr(O=1 | X = 1)
            r += 1

    return res


def slug_bwd_p_o_given_x_simple(e, b, res=None):
    """Caluclate Pr(O | X, e, b)

    return [2 x 2]
    entry[k, o] is Pr(O_{lrj} = o |X_{lr} = k)


    """
    if res is None:
        res = np.empty((2, 2))

    res[0, 0] = 1 - e  # Pr(O_lrj = 0 | X=0)
    res[0, 1] = e  # Pr(O_lrj = 1 | X=0)
    res[1, 0] = b  # Pr(O_lrj = 0 | X=1)
    res[1, 1] = 1 - b  # Pr(O_lrj = 1 | X=1)

    return res


def slug_bwd_p_o_given_gac_simple(bx):
    """calculate Pr(O_lrj | G_l, A_l, C_lr)
        = Pr(O|X=0) Pr(X=0 | G, A, C) + Pr(O|X=1)Pr(X=1, G, A, C)
        = (1-BX) Pr(X=0 | GAC) + B Pr(X=1 | GAC)

    return P[o, g] = Pr(O=o | G=g)

    """
    res = np.empty((3, 2))
    res[0, 0] = (
        bx[0, 0] * 1.0 + bx[1, 0] * 0.0
    )  # Pr(O = 0 | G=0, C=0) = Pr(O=0 | A=0, C=1)
    res[0, 1] = (
        bx[0, 1] * 1.0 + bx[1, 1] * 0.0
    )  # Pr(O = 1 | G=0, C=0) = Pr(O=1 | A=0, C=1)

    res[1, 0] = bx[0, 0] * 0.5 + bx[1, 0] * 0.5  # Pr(O = 0 | G=1, C=0)
    res[1, 1] = bx[0, 1] * 0.5 + bx[1, 1] * 0.5  # Pr(O = 1 | G=1, C=0)

    res[2, 0] = (
        bx[0, 0] * 0.0 + bx[1, 0] * 1.0
    )  # Pr(O = 0 | G=2, C=0) = Pr(O=0 | A=1, C=1)
    res[2, 1] = (
        bx[0, 1] * 0.0 + bx[1, 1] * 1.0
    )  # Pr(O = 1 | G=2, C=0) = Pr(O=1 | A=1, C=1)

    return res.T


def slug_bwd_p_o_given_gac(BX, IX, res):
    """calculate Pr(O_lrj | G_l, A_l, C_lr)
        = Pr(O|X=0) Pr(X=0 | G, A, C) + Pr(O|X=1)Pr(X=1, G, A, C)
        = (1-BX) Pr(X=0 | GAC) + B Pr(X=1 | GAC)

    return [R x 5] array

    BX [R x 2] array of Pr(O | X = 0, 1)
    """
    res[0] = BX[r, 0] * 1.0 + BX[r, 1] * 0.0
    res[1] = BX[r, 0] * 0.5 + BX[r, 1] * 0.5
    res[2] = BX[r, 0] * 0.0 + BX[r, 1] * 1.0
    res[3] = res[o, 0]
    res[4] = res[o, 2]

    return res


def slug_post_x(FX, BX):
    """calculate Pr(X | 0)

    given Pr(X | Z, psi), Pr(O | X) and N, O


    PX[i, j, k] = Pr(X_i = k | O_i = j)
    """
    post_x = np.empty((FX.shape[0], 2, 2))

    # Pr(X_i = 0 , O_i = 0) = Pr(X_i = 0) Pr(O_i | X_i)
    post_x[:, 0, 0] = (1 - FX) * BX[0, 0]
    # Pr(X_i = 1 , O_i = 0) = Pr(X_i = 1) Pr(O_i | X_i)
    post_x[:, 0, 1] = FX * BX[1, 0]

    # Pr(X_i = 0 , O_i = 1) = Pr(X_i = 0) Pr(O_i | X_i)
    post_x[:, 1, 0] = (1 - FX) * BX[0, 1]
    # Pr(X_i = 1 , O_i = 1) = Pr(X_i = 1) Pr(O_i | X_i)
    post_x[:, 1, 1] = FX * BX[1, 1]

    post_x = post_x / np.expand_dims(np.sum(post_x, 2), 2)

    return post_x


def test_post_gac():
    fwd_g = np.array([[1, 0.0, 0.0], [1.0, 0.0, 0], [0.25, 0.5, 0.25]])
    fwd_a = np.array([1.0, 1, 0.5])
    fwd_c = np.array([0.1, 0.4, 0.1, 0.6662])
    O = np.array([1, 0, 3, 2, 3])
    N = np.array([1, 1, 5, 2, 6])
    e, b = 0.0, 0.00

    class _IX:
        n_obs = 5
        n_rg = 4
        n_snps = 3
        OBS2RG = np.array([0, 0, 1, 2, 3])
        OBS2SNP = np.array([0, 0, 1, 1, 2])

    IX = _IX()

    bwd_x = slug_bwd_p_o_given_x_simple(e, b)
    bwd_gac = slug_bwd_p_o_given_gac_simple(bwd_x)

    return c


def slug_bwd_p_o_given_g(bwd_gac, fwd_a, fwd_c, OBS2SNP, OBS2RG, n_obs):
    # [i, j, k] Pr(G_i=j | O = k )
    bwd_g = np.empty((n_obs, 3, 2))

    # contamination stuff, independent of g
    bwd_g[:, :, 0] = np.expand_dims(
        fwd_c[OBS2RG] * fwd_a[OBS2SNP] * bwd_gac[0, 2], 1
    )  # alt cont
    bwd_g[:, :, 0] += np.expand_dims(
        fwd_c[OBS2RG] * (1 - fwd_a[OBS2SNP]) * bwd_gac[0, 1], 1
    )  # ref cont
    bwd_g[:, :, 1] = np.expand_dims(
        fwd_c[OBS2RG] * fwd_a[OBS2SNP] * bwd_gac[1, 2], 1
    )  # alt cont
    bwd_g[:, :, 1] += np.expand_dims(
        fwd_c[OBS2RG] * (1 - fwd_a[OBS2SNP]) * bwd_gac[1, 1], 1
    )  # ref cont

    # Pr(O = 0 | G=0)
    bwd_g[:, 0, 0] += (1 - fwd_c[OBS2RG]) * bwd_gac[0, 0]  # ref endo

    # Pr(O = 1 | G=0)
    bwd_g[:, 0, 1] += (1 - fwd_c[OBS2RG]) * bwd_gac[1, 0]  # ref endo

    # Pr(O = 0 | G=1)
    bwd_g[:, 1, 0] += (1 - fwd_c[OBS2RG]) * bwd_gac[0, 1]  # ref endo

    # Pr(O = 1 | G=1)
    bwd_g[:, 1, 1] += (1 - fwd_c[OBS2RG]) * bwd_gac[1, 1]  # ref endo

    # Pr(O = 0 | G=2)
    bwd_g[:, 2, 0] += (1 - fwd_c[OBS2RG]) * bwd_gac[0, 2]  # ref endo

    # Pr(O = 1 | G=2)
    bwd_g[:, 2, 1] += (1 - fwd_c[OBS2RG]) * bwd_gac[1, 2]  # ref endo
    return bwd_g


@njit
def slug_post_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, OBS2SNP, OBS2RG, n_obs):
    # calculate for single read [i x j x k] Pr(C_i = j , O = k)
    post_c = np.empty((n_obs, 2, 2))
    # Pr(C_lrj=0 , O_lrj=0)
    post_c[:, 0, 0] = bwd_x[0, 0] * (1 - fwd_x_nocont[OBS2SNP]) * (1 - fwd_c[OBS2RG])
    post_c[:, 0, 0] += bwd_x[0, 1] * fwd_x_nocont[OBS2SNP] * (1 - fwd_c[OBS2RG])
    # Pr(C_lrj=0 , O_lrj=1)
    post_c[:, 0, 1] = bwd_x[0, 1] * (1 - fwd_x_nocont[OBS2SNP]) * (1 - fwd_c[OBS2RG])
    post_c[:, 0, 1] += bwd_x[1, 1] * fwd_x_nocont[OBS2SNP] * (1 - fwd_c[OBS2RG])
    # Pr(C_lrj=1 , O_lrj=0)
    post_c[:, 1, 0] = bwd_x[0, 0] * (1 - fwd_x_cont[OBS2SNP]) * fwd_c[OBS2RG]
    post_c[:, 1, 0] += bwd_x[1, 0] * fwd_x_cont[OBS2SNP] * fwd_c[OBS2RG]

    # Pr(C_lrj=1 , O_lrj=1)
    post_c[:, 1, 1] = bwd_x[0, 1] * (1 - fwd_x_cont[OBS2SNP]) * fwd_c[OBS2RG]
    post_c[:, 1, 1] += bwd_x[1, 1] * fwd_x_cont[OBS2SNP] * fwd_c[OBS2RG]

    # Pr(C | O)
    post_c = post_c / np.expand_dims(np.sum(post_c, 1), 1)
    na_to_zero(post_c)

    return post_c


def slug_post_g(bwd_g, fwd_g, N, O, OBS2SNP, n_snps):
    """calculates marginalized posteriors for g, a and c
    bwd_gac [2 x 3] array of Pr(O = i | G = j, c=0) which is the same
        as for c=1 with het removed
    fwd_g : forward probs of genotype [n_snps x 3]

    calc bwd_g: Pr(O | G) = Pr(O | G, C=0) Pr(C=0) +
        Pr(O | C=1, A=0) Pr(C=1) Pr(A=0)
        Pr(O | C=1, A=1) Pr(C=1) Pr(A=1)

    as events are disjoint

    return:
    post_g : posterior \prod_rj Pr(O_lrj | G_l) Pr(G_l | Z_l)   [L x 3]

    """

    # as we might have many reads per SNP, we do this in log scale
    log_post_g = np.zeros((n_snps, 3))
    for obs, snp in enumerate(OBS2SNP):
        print(obs, snp)
        log_post_g[snp] += np.log(bwd_g[obs, :, 0] + 1e-10) * (N[obs] - O[obs])
        log_post_g[snp] += np.log(bwd_g[obs, :, 1] + 1e-10) * O[obs]

    # Pr(O, G)
    log_post_g += np.log(fwd_g)

    g_scaling = np.max(log_post_g, 1)

    log_post_g = np.exp((log_post_g.T - g_scaling.T).T)

    # Pr(G | O)
    log_post_g /= np.expand_dims(np.sum(log_post_g, 1), 1)

    return log_post_g


def update_ftau(F, tau, post_g, IX):
    """updates the parameters f and tau
    F, tau [n_sfs x 1] : parameter arrays to be updated
    bwd_gac [2 x 3] array of Pr(O = i | G = j, c=0) which is the same
        as for c=1 with het removed
    fwd_g : forward probs of genotype [n_snps x 3]
    fwd_c : forward probs of contamination [n_obs x 1]

    calc bwd_g: Pr(O | G) = Pr(O | G, C=0) Pr(C=0) +
        Pr(O | C=1, A=0) Pr(C=1, A=0)
        Pr(O | C=1, A=1) Pr(C=1, A=1)

    as events are disjoint


    """

    for k in range(IX.n_sfs):
        g0, g1, g2 = np.mean(post_g[IX.SNP2SFS == k], 0)

        F[k] = 2 * g0 / (2 * g0 + g1 + 1e-8) - g1 / (2 * g2 + g1 + 1e-8)

        if np.isnan(F[k]):
            breakpoint()

        F[:] = np.minimum(1, np.maximum(0, F))

        tau[k] = g1 / 2.0 + g2


# GOOD
@njit
def na_to_zero(x):
    y = x.reshape(np.prod(np.array((x.shape))))
    y[np.isnan(y)] = 0.0


@njit
def update_c(c, post_c, O, N, OBS2SNP, OBS2RG, n_obs, n_rg):
    """update c
    parameters:
    c - vector of contamination rates to updates
    bwd_x : Pr(O | X)    [2 x 2] array
    fwd_x_cont \sum_A Pr(X | C=1, A)Pr(A)   [L x 1] array
    fwd_x_nocont \sum_G Pr(X | C=1, G) Pr(G) [L x 1] array
    O, N observed alt and tot reads; [O x 1]
    """

    # np.nan_to_num(post_c, copy= False)
    q = np.vstack((N - O, O)).T

    post_c *= np.expand_dims(q, 1)

    for r in range(n_rg):
        x0 = np.sum(post_c[OBS2RG == r, 0])
        x1 = np.sum(post_c[OBS2RG == r, 1])
        c[r] = x1 / (x0 + x1)
        # print(r, c[r], np.sum(OBS2RG==r))

    return c


def test_update_c():
    pg = np.array([[1, 0.0, 0.0], [1.0, 0.0, 0], [0.25, 0.5, 0.25]])
    pa = np.array([1.0, 1, 0.5])
    O = np.array([1, 0, 3, 2, 3])
    N = np.array([1, 1, 5, 2, 6])
    e, b = 0.0, 0.00

    pc = np.array([0.1, 0.4, 0.1, 0.6662])

    class _IX:
        n_obs = 5
        n_rg = 4
        n_snps = 3
        OBS2RG = np.array([0, 0, 1, 2, 3])
        OBS2SNP = np.array([0, 0, 1, 1, 2])

    IX = _IX()

    bwd_x = slug_bwd_p_o_given_x_simple(e, b)

    fwd_x_nocont = slug_fwd_p_x_nocont(pg)
    fwd_x_cont = slug_fwd_p_x_cont(pa)

    c = update_c(
        pc,
        bwd_x,
        fwd_x_cont,
        fwd_x_nocont,
        O,
        N,
        IX.OBS2SNP,
        IX.OBS2RG,
        IX.n_obs,
        IX.n_rg,
    )

    assert c[3] == pc[3]  # uninformative data
    assert c[0] == 0.5  # 50/50 strict informative
    assert c[2] == 1  # data only compatible with contamination
    assert c[1] == 3 / 5  # data only compatible with contamination

    return c


class SlugPars(object):
    def __init__(self):
        self.cont = np.array((0, 0.1, 1, 0.4))
        self.e = 0
        self.b = 0.1
        self.tau = np.array((1, 0.1))
        self.F = np.array((0, 1))


class SlugData(object):
    def __init__(self):
        self.N = np.array((1, 1, 10, 3, 4))
        self.O = np.array((0, 1, 4, 2, 4))
        self.psi = np.array([0, 0, 0.1])


class SlugIndex(object):
    def __init__(self):
        self.n_obs = 5
        self.n_rg = 4
        self.n_snps = 3
        self.n_sfs = 2
        self.OBS2RG = np.array([0, 0, 1, 2, 3])
        self.OBS2SNP = np.array([0, 0, 1, 1, 2])
        self.SNP2SFS = np.array([0, 0, 1])


def test_updates():
    data = SlugData()
    pars = SlugPars()
    IX = SlugIndex()


def update_pars(pars, data):
    """estimate tau, F, e, b, c; m-step of em algorithm"""

    """update tau, F
    need fwd_g, bwd_g,  IX.n_sfs, IX.SNP2SFS

    updates F, tau
    """

    fwd_g = slug_p_gt_diploid(pars.tau, pars.F)[IX.SNP2SFS]  # size [L x 3]
    fwd_a = data.psi  # size [L x 1]
    fwd_c = pars.cont  # size [O x 1]

    bwd_x = slug_bwd_p_o_given_x_simple(pars.e, pars.b)  # size [2 x 2]

    # est e, b here
    if do_update_eb:
        fwd_x_cont = slug_fwd_p_x_cont(fwd_a)
        fwd_x_nocont = slug_fwd_p_x_nocont(fwd_g)
        fwd_x = slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, IX)
        post_x = slug_post_x(fwd_x, bwd_x)
        tmp = np.sum(post_x, 0)
        pars.e = tmp[0, 1] / np.sum(tmp[0])
        pars.b = tmp[1, 0] / np.sum(tmp[1])

    if do_update_cont:
        post_c = slug_post_c(
            bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, IX.OBS2SNP, IX.OBS2RG, IX.n_obs
        )
        update_c(
            pars.cont, post_c, data.O, data.N, IX.OBS2SNP, IX.OBS2RG, IX.n_obs, IX.n_rg
        )

    if do_update_ftau:
        bwd_gac = slug_bwd_p_o_given_gac_simple(bwd_x)  # size [3 x 2]
        bwd_g = slug_bwd_p_o_given_g(
            bwd_gac, fwd_a, fwd_c, IX.OBS2SNP, IX.OBS2RG, IX.n_obs
        )  # size [O x 3]
        post_g = slug_post_g(bwd_g, fwd_g, data.N, data.O, IX.OBS2SNP, IX.n_snps)
        update_ftau(pars.F, pars.tau, post_g, IX)

    """update e, b
    need O_flat, fwd_x, bwd_x
    updates e, b
    """

    # post_x is Pr(X = 1 | \theta')
    post_x = fwd_x * bwd_x / (fwd_x * bwd_x + (1 - fwd_x) * (1 - bwd_x))
    e = np.mean((1 - post_x)[O_flat == 1])  # Pr(O= 1 | X = 0)
    b = np.mean(post_x[O_flat == 0])  # Pr(O= 0 | X = 1)
