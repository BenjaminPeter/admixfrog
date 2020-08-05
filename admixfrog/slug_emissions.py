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

    if not np.allclose(np.sum(res, 1), 1):
        breakpoint()
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

    res[0, 0] = 1 - e  # Pr(O_lrj = 0 | X=0)
    res[0, 1] = e  # Pr(O_lrj = 1 | X=0)
    res[1, 0] = b  # Pr(O_lrj = 0 | X=1)
    res[1, 1] = 1 - b  # Pr(O_lrj = 1 | X=1)

    return res


def slug_bwd_p_o_given_gac(bx):
    """calculate Pr(O_lrj | G_l, A_l, C_lr) 
            = Pr(O|X=0) Pr(X=0 | G, A, C) + Pr(O|X=1)Pr(X=1, G, A, C)
            = (1-BX) Pr(X=0 | GAC) + B Pr(X=1 | GAC)

        return P[g, o] = Pr(O=o | G=g)
        
    """
    res = np.empty((3, 2))
    res[0, 0] = (
        bx[0, 0] * 1.0 + bx[1, 0] * 0.0
    )  # Pr(O = 0 | G=0, C=0) = Pr(O=0 | A=0, C=1)
    res[0, 1] = (
        bx[0, 1] * 1.0 + bx[1, 1] * 0.0
    )  # Pr(O = 1 | G=0, C=0) = Pr(O=1 | A=0, C=1)

    # Pr(O=0 | G = 1) = Pr(O=0 | X=0) Pr(X=0 | G=1) + Pr(O=0 | X=1) Pr(X=1 | G=1)
    # res[1, 0]       = bx[0, 0]      Pr(X=0 | G=1) + bx[1, 0] Pr(X=1 | G=1)
    res[1, 0] = bx[0, 0] * 0.5 + bx[1, 0] * 0.5  # Pr(O = 0 | G=1, C=0)
    res[1, 1] = bx[0, 1] * 0.5 + bx[1, 1] * 0.5  # Pr(O = 1 | G=1, C=0)

    res[2, 0] = (
        bx[0, 0] * 0.0 + bx[1, 0] * 1.0
    )  # Pr(O = 0 | G=2, C=0) = Pr(O=0 | A=1, C=1)
    res[2, 1] = (
        bx[0, 1] * 0.0 + bx[1, 1] * 1.0
    )  # Pr(O = 1 | G=2, C=0) = Pr(O=1 | A=1, C=1)

    return res


def slug_post_x(FX, BX):
    """calculate Pr(X | 0) 

    given Pr(X | Z, psi), Pr(O | X) and N, O


    PX[i, j, k] = Pr(X_i = k | O_i = j)
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
    na_to_zero(post_x)

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
def slug_bwd_p_o_given_g_nocont(bwd_gac, fwd_a, fwd_c, OBS2SNP, OBS2RG, n_obs):
    # [i, j, k] Pr(O_i= k | G_i = j )
    bwd_g = np.zeros((n_obs, 3, 2))

    # contamination stuff, independent of g
    # Pr(G_i | O = 0, C=1) = Pr(C=1) Pr(A=1) Pr(O =0 | A =1, C=1)

    if False: #also marginalizing contamination
        bwd_g[:, :, 0] = np.expand_dims(
            fwd_c[OBS2RG] * fwd_a[OBS2SNP] * bwd_gac[2, 0], 1
        )  # alt cont
        # Pr(G_i | O = 0, C=1) += Pr(C=1) Pr(A=0) Pr(O =0 | A =0, C=1)
        bwd_g[:, :, 0] += np.expand_dims(
            fwd_c[OBS2RG] * (1 - fwd_a[OBS2SNP]) * bwd_gac[0, 0], 1
        )  # ref cont
        bwd_g[:, :, 1] = np.expand_dims(
            fwd_c[OBS2RG] * fwd_a[OBS2SNP] * bwd_gac[2, 1], 1
        )  # alt cont
        bwd_g[:, :, 1] += np.expand_dims(
            fwd_c[OBS2RG] * (1 - fwd_a[OBS2SNP]) * bwd_gac[0, 1], 1
        )  # ref cont

    # Pr(O = 0 | G=0)
    bwd_g[:, 0, 0] += (1 - fwd_c[OBS2RG]) * bwd_gac[0, 0]  # ref endo

    # Pr(O = 1 | G=0)
    bwd_g[:, 0, 1] += (1 - fwd_c[OBS2RG]) * bwd_gac[0, 1]  # ref endo

    # Pr(O = 0 | G=1)
    bwd_g[:, 1, 0] += (1 - fwd_c[OBS2RG]) * bwd_gac[1, 0]  # ref endo

    # Pr(O = 1 | G=1)
    bwd_g[:, 1, 1] += (1 - fwd_c[OBS2RG]) * bwd_gac[1, 1]  # ref endo

    # Pr(O = 0 | G=2)
    bwd_g[:, 2, 0] += (1 - fwd_c[OBS2RG]) * bwd_gac[2, 0]  # ref endo

    # Pr(O = 1 | G=2)
    bwd_g[:, 2, 1] += (1 - fwd_c[OBS2RG]) * bwd_gac[2, 1]  # ref endo
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


@njit
def slug_post_g(bwd_g, fwd_g, ALT, REF, OBS2SNP, n_snps):
    """calculates marginalized posteriors for g, a and c
        bwd_gac [2 x 3] array of Pr(O = i | G = j, c=0) which is the same
            as for c=1 with het removed
        fwd_g : forward probs of genotype [n_snps x 3]

        calc bwd_g: Pr(O | G) = Pr(O | G, C=0) Pr(C=0) + 
            Pr(O | C=1, A=0) Pr(C=1) Pr(A=0)
            Pr(O | C=1, A=1) Pr(C=1) Pr(A=1)

        as events are disjoint

        return:
        post_g : posterior 𝚷_rj Pr(O_lrj | G_l) Pr(G_l | Z_l)   [L x 3]

    """

    # as we might have many reads per SNP, we do this in log scale
    log_post_g = np.zeros((n_snps, 3))
    for obs, snp in enumerate(OBS2SNP):
        log_post_g[snp] += np.log(bwd_g[obs, :, 0] + 1e-100) * REF[obs]
        log_post_g[snp] += np.log(bwd_g[obs, :, 1] + 1e-100) * ALT[obs]

    # Pr(O, G)
    log_post_g += np.log(fwd_g + 1e-100)

    x = log_post_g
    for i in range(x.shape[0]):
        x[i] = np.exp(x[i] - np.max(x[i]))

    # g_scaling = np.max(log_post_g, 1)
    # log_post_g = np.exp((log_post_g.T - g_scaling.T).T)

    # Pr(G | O)
    log_post_g /= np.expand_dims(np.sum(log_post_g, 1), 1)

    return log_post_g


@njit
def slug_post_a(bwd_a, fwd_a, REF, ALT, OBS2SNP, n_snps):
    """calculates marginalized posteriors a and c
        fwd_g : forward probs of genotype [n_snps x 3]

        calc bwd_a: Pr(O | G) = Pr(O | A, C=1) Pr(C=1) + 
            Pr(O | C=0, G=0) Pr(C=0) Pr(G=0)
            Pr(O | C=0, G=1) Pr(C=0) Pr(G=1)
            Pr(O | C=0, G=2) Pr(C=0) Pr(G=2)

        as events are disjoint

        return:
        post_a : posterior 𝚷_rj Pr(O_lrj | A_l) Pr(A_l | delta, psi)   [L x 3]

    """

    # as we might have many reads per SNP, we do this in log scale
    # log Pr(A_i = k | O_i = j)
    log_post_a = np.zeros((n_snps, 2))

    # get log Pr(O | A)
    for obs, snp in enumerate(OBS2SNP):
        log_post_a[snp] += np.log(bwd_a[obs, :, 0] + 1e-100) * REF[obs]
        log_post_a[snp] += np.log(bwd_a[obs, :, 1] + 1e-100) * ALT[obs]

    # log Pr(O, A)
    log_post_a += np.log(fwd_a + 1e-100)

    x = log_post_a
    for i in range(x.shape[0]):
        x[i] = np.exp(x[i] - np.max(x[i]))

    # g_scaling = np.max(log_post_g, 1)
    # log_post_g = np.exp((log_post_g.T - g_scaling.T).T)

    # Pr(G | O)
    log_post_a /= np.expand_dims(np.sum(log_post_a, 1), 1)

    return log_post_a


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
        if IX.haploid_snps is None:
            h0, h1, h2 = 0, 0 ,0
        else:
            ix1 = IX.SNP2SFS == k
            ix1[IX.haploid_snps] = False
            if np.all(~ix1):
                h0, h1, h2 = 0, 0 ,0
            else:
                h0, h1, h2 = np.mean(post_g[ix1], 0)

        f0, f1, f2 = g0 - h0, g1 - h1, g2 - h2
        #F[k] = 2 * f0 / (2 * f0 + f1 + 1e-100) - f1 / (2 * f2 + f1 + 1e-100)

        #np.clip(F, 0, 1, out=F)  # round to [0, 1]

        tau[k] = g1 / 2.0 + g2

def update_ftau_debug(data, pars, post_g):
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
    ll = - np.inf

    for k in range(data.n_sfs):
        g0, g1, g2 = np.mean(post_g[data.SNP2SFS == k], 0)
        if data.haploid_snps is None:
            h0, h1, h2 = 0, 0 ,0
        else:
            ix1 = data.SNP2SFS == k
            ix1[data.haploid_snps] = False
            if np.all(~ix1):
                h0, h1, h2 = 0, 0 ,0
            else:
                h0, h1, h2 = np.mean(post_g[ix1], 0)

        f0, f1, f2 = g0 - h0, g1 - h1, g2 - h2
        #F[k] = 2 * f0 / (2 * f0 + f1 + 1e-100) - f1 / (2 * f2 + f1 + 1e-100)

        #np.clip(F, 0, 1, out=F)  # round to [0, 1]

        pars.tau[k] = g1 / 2.0 + g2
        ll, prev_ll = calc_full_ll(data, pars), ll
        print(f'{k} : ll = {ll}')
        if ll - prev_ll < 0:
            breakpoint()

@njit
def na_to_zero(x):
    y = x.reshape(np.prod(np.array((x.shape))))
    y[np.isnan(y)] = 0.0


@njit
def update_c(c, post_c, REF, ALT, OBS2SNP, OBS2RG, n_obs, n_rgs):
    """update c
        parameters:
        c - vector of contamination rates to updates
        bwd_x : Pr(O | X)    [2 x 2] array
        fwd_x_cont Σ_A Pr(X | C=1, A)Pr(A)   [L x 1] array
        fwd_x_nocont Σ_G Pr(X | C=1, G) Pr(G) [L x 1] array
        O, N observed alt and tot reads; [O x 1]
    """

    # np.nan_to_num(post_c, copy= False)
    q = np.vstack((REF, ALT)).T

    post_c *= np.expand_dims(q, 1)

    for r in range(n_rgs):
        x0 = np.sum(post_c[OBS2RG == r, 0])
        x1 = np.sum(post_c[OBS2RG == r, 1])
        c[r] = x1 / (x0 + x1)
        # print(r, c[r], np.sum(OBS2RG==r))

    return c

def calc_full_ll(data, pars):
    fwd_g = slug_p_gt_diploid(pars.tau, pars.F)[data.SNP2SFS]  # size [L x 3]
    if data.haploid_snps is not None:
        fwd_g[data.haploid_snps] = slug_p_gt_haploid(pars.tau)[data.SNP2SFS[data.haploid_snps]]  # size [L x 3]

    fwd_x_cont = slug_fwd_p_x_cont(data.psi)  # size [L x 1]
    fwd_x_nocont = slug_fwd_p_x_nocont(fwd_g)  # size [L x 1]
    fwd_x = slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, pars.cont, data)

    bwd_x = slug_bwd_p_o_given_x(pars.e, pars.b)  # size [2 x 2]
    return calc_ll(fwd_x, bwd_x, data.REF, data.ALT)

def calc_ll(fwd_x, bwd_x, REF, ALT):
    ll = np.sum(
        np.log((1 - fwd_x) * bwd_x[0, 0] + fwd_x * bwd_x[1, 0] + 1e-100) * REF
    )
    ll += np.sum(np.log((1 - fwd_x) * bwd_x[0, 1] + fwd_x * bwd_x[1, 1] + 1e-100) * ALT)
    return ll


def update_eb(post_x, ALT, REF, two_errors=True):
    post_x[:, 0] *= np.expand_dims(REF, 1)
    post_x[:, 1] *= np.expand_dims(ALT, 1)
    tmp = np.sum(post_x, 0)

    if two_errors:
        # ΣPr(X=1 | O = 0) / ΣPr(X=0|O=0) + ΣPr(X=1|O=0)
        e = tmp[1, 0] / np.sum(tmp[:, 0])
        b = tmp[0, 1] / np.sum(tmp[:, 1])
    else:
        e = (tmp[1, 0] + tmp[0, 1]) / np.sum(tmp)
        b = e

    if np.isnan(e):
        e = 0
    if np.isnan(b):
        b = 0
    return e, b


def em(pars, data, controller):
    for i in range(controller.n_iter):
        update_pars(pars, data, controller)
        s = f"iter {i}: Δll : {pars.delta_ll:.6f} | e={pars.e:.4f} | b={pars.b:.4f}"
        s += f" | Δc : {pars.delta_cont:.4f} | Δtau : {pars.delta_tau:.4f}"
        print(s)
        if pars.ll - pars.prev_ll < controller.ll_tol:
            break
    return pars


def update_pars(pars, data, controller):
    """estimate tau, F, e, b, c;  one step of em algorithm
    """

    """update tau, F
    need fwd_g, bwd_g,  IX.n_sfs, IX.SNP2SFS

    updates F, tau
    """
    O = controller
    IX = data

    fwd_g = slug_p_gt_diploid(pars.tau, pars.F)[IX.SNP2SFS]  # size [L x 3]
    if IX.haploid_snps is not None:
        fwd_g[IX.haploid_snps] = slug_p_gt_haploid(pars.tau)[IX.SNP2SFS[IX.haploid_snps]]  # size [L x 3]

    fwd_a = data.psi  # size [L x 1]
    fwd_c = pars.cont  # size [O x 1]

    bwd_x = slug_bwd_p_o_given_x(pars.e, pars.b)  # size [2 x 2]

    fwd_x_cont = slug_fwd_p_x_cont(fwd_a)  # size [L x 1]
    fwd_x_nocont = slug_fwd_p_x_nocont(fwd_g)  # size [L x 1]
    fwd_x = slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, IX)

    if O.do_ll:
        pars.ll, pars.prev_ll = calc_ll(fwd_x, bwd_x, data.REF, data.ALT), pars.ll

    if O.do_update_eb:
        post_x = slug_post_x(fwd_x, bwd_x)  # [O x 2 x 2]
        pars.prev_e, pars.prev_b = pars.e, pars.b
        pars.e, pars.b = update_eb(post_x, data.ALT, data.REF)
        del post_x
        if np.abs(pars.e + pars.prev_e) + np.abs(pars.b + pars.prev_b) < 0.0001:
            print("stopping error updates")
            O.do_update_eb = False

    if O.do_update_cont:
        post_c = slug_post_c(
            bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, IX.OBS2SNP, IX.OBS2RG, IX.n_obs
        )  # [O x 2 x 2]
        pars.prev_cont[:] = pars.cont
        update_c(
            pars.cont, post_c, data.REF, data.ALT, IX.OBS2SNP, IX.OBS2RG, IX.n_obs, IX.n_rgs
        )
        del post_c
        # if np.sum(np.abs(pars.cont - pars.prev_cont)) < 0.0001:
        #    print("stopping cont updates")
        #    O.do_update_cont = False

    if O.do_update_ftau:
        bwd_gac = slug_bwd_p_o_given_gac(bwd_x)  # size [3 x 2]
        bwd_g = slug_bwd_p_o_given_g_nocont(
            bwd_gac, fwd_a, fwd_c, IX.OBS2SNP, IX.OBS2RG, IX.n_obs
        )  # size [O x 3]
        breakpoint()
        post_g = slug_post_g(
            bwd_g, fwd_g, data.ALT, data.REF, IX.OBS2SNP, IX.n_snps
        )  # size [L x 3]
        del bwd_gac, bwd_g
        pars.prev_F[:], pars.prev_tau[:] = pars.F, pars.tau
        #update_ftau(pars.F, pars.tau, post_g, IX)
        update_ftau_debug(data, pars, post_g)
        del post_g

        if np.sum(np.abs(pars.cont - pars.prev_cont)) < 0.0001:
            print("stopping ftau updates")
            O.do_update_ftau = False
