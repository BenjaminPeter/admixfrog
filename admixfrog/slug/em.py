from numba import njit
import numpy as np
from copy import deepcopy
from .emissions import fwd_p_g
from .emissions import bwd_p_o_given_x
from .emissions import bwd_p_one_o_given_g, bwd_p_all_o_given_g
from .emissions import fwd_p_x, fwd_p_x_cont, message_fwd_p_x_nocont
from .emissions import posterior_g, posterior_x, posterior_c
from .emissions import full_posterior_genotypes, calc_ll, calc_full_ll
from ..utils.log import log_


def update_ftau(old_F, old_tau, data, post_g, update_F=True):
    """updates the SFS parameters

    tau represents the conditional SFS entry, i.e. what proportion of sites are
    derived in any particular SFS category

    F is the inbreeding coefficient, i.e. what proportion of sites coalesce
    before the SFS entry

    Parameters
    ----------
    old_F : array[S x 1]
        old F
    old_tau : array[S x 1]
        old tau
    data : SlugData
        data object, used for various indices
    post_g : array[L x 3]
        Pr(G | O)
    update_F : bool, optional
        should F be updated

    Returns
    -------
    new_F : array[S x 1]
    new_tau : array[Sx 1]
    """

    tau, F = np.zeros(data.n_sfs), np.zeros(data.n_sfs)
    tau[:], F[:] = old_tau, old_F

    for k in range(data.n_sfs):
        g0, g1, g2 = np.sum(post_g[(data.SNP2SFS == k) & (~data.FLIPPED)], 0)
        f2, f1, f0 = np.sum(post_g[(data.SNP2SFS == k) & (data.FLIPPED)], 0)

        G0, G1, G2 = g0 + f0, g1 + f1, g2 + f2

        # update tau
        if G0 + G1 + G2 > 0:
            tau[k] = (G1 / 2.0 + G2) / (G0 + G1 + G2)
        else:
            tau[k] = 0.0

        # for F, exclude X-chromosome stuff
        if update_F:
            if data.haploid_snps is None:
                pass
            else:
                ix1 = data.SNP2SFS == k
                ix1[data.haploid_snps] = False
                if np.all(~ix1):
                    pass
                else:
                    g0, g1, g2 = np.sum(post_g[(ix1) & (~data.FLIPPED)], 0)
                    f2, f1, f0 = np.sum(post_g[(ix1) & (data.FLIPPED)], 0)
                    G0, G1, G2 = g0 + f0, g1 + f1, g2 + f2

            F[k] = 2 * G0 / (2 * G0 + G1 + 1e-300) - G1 / (2 * G2 + G1 + 1e-300)
        np.clip(F, 0, 1, out=F)  # round to [0, 1]

    return F, tau


def update_c(post_c, READ2RG, n_rgs):
    """update c
    parameters:
    c - vector of contamination rates to updates
    post_c : Pr( C_i = k| O_i = j)
    REF, ALT: number of reference and alt alleles
    OBS2RG: indices for readgroup
    n_rgs: number of rgs
    """

    c = np.empty(n_rgs)

    for r in range(n_rgs):
        c[r] = np.mean(post_c[READ2RG == r])

    return c


def update_eb(post_x, R, two_errors=False):
    not_bias = np.sum(post_x[:, 1] * (R == 1))
    not_errors = np.sum(post_x[:, 0] * (R == 0))
    bias = np.sum(post_x[:, 1] * (R == 0))
    errors = np.sum(post_x[:, 0] * (R == 1))

    e = errors / (errors + not_errors)
    b = bias / (bias + not_bias)

    e = 0 if np.isnan(e) else e
    b = 0 if np.isnan(b) else b

    return e, b


def update_pars(pars, data, controller, latents=None):
    """update all parameters; 1 EM step"""
    O = controller

    pars = deepcopy(pars) if controller.copy_pars else pars

    """ calc unconditional forward probs Pr(G), Pr(C), Pr(A)"""
    fwd_g = fwd_p_g(data, pars)
    fwd_a = data.psi  # size [L x 1]
    fwd_c = pars.cont  # size [O x 1]

    """run backward algorithm to calculate Pr(O | .)"""
    bwd_x = bwd_p_o_given_x(data.READS, pars.e, pars.b)
    bwd_g1 = bwd_p_one_o_given_g(
        bwd_x, fwd_a, fwd_c, data.READ2SNP, data.READ2RG, data.n_reads
    )  # size [O x 3]
    bwd_g = bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    """remaining forward probs Pr(X| C=0, G), Pr(X | C=1, A), Pr(X | C, A G) """
    fwd_x_cont = fwd_p_x_cont(fwd_a, data.READ2SNP)  # size [L x 1]
    fwd_x_nocont = message_fwd_p_x_nocont(
        fwd_g, bwd_g, bwd_g1, data.READ2SNP
    )  # size [L x 1]

    post_c = posterior_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, data.READ2RG)
    post_x = posterior_x(bwd_x, fwd_x_cont, fwd_x_nocont, fwd_c, data.READ2RG)

    fwd_x = fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, data.READ2RG)
    post_g = posterior_g(bwd_g, fwd_g)

    if O.update_ftau:
        post_g = posterior_g(bwd_g, fwd_g)
        pars.prev_F[:], pars.prev_tau[:] = pars.F, pars.tau
        pars.F, pars.tau = update_ftau(
            pars.F, pars.tau, data, post_g, update_F=controller.update_F
        )

    if O.update_cont:
        post_c = posterior_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, data.READ2RG)
        pars.prev_cont[:] = pars.cont
        pars.cont = update_c(post_c, data.READ2RG, data.n_rgs)

    if O.update_eb:
        post_x = posterior_x(bwd_x, fwd_x_cont, fwd_x_nocont, fwd_c, data.READ2RG)
        pars.prev_e, pars.prev_b = pars.e, pars.b
        pars.e, pars.b = update_eb(post_x, data.READS, two_errors=O.update_bias)

    if O.do_ll:
        pars.ll, pars.prev_ll = calc_full_ll(data, pars), pars.ll

    return pars


def em(pars, data, controller):
    for i in range(controller.n_iter):
        update_pars(pars, data, controller)
        s = f"iter {i}: ll: {pars.ll:4f} | Δll : {pars.delta_ll:4f} | e={pars.e[0]:.4f} | b={pars.b[0]:.4f}"
        s += f" | Δc : {pars.delta_cont:.4f} | Δtau : {pars.delta_tau:.4f} | ΔF : {pars.delta_F:.4f}"
        log_._info(s)
        if pars.delta_ll < controller.ll_tol:
            break

    posterior_gt = full_posterior_genotypes(data, pars)
    return posterior_gt
