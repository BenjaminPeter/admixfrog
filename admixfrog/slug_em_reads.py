from numba import njit
import numpy as np
from copy import deepcopy
from .slug_emissions_reads import slug_fwd_p_g
from .slug_emissions_reads import slug_bwd_p_o_given_x
from .slug_emissions_reads import slug_bwd_p_one_o_given_g,slug_bwd_p_all_o_given_g
from .slug_emissions_reads import slug_fwd_p_x, slug_fwd_p_x_cont, message_fwd_p_x_nocont
from .slug_emissions_reads import slug_post_g, slug_post_x, slug_post_c
from .slug_emissions_reads import full_posterior_genotypes, calc_ll, calc_full_ll_reads


def norm(self):
    """simple L2 - norm function to test parameter vector convergence and squarem steps"""
    return np.sqrt(np.sum(np.power(self, 2)))

def update_ftau(old_F, old_tau, data, post_g, update_F = True):
    """updates the SFS parameters

    tau represents the conditional SFS entry, i.e. what proportion of sites are derived
    in any particular SFS category

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
        g0, g1, g2 = np.mean(post_g[data.SNP2SFS == k], 0)
        tau[k] = (g1 / 2.0 + g2) / (g0 + g1 + g2)
    
    return F, tau

@njit
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
    not_bias = np.sum(post_x[:,1] * (R == 1))
    not_errors = np.sum(post_x[:,0] * (R == 0))
    bias = np.sum(post_x[:,1] * (R == 0))
    errors = np.sum(post_x[:,0] * (R == 1))

    e = errors / (errors + not_errors)
    b = bias / (bias + not_bias)

    e = 0 if np.isnan(e) else e
    b = 0 if np.isnan(b) else b


    return e, b

def update_pars_reads(pars, data, controller):
    """update all parameters; 1 EM step
    """
    O = controller

    pars = deepcopy(pars) if controller.copy_pars else pars

    """ calc unconditional forward probs Pr(G), Pr(C), Pr(A)"""
    fwd_g = slug_fwd_p_g(data, pars)
    fwd_a = data.psi  # size [L x 1]
    fwd_c = pars.cont  # size [O x 1]

    """run backward algorithm to calculate Pr(O | .)"""
    bwd_x = slug_bwd_p_o_given_x(data.READS, pars.e, pars.b)  
    bwd_g1 = slug_bwd_p_one_o_given_g(
        bwd_x, fwd_a, fwd_c, 
        data.READ2SNP, data.READ2RG, data.n_reads
    )  # size [O x 3]
    bwd_g = slug_bwd_p_all_o_given_g(bwd_g1, data.READ2SNP, data.n_snps)

    """remaining forward probs Pr(X| C=0, G), Pr(X | C=1, A), Pr(X | C, A G) """
    fwd_x_cont = slug_fwd_p_x_cont(fwd_a, data.READ2SNP)  # size [L x 1]
    fwd_x_nocont = message_fwd_p_x_nocont(fwd_g, bwd_g, bwd_g1, data.READ2SNP)  # size [L x 1]
    
    post_c = slug_post_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, data.READ2RG)
    post_x = slug_post_x(bwd_x, fwd_x_cont, fwd_x_nocont, fwd_c, data.READ2RG) 
    from .slug_emissions_reads import slug_fwd_p_x
    fwd_x = slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, data.READ2RG)
    post_g = slug_post_g(bwd_g, fwd_g)
    #breakpoint()

    if O.do_update_ftau:
        post_g = slug_post_g(bwd_g, fwd_g)
        pars.prev_F[:], pars.prev_tau[:] = pars.F, pars.tau
        pars.F, pars.tau = update_ftau(pars.F, pars.tau, data, post_g)

    if O.do_update_cont:
        post_c = slug_post_c(bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, data.READ2RG)
        pars.prev_cont[:] = pars.cont
        pars.cont = update_c(post_c, data.READ2RG, data.n_rgs)

    if O.do_update_eb:
        post_x = slug_post_x(bwd_x, fwd_x_cont, fwd_x_nocont, fwd_c, data.READ2RG) 
        pars.prev_e, pars.prev_b = pars.e, pars.b
        pars.e, pars.b = update_eb(post_x, data.READS, two_errors = True)


    if O.do_ll:
        pars.ll, pars.prev_ll = calc_full_ll_reads(data, pars), pars.ll
        if pars.ll < pars.prev_ll:
            pass


    return pars

def squarem(pars0, data, controller):
    EPS = 1e-5
    MIN_STEP_0, MAX_STEP_0 = 1., 1.
    MSTEP = 4.
    """squarem port from R"""

    controller.copy_pars = True #required for squarem
    controller.n_iter = 50
    min_step, max_step = MIN_STEP_0, MAX_STEP_0
    pars = pars0
    controller.do_update_ftau = True
    controller.do_update_cont = True
    controller.do_update_eb = True

    for i in range(controller.n_iter):
        pars1 = update_pars_reads(pars, data, controller)
        Δp1 = pars1 - pars
        if norm(Δp1) < EPS: #  or pars1.ll - pars1.prev_ll < EPS: 
            if pars1.ll < pars1.prev_ll:
                breakpoint()
            pars = pars1
            break

        pars2 = update_pars_reads(pars1, data, controller)
        Δp2 = pars2 - pars1
        if norm(Δp2) < EPS: #  or pars2.ll - pars2.prev_ll < EPS: 
            if pars2.ll < pars2.prev_ll:
                breakpoint()
            pars = pars2
            break

        Δp3 = Δp2 - Δp1

        step_size = norm(Δp1) / norm(Δp3)
        step_size = np.clip(step_size, min_step, max_step)

        pars_sq = deepcopy(pars2)
        pars_sq.pars[:] = pars.pars + 2 * step_size * Δp1 + step_size * step_size * Δp3

        #"stabilization step if all parameters are in bounds
        if np.all(0  <= pars_sq.pars) and np.all(pars_sq.pars <= 1):
            pass
        else:
            pars_sq.pars[pars_sq.pars < 0] = EPS
            pars_sq.pars[pars_sq.pars > 1] = 1 - EPS

        pars_sq = update_pars_reads(pars_sq, data, controller)

        print(f"LLs p0 {pars.ll:.4f} | p1 {pars1.ll-pars.ll:.4f} | p2 {pars2.ll-pars1.ll:.4f} | psq {pars_sq.ll-pars2.ll:.4f}")

        # ll did not improve
        if pars_sq.ll <= pars2.ll:
            pars = pars2
            if step_size >= max_step: 
                max_step = np.maximum(MAX_STEP_0, max_step / MSTEP)
        else:
            pars = pars_sq
            if step_size == max_step:
                max_step *= MSTEP

        

        s =  f"iter {i}: step: {step_size:.3f} | ll: {pars.ll:4f}"
        s += f"Δll : {pars.delta_ll:.6f} | e={pars.e[0]:.4f} | b={pars.b[0]:.4f}"
        s += f" | Δc : {pars.delta_cont:.4f} | Δtau : {pars.delta_tau:.4f} | ΔF : {pars.delta_F:.4f}"
        s += f" | Δ1 : {norm(Δp1)}| Δ2:{norm(Δp2)} "
        print(s)

    posterior_gt = full_posterior_genotypes(data, pars)
    return pars, posterior_gt


def em(pars, data, controller):
    for i in range(controller.n_iter):
        update_pars_reads(pars, data, controller)
        s = f"iter {i}: ll: {pars.ll:4f} | Δll : {pars.delta_ll:4f} | e={pars.e[0]:.4f} | b={pars.b[0]:.4f}"
        s += f" | Δc : {pars.delta_cont:.4f} | Δtau : {pars.delta_tau:.4f} | ΔF : {pars.delta_F:.4f}"
        print(s)
        if pars.ll - pars.prev_ll < controller.ll_tol:
            break


    posterior_gt = full_posterior_genotypes(data, pars)
    return posterior_gt


