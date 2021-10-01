from numba import njit
import numpy as np
from .slug_emissions import *

def update_ftau(old_F, old_tau, data, post_g, update_F = True):
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

    tau, F = np.empty(data.n_sfs), np.empty(data.n_sfs)
    tau[:], F[:] = old_tau, old_F

    for k in range(data.n_sfs):
        g0, g1, g2 = np.mean(post_g[data.SNP2SFS == k], 0)

        tau[k] = (g1 / 2.0 + g2) / (g0 + g1 + g2)

        if update_F:
            if data.haploid_snps is None:
                f0, f1, f2 = 0, 0 ,0
            else:
                ix1 = data.SNP2SFS == k
                ix1[data.haploid_snps] = False
                if np.all(~ix1):
                    f0, f1, f2 = 0, 0 ,0
                else:
                    f0, f1, f2 = np.mean(post_g[ix1], 0)

            F[k] = 2 * f0 / (2 * f0 + f1 + 1e-300) - f1 / (2 * f2 + f1 + 1e-300)
            np.clip(F, 0, 1, out=F)  # round to [0, 1]
    
    return F, tau

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
    np.set_printoptions(suppress=True, precision=4)
    ll = - np.inf

    for k in range(data.n_sfs):
        g0, g1, g2 = np.mean(post_g[data.SNP2SFS == k], 0)
        if data.haploid_snps is None:
            h0, h1, h2 = 0, 0 ,0
        else:
            ix1 = data.SNP2SFS == k
            ix1[data.haploid_snps] = False
            if np.all(~ix1):
                f0, f1, f2 = 0, 0 ,0
            else:
                f0, f1, f2 = np.mean(post_g[ix1], 0)

        #f0, f1, f2 = g0 - h0, g1 - h1, g2 - h2
        pars.F[k], prev_F = 2 * f0 / (2 * f0 + f1 + 1e-300) - f1 / (2 * f2 + f1 + 1e-300), pars.F[k]

        np.clip(pars.F, 0, 1, out=pars.F)  # round to [0, 1]

        pars.tau[k], prev_tau = (g1 / 2.0 + g2 ) / (g0 + g1 + g2), pars.tau[k]
        ll, prev_ll = calc_full_ll(data, pars), ll
        print(f'{k} : ll = {ll}: Delta: {ll-prev_ll:.4f} | tau : {prev_tau:.4f} -> {pars.tau[k]:.4f}')
        if ll - prev_ll < 0:
            breakpoint()

#@njit
def update_c(post_c, REF, ALT, OBS2RG, n_rgs):
    """update c
        parameters:
        c - vector of contamination rates to updates
        post_c : Pr( C_i = k| O_i = j) 
        REF, ALT: number of reference and alt alleles
        OBS2RG: indices for readgroup
        n_rgs: number of rgs
    """

    c = np.empty(n_rgs)

    reads = np.vstack((REF, ALT)).T

    v = post_c * np.expand_dims(reads, 2)

    breakpoint()

    for r in range(n_rgs):
        #x = np.sum(post_c[OBS2RG == r], (0, 1))
        #x = np.sum(np.sum(post_c[OBS2RG == r], 0), 0)
        x = np.sum(np.sum(v[OBS2RG == r], 0), 0)
        c[r] = x[1] / (x[0] + x[1])
        #print(r, c[r], np.sum(OBS2RG==r))

    return c

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

def squarem(pars0, data, controller):
    EPS = 1e-6
    MIN_STEP_0, MAX_STEP_0 = 1., 1.
    MSTEP = 4.
    """squarem port from R"""

    controller.copy_pars = True #required for squarem
    controller.n_iter = 50
    min_step, max_step = MIN_STEP_0, MAX_STEP_0
    pars = pars0
    controller.do_update_ftau = False
    controller.do_update_cont = True
    controller.do_update_eb = True
    pars.e, pars.b = 0, 0
    for i in range(controller.n_iter):
        pars1 = update_pars(pars, data, controller)
        Δp1 = pars1 - pars
        if pars1.ll < pars.ll:
            breakpoint()
        if norm(Δp1) < EPS : # or pars1.ll < pars1.prev_ll: 
            pars = pars1
            break

        pars2 = update_pars(pars1, data, controller)
        Δp2 = pars2 - pars1
        if norm(Δp2) < EPS : # or pars2.ll < pars2.prev_ll: 
            pars = pars2
            break

        Δp3 = Δp2 - Δp1

        step_size = norm(Δp1) / norm(Δp3)
        print(f'basic step: {step_size} = {norm(Δp1)} / {norm(Δp3)}')
        step_size = np.clip(step_size, min_step, max_step)

        pars_sq = deepcopy(pars2)
        pars_sq.pars[:] = pars.pars + 2 * step_size * Δp1 + step_size * step_size * Δp3

        #"stabilization step if all parameters are in bounds
        if np.all(0  <= pars_sq.pars) and np.all(pars_sq.pars <= 1):
            pars_sq = update_pars(pars_sq, data, controller)
            print(f"stab step {pars2.ll} -> {pars_sq.ll}")
        else:
            pars_sq.pars[pars_sq.pars < 0] = EPS
            pars_sq.pars[pars_sq.pars > 1] = 1 - EPS
            pars_sq = update_pars(pars_sq, data, controller)
            print(f"out of bounds {pars2.ll} -> {pars_sq.ll} | {np.sum(pars_sq.pars < 0)}")
            #pars_sq = pars2

        print(f"LLs p0 {pars.ll:.4f} | p1 {pars1.ll-pars.ll:.4f} | p2 {pars2.ll-pars1.ll:.4f} | psq {pars_sq.ll-pars2.ll:.4f}")

        # ll did not improve
        if pars_sq.ll <= pars2.ll:
            print(f"bad improvement'{min_step} - {max_step} | {step_size}")
            pars = pars2
            if step_size >= max_step: 
                max_step = np.maximum(MAX_STEP_0, max_step / MSTEP)
                step_size = 1.
        else:
            pars = pars_sq
            print(f"good step!!!'{min_step} - {max_step} | {step_size}")

        if step_size == max_step:
            max_step *= MSTEP

        

        s =  f"iter {i}: step: {step_size:.3f} | ll: {pars.ll:4f}"
        s += f"Δll : {pars.delta_ll:.6f} | e={pars.e[0]:.4f} | b={pars.b[0]:.4f}"
        s += f" | Δc : {pars.delta_cont:.4f} | Δtau : {pars.delta_tau:.4f} | ΔF : {pars.delta_F:.4f}"
        print(s)

    posterior_gt = full_posterior_genotypes(data, pars)
    return pars, posterior_gt

def norm(self):
    return np.sqrt(np.sum(np.power(self, 2)))

def em(pars, data, controller):
    for i in range(controller.n_iter):
        update_pars(pars, data, controller)
        s = f"iter {i}: ll: {pars.ll:4f} | Δll : {pars.delta_ll:4f} | e={pars.e[0]:.4f} | b={pars.b[0]:.4f}"
        s += f" | Δc : {pars.delta_cont:.4f} | Δtau : {pars.delta_tau:.4f} | ΔF : {pars.delta_F:.4f}"
        print(s)
        if pars.ll - pars.prev_ll < controller.ll_tol:
            break


    posterior_gt = full_posterior_genotypes(data, pars)
    return posterior_gt


def update_pars(pars, data, controller):
    """update all parameters; 1 EM step
    """
    O = controller
    IX = data

    pars = deepcopy(pars) if controller.copy_pars else pars

    fwd_g = slug_p_gt_diploid(pars.tau, pars.F)[IX.SNP2SFS]  # size [L x 3]
    if IX.haploid_snps is not None:
        fwd_g[IX.haploid_snps] = slug_p_gt_haploid(pars.tau)[IX.SNP2SFS[IX.haploid_snps]]  # size [L x 3]

    fwd_a = data.psi  # size [L x 1]

    fwd_c = np.empty_like(pars.cont) #copy just to make sure no update effects
    fwd_c[:] = pars.cont  # size [O x 1]

    bwd_x = slug_bwd_p_o_given_x(pars.e, pars.b)  # size [2 x 2]

    fwd_x_cont = slug_fwd_p_x_cont(fwd_a)  # size [L x 1]
    fwd_x_nocont = slug_fwd_p_x_nocont(fwd_g)  # size [L x 1]
    fwd_x = slug_fwd_p_x(fwd_x_cont, fwd_x_nocont, fwd_c, IX)
    

    if O.do_update_ftau:

        bwd_g1 = slug_bwd_p_one_o_given_g(
            bwd_x, fwd_a, fwd_c, IX.OBS2SNP, IX.OBS2RG, IX.n_obs
        )  # size [O x 3]

        bwd_g = slug_bwd_p_all_o_given_g(bwd_g1, data.REF, data.ALT, data.OBS2SNP, data.n_snps)
        post_g = slug_post_g(bwd_g, fwd_g)
        pars.prev_F[:], pars.prev_tau[:] = pars.F, pars.tau
        pars.F, pars.tau = update_ftau(pars.F, pars.tau, data, post_g)

    if O.do_update_eb:
        post_x = slug_post_x(fwd_x, bwd_x)  # [O x 2 x 2]
        pars.prev_e, pars.prev_b = pars.e, pars.b
        pars.e, pars.b = update_eb(post_x, data.ALT, data.REF)

    if O.do_update_cont:
        post_c = slug_post_c(
            bwd_x, fwd_x_nocont, fwd_x_cont, fwd_c, IX.OBS2SNP, IX.OBS2RG, IX.n_obs
        )  # [O x 2 x 2]
        pars.prev_cont[:] = pars.cont
        pars.cont = update_c(post_c, data.REF, data.ALT, data.OBS2RG, data.n_rgs)

    if O.do_ll:
        pars.ll, pars.prev_ll = calc_full_ll(data, pars), pars.ll
        if pars.ll < pars.prev_ll:
            pass
            #breakpoint()


    return pars

