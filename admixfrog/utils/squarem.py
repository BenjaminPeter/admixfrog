import numpy as np
from copy import deepcopy
import logging


def norm(self):
    """simple L2 - norm function to test parameter vector convergence and squarem steps"""
    return np.sqrt(np.nansum(np.power(self, 2)))


def squarem(pars0, data, options, updater, latents=None):
    """squarem port from R
    updater is a  function(pars, data, latents, options) that will calculate one EM iteration
    """

    EPS = options.param_tol
    min_step = options.squarem_min
    max_step = options.squarem_max
    pars = pars0

    for i in range(options.n_iter):
        pars1 = updater(pars, data, options, latents)
        Δp1 = pars1 - pars
        if norm(Δp1) < EPS:  #  or pars1.ll - pars1.prev_ll < EPS:
            pars = pars1
            break

        pars2 = updater(pars1, data, options, latents)
        Δp2 = pars2 - pars1
        if norm(Δp2) < EPS or pars2.ll - pars1.ll < EPS:
            pars = pars2
            break

        Δp3 = Δp2 - Δp1

        step_size = norm(Δp1) / norm(Δp3)
        step_size = np.clip(step_size, min_step, max_step)

        pars_sq = deepcopy(pars2)
        pars_sq.pars[:] = pars.pars + 2 * step_size * Δp1 + step_size * step_size * Δp3

        # "stabilization step if all parameters are in bounds
        if np.all(0 <= pars_sq.pars) and np.all(pars_sq.pars <= 1):
            pass
        else:
            pars_sq.pars[pars_sq.pars < 0] = EPS
            pars_sq.pars[pars_sq.pars > 1] = 1 - EPS

        pars_sq = updater(pars_sq, data, options, latents)

        logging.info(
            f"LLs p0 {pars.ll:.4f} | p1 {pars1.ll-pars.ll:.4f} | p2 {pars2.ll-pars1.ll:.4f} | psq {pars_sq.ll-pars2.ll:.4f}"
        )

        # ll did not improve
        if pars_sq.ll <= pars2.ll:
            pars = pars2
            if step_size >= max_step:
                max_step = np.maximum(
                    options.squarem_max, max_step / options.squarem_mstep
                )
        else:
            pars = pars_sq
            if step_size == max_step:
                max_step *= options.squarem_mstep

        s = f"iter {i}: step: {step_size:.3f} | ll: {pars.ll:4f} "
        s += f"Δll : {pars.delta_ll:.6f} | "
        s += f" | Δ1 : {norm(Δp1):.4f}| Δ2:{norm(Δp2):.4f} "
        logging.info(s)

    return pars
