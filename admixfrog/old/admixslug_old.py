"""hackish implementaion of admixslug,
using stuff from admixfrog"""


def update_sfs_genotypes(G, P, IX, F, tau):
    G[:, 0] = F * (1 - tau) + (1 - F) * (1 - tau) * (1 - tau)
    G[:, 1] = (1 - F) * 2 * (1 - tau) * tau
    G[:, 2] = F * tau + (1 - F) * tau * tau
    return 0


def update_snp(E, G, P, IX, cont, error, gt_mode):
    cflat = np.array([cont[lib] for lib in P.lib])
    eflat = np.array([error[lib] for lib in P.lib])

    E[:] = p_snps_given_gt(P, cflat, eflat, IX, gt_mode)
    E[:] = G[IX.SNP2BIN] * E
    return 0


def update_ftau(PG, IX, F, tau):
    for i in range(IX.n_sfs):
        g0, g1, g2 = np.mean(PG[IX.SNP2BIN == i], 0)
        F[i] = 2 * g0 / (2 * g0 + g1 + 1e-8) - g1 / (2 * g2 + g1 + 1e-8)
        if np.isnan(F[i]):
            breakpoint()
        F[:] = np.minimum(1, np.maximum(0, F))
        tau[i] = g1 / 2 + g2


def em_sfs(
    P,
    IX,
    pars,
    est_options,
    max_iter=1000,
    ll_tol=1e-1,
    gt_mode=False,
    do_hmm=True,
    **kwargs
):
    O = est_options
    cont, error, F, tau, sex = pars
    ll = -np.inf
    n_gt = 3

    # create arrays for posterior, emissions
    G = np.zeros((IX.n_sfs, n_gt))  # Pr(G | F, tau)
    E = np.zeros((IX.n_snps, n_gt))  # P(O, G | tau, F)

    # update G: given F, tau, calculate Pr(G | F, tau) for each SFS class
    update_sfs_genotypes(G, P, IX, F, tau)

    # update E: Pr(O, G | F, tau)
    update_snp(E, G, P, IX, cont, error, gt_mode)

    for it in range(max_iter):
        ll, old_ll = np.sum(np.log(np.sum(E, 1))), ll
        if np.isnan(ll):
            breakpoint()
        assert not np.isnan(ll)

        tpl = (
            IX.n_reads / 1000,
            IX.n_obs / 1000,
            IX.n_snps / 1000,
            IX.n_bins / 1000,
            it,
            np.mean(tau),
            ll,
            ll - old_ll,
        )
        log_.info("[%dk|%dk|%dk|%dk]: iter:%d |tavg:%.3f\tLL:%.4f\tΔLL:%.4f" % tpl)
        if ll - old_ll < ll_tol:
            break

        if np.any(np.isnan(E)):
            raise ValueError("nan observed in emissions")

        E[:] = E / np.sum(E, 1)[:, np.newaxis]
        update_ftau(E, IX, F, tau)
        PG = E.reshape(E.shape[0], 1, E.shape[1])
        delta = update_contamination(cont, error, P, PG, IX, est_options)

        update_sfs_genotypes(G, P, IX, F, tau)
        update_snp(E, G, P, IX, cont, error, gt_mode)

    pars = SlugPars(cont, error, F, tau, sex)

    return G, PG, pars, ll
