import numpy as np
def baum_welch(
    alpha_0,
    trans_mat,
    bins,
    bin_data,
    freqs,
    # bed, data, bin_size,
    cont,
    error=1e-2,
    max_iter=2000,
    ll_tol=1e-1,
    bad_snp_cutoff=1e-10,
    optimize_cont=True,
    garbage_state=True,
    gamma_names=None,
):

    print("BNS %s" % bad_snp_cutoff, file=sys.stderr)
    n_states = trans_mat.shape[0]
    ll = -np.inf

    assert len(alpha_0) == n_states
    assert trans_mat.shape[1] == n_states
    assert np.allclose(np.sum(trans_mat, 1), 1)

    emissions = get_emissions(
        cont,
        bins,
        bin_data,
        freqs,
        e=error,
        bad_snp_cutoff=bad_snp_cutoff,
        garbage_state=garbage_state,
    )
    n_seqs = len(emissions)

    if optimize_cont:
        libs = pd.unique(freqs.lib)
        split_ids = split_freqs(libs, freqs)
    for it in range(max_iter):
        alpha, beta, gamma, n = fwd_bwd_algorithm(alpha_0, emissions, trans_mat)
        ll, old_ll = np.sum([np.sum(np.log(n_i)) for n_i in n]), ll

        print("iter %d [%d/%d]: %s -> %s" % (it, n_seqs, n_states, ll, ll - old_ll))
        if ll - old_ll < ll_tol:
            break

        # update stuff
        alpha_0 = np.linalg.matrix_power(trans_mat, 10000)[0]
        trans_mat = update_transitions(trans_mat, alpha, beta, gamma, emissions, n)
        if optimize_cont:
            cont = update_contamination_cy(
                cont,
                error,
                bin_data,
                freqs,
                gamma,
                split_ids,
                libs,
                garbage_state=garbage_state,
            )
            # cont = update_contamination_py(cont, error, bin_data, freqs, gamma,bad_snp_cutoff, garbage_state=True)
            emissions = get_emissions(
                cont,
                bins,
                bin_data,
                freqs,
                bad_snp_cutoff=bad_snp_cutoff,
                e=error,
                garbage_state=garbage_state,
            )

        if gamma_names is not None:
            print(*gamma_names, sep="\t")
        print(*["%.3f" % a for a in alpha_0], sep="\t")

    return gamma, ll, trans_mat, alpha_0, dict(cont), emissions

def get_emissions_py(
    cont, bins, bin_data, freqs, e=1e-2, bad_snp_cutoff=1e-10, garbage_state=True
):
    n_snps = len(freqs.O)
    n_steps = bins.shape[0]
    n_states = freqs.P.shape[1] + garbage_state

    emissions = np.ones((n_steps, n_states))
    c = np.array([cont[l] for l in freqs.lib])

    for s in range(n_states):

        if garbage_state and s == n_states - 1:
            break

        p = freqs.P_cont * c + freqs.P[:, s] * (1. - c)
        p = p * (1. - e) + (1. - p) * e
        em = pbinom(freqs.O, freqs.N, p)

        for i in range(n_snps):
            row = bin_data[i, 1]
            emissions[row, s] *= em[i]

    if garbage_state:
        emissions[:, n_states - 1] = bad_snp_cutoff
    else:
        emissions = np.maximum(emissions, bad_snp_cutoff)

    chroms = np.unique(bins[:, 0])
    emissions = [emissions[bins[:, 0] == chrom] for chrom in chroms]
    return emissions



def update_transitions(old_trans_mat, alpha, beta, gamma, emissions, n):
    new_trans_mat = np.zeros_like(old_trans_mat)
    n_states = old_trans_mat.shape[0]
    # update transition
    for i in range(n_states):
        for j in range(n_states):
            for a, b, e, n_ in zip(alpha, beta, emissions, n):
                new_trans_mat[i, j] += np.sum(
                    a[:-1, i] * old_trans_mat[i, j] * b[1:, j] * e[:, j] / n_[1:]
                )

    gamma_sum = np.sum([np.sum(g[:-1], 0) for g in gamma], 0)
    new_trans_mat /= gamma_sum[:, np.newaxis]

    # underflow due to absence of state
    if not np.allclose(np.sum(new_trans_mat, 1), 1):
        for i in range(n_states):
            if np.any(np.isnan(new_trans_mat[i])):
                new_trans_mat[i] = 0.
                new_trans_mat[i, i] - 1.

    return new_trans_mat
