
def get_po_given_zc(c, e, freqs, bad_snp_cutoff=1e-10):
    """
    calculate Pr(O | Z, c, p), the probability of all observations given 
    the hidden states, allele frequencies in all potential donors and contamination rate

    - this is currenlty used to maximize the likelihood of c
    - the binomial coefficient is omitted
    - error is ignored

    cont : contamination estimate, by library
    freqs : allele frequencies, reads, by library/chromosome
    """
    n_states = freqs.P.shape[1]
    p = c * freqs.P_cont[:, np.newaxis] + (1. - c) * freqs.P
    p = p * (1 - e) + (1 - p) * e
    return (
        p ** freqs.O[:, np.newaxis] * (1. - p) ** (freqs.N - freqs.O)[:, np.newaxis],
    )
    return np.maximum(
        p ** freqs.O[:, np.newaxis] * (1. - p) ** (freqs.N - freqs.O)[:, np.newaxis],
        bad_snp_cutoff / n_states,
    )

def update_contamination_py(
    cont, error, bin_data, freqs, gamma, bad_snp_cutoff=1e-10, garbage_state=True
):
    """
    update emissions by maximizing contamination parameter

    cont: list of contamination rates (by library)
    gamma : Pr(Z | O, c_prev)
    snp_to_z : data structure connection snp to hidden state/ chrom coordinate
    data:

    we split up SNP 1. by library, 2. by sequence
    data is organized as data[lib_id][seq_id][state_id],
    i.e. the entry of data[lib_id][seq_id][state_id] is a Freq(O, N, p_cont, p_state)
        object


    """
    libs = pd.unique(freqs.lib)
    for i, lib in enumerate(libs):
        f_ = freqs.lib == lib
        O, N, P_cont, P = freqs.O[f_], freqs.N[f_], freqs.P_cont[f_], freqs.P[f_]

        G = np.array([gamma[i][j + 1] for i, _, j in bin_data[f_]])
        # print(lib, np.sum(G, 0), end = "\t")

        if garbage_state:
            G = G[:, :-1]

        # numeric minimization
        def get_po_given_zc_all(c):
            f = Freqs(O, N, P_cont, P, lib)
            po_given_zc = get_po_given_zc(c, error, f, bad_snp_cutoff)
            prob = np.sum(np.log(np.sum(G * po_given_zc, 1)))
            # print("[%s]minimizing c:" % lib, c, prob)
            return -prob

        p0 = get_po_given_zc_all(cont[lib])
        OO = minimize(get_po_given_zc_all, [cont[lib]], bounds=[(0., 1)])
        print(
            "[%s/%s]minimizing c: [%.3f->%.3f]: %4f, %4f"
            % (lib, len(O), cont[lib], OO.x[0], p0, OO.fun)
        )
        cont[lib] = OO.x[0]
    return dict(cont)
