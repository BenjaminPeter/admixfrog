from numba import jit, prange
import numpy as np
import pandas as pd
import sys
from collections import namedtuple, defaultdict
from scipy.stats import binom
from scipy.optimize import minimize
from .hmm import get_emissions_cy, split_freqs, update_contamination_cy

np.set_printoptions(suppress=True, precision=4)
pbinom = binom.pmf
Freqs = namedtuple("Freqs", ("O", "N", "P_cont", "P", "lib"))

# @jit(nopython=True)
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


@jit(nopython=True)
def calc_ll(alpha0, trans_mat, emissions):
    """likelihood using forward algorithm"""
    _, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
    return np.sum([np.sum(np.log(n_)) for n_ in n])


@jit(nopython=True)
def fwd_step(alpha_prev, E, trans_mat):
    alpha_new = (alpha_prev @ trans_mat) * E
    n = np.sum(alpha_new)
    return alpha_new / n, n


@jit(nopython=True)
def fwd_algorithm_single_obs(alpha0, emission, trans_mat):
    """
    calculate P(X_t | o_[1..t], a0)
    =P(X_t , o_[1..t], a0 | o_[1..t])
    """
    n_steps, n_states = emission.shape
    alpha = np.empty((n_steps + 1, n_states))
    n = np.empty((n_steps + 1))
    alpha[0] = alpha0
    n[0] = np.sum(alpha0)
    for i in range(n_steps):
        alpha[i + 1], n[i + 1] = fwd_step(alpha[i], emission[i], trans_mat)
    return alpha, n


# @jit(nopython=True)
def fwd_algorithm(alpha0, emissions, trans_mat):
    """
    calculate P(X_t | o_[1..t], a0)
    =P(X_t , o_[1..t], a0 | o_[1..t])
    """
    # alpha, n = [], []
    n_seqs = len(emissions)
    alpha = [np.empty((2, 2)) for _ in range(n_seqs)]
    n = [np.empty((2)) for _ in range(n_seqs)]
    for i in range(n_seqs):
        alpha[i], n[i] = fwd_algorithm_single_obs(alpha0, emissions[i], trans_mat)
    return alpha, n


@jit(nopython=True)
def bwd_step(beta_next, E, trans_mat, n):
    beta = (trans_mat * E) @ beta_next
    return beta / n


@jit(nopython=True)
def bwd_algorithm_single_obs(emission, trans_mat, n):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n])
    """

    n_steps, n_states = emission.shape

    beta = np.empty((n_steps + 1, n_states))
    beta[n_steps] = 1

    # for i, e in zip(range(n_steps-1, -1, -1), reversed(em)):
    for i in range(n_steps - 1, -1, -1):
        beta[i] = bwd_step(beta[i + 1], emission[i], trans_mat, n[i + 1])
    return beta


# @jit(nopython=True)
def bwd_algorithm(emissions, trans_mat, n):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n])

    emissions : list[np.array] 
        list of emission probabilities, one per observed sequence (i.e. chromosome)
    n : list[np.array]
        list of normalization constants
    trans_mat : np.array<n_states x n_states>
        transition matrix
    """

    n_seqs = len(emissions)
    beta = [np.empty((2, 2)) for _ in range(n_seqs)]
    for i in range(n_seqs):
        n_i, em = n[i], emissions[i]
        n_steps, n_states = em.shape
        beta_i = bwd_algorithm_single_obs(em, trans_mat, n_i)
        beta[i] = beta_i
    return beta


# @jit(nopython=True)
def fwd_bwd_algorithm(alpha0, emissions, trans_mat):
    alpha, n = fwd_algorithm(alpha0=alpha0, emissions=emissions, trans_mat=trans_mat)
    beta = bwd_algorithm(emissions=emissions, n=n, trans_mat=trans_mat)
    gamma = [a * b for (a, b) in zip(alpha, beta)]
    # for g in gamma:
    #    assert np.allclose(np.sum(g, 1), 1)
    return alpha, beta, gamma, n


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


def data2freqs(data, state_ids, cont_id, do_hets=True):
    P_homo = [data[s] for s in state_ids]
    if do_hets:
        P_het = []
        for i, s in enumerate(state_ids):
            for s2 in state_ids[i + 1 :]:
                P_het.append((data[s] + data[s2]) / 2)
        P = np.vstack((P_homo, np.vstack(P_het))).T
    else:
        P = P_homo

    f = Freqs(
        O=np.array(data.alt),
        N=np.array(data.ref + data.alt),
        P_cont=np.array(data[cont_id]),
        lib=np.array(data.lib),
        P=P,
    )

    return f


def get_emissions(*args, **kwargs):
    return get_emissions_cy(*args, **kwargs)


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


def bins_from_bed(bed, data, bin_size):
    chroms = np.unique(bed.chrom)
    bin_loc, data_loc = [], []
    bin0 = 0

    for i, chrom in enumerate(chroms):
        pos = bed.pos[bed.chrom == chrom]
        pos_data = data.pos[data.chrom == chrom]

        min_ = int(np.floor(pos.head(1) / bin_size) * bin_size)
        max_ = int(np.ceil(pos.tail(1) / bin_size) * bin_size)

        bins = np.arange(min_, max_, bin_size, dtype="int")
        # dig = np.digitize(pos, bins, right=False) - 1
        dig_data = np.digitize(pos_data, bins, right=False) - 1

        bin_ids = range(bin0, bin0 + len(bins))

        bin_loc.append(np.vstack((np.zeros_like(bins) + i, bins, bin_ids)).T)
        # snp_id.append(np.vstack((np.zeros_like(dig)+i, dig)).T)
        data_loc.append(
            np.vstack((np.zeros_like(dig_data) + i, dig_data + bin0, dig_data)).T
        )

        bin0 += len(bins)

    bins = np.vstack(bin_loc)
    # snps = np.hstack((bed, np.vstack(snp_id)))
    data_bin = np.vstack(data_loc)
    return bins, data_bin


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


def viterbi(alpha0, trans_mat, emissions):
    return [viterbi_single_obs(alpha0, trans_mat, e) for e in emissions]


def viterbi_single_obs(alpha0, trans_mat, emissions):
    n_steps, n_states = emissions.shape

    ll = np.ones_like(emissions)
    backtrack = np.zeros_like(emissions, int)

    log_e = np.log(emissions)
    log_t = np.log(trans_mat)

    for i in range(n_steps):
        if i == 0:
            aux_mat = np.log(alpha0) + log_t + log_e[0]
        else:
            aux_mat = ll[i - 1] + log_t + log_e[i]
        ll[i] = np.max(aux_mat, 1)
        backtrack[i] = np.argmax(aux_mat, 1)

    path = np.empty(n_steps)
    cursor = np.argmax(ll[-1])
    for i in range(n_steps - 1, -1, -1):
        cursor = backtrack[i, cursor]
        path[i] = cursor

    return path


def run_hmm(
    infile,
    bedfile,
    split_lib=True,
    state_ids=("AFR", "NEA", "DEN"),
    cont_id="AFR",
    bin_size=1e4,
    garbage_state=True,
    do_hets=True,
    **kwargs
):
    """ run baum welch to find introgressed tracts
    infile formatting:
        - chrom, map: genomic coordinates
        - tref, talt: number of ref, alt reads
        - f_nea, f_afr: frequency of african and neandertal
    """
    data = pd.read_csv(infile)
    try:
        data.chrom[data.chrom == "X"] = 23
        data.chrom[data.chrom == "Y"] = 24
        data.chrom = data.chrom.astype(int)
    except TypeError:
        pass
    data = data[data.ref + data.alt > 0]
    states = list(set(list(state_ids) + [cont_id]))
    cols = ["chrom", "pos", "map", "lib", "ref", "alt"] + states
    data = data[cols]
    data = data.dropna()
    q = np.quantile(data.ref + data.alt, .99)
    data = data[data.ref + data.alt <= q]

    if "lib" not in data or (not split_lib):
        data = data.groupby(
            ("chrom", "pos", "NEA", "DEN", "AFR", "PAN"), as_index=False
        ).agg({"ref": sum, "alt": sum})
        q = np.quantile(data.ref + data.alt, .999)
        data = data[data.ref + data.alt <= q]
        data["lib"] = "lib0"

    # load bed
    bed = pd.read_table(bedfile, header=None)[[0, 2]]
    bed.columns = ["chrom", "pos"]
    try:
        bed.chrom[bed.chrom == "X"] = "23"
        bed.chrom[bed.chrom == "Y"] = "24"
        bed.chrom = bed.chrom.astype(int)
    except TypeError:
        pass

    n_states = len(state_ids)
    n_states = int(n_states + n_states * (n_states - 1) / 2) + garbage_state
    alpha_0 = np.array([1 / n_states] * n_states)
    trans_mat = np.zeros((n_states, n_states)) + 2e-2
    np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
    cont = defaultdict(lambda: 1e-2)

    bins, bin_data = bins_from_bed(bed, data, bin_size)
    freqs = data2freqs(data, state_ids, cont_id, do_hets=do_hets)

    # get state names
    n_homo = [s for s in state_ids]
    if do_hets:
        n_het = []
        for i, s in enumerate(state_ids):
            for s2 in state_ids[i + 1 :]:
                n_het.append(s + s2)
        gamma_names = n_homo + n_het
    else:
        gamma_names = n_homo

    if garbage_state:
        gamma_names += ["GARBAGE"]

    gamma, ll, trans_mat, alpha_0, cont, emissions = baum_welch(
        alpha_0=alpha_0,
        trans_mat=trans_mat,
        bins=bins,
        bin_data=bin_data,
        freqs=freqs,
        garbage_state=garbage_state,
        cont=cont,
        gamma_names=gamma_names,
        **kwargs
    )

    viterbi_path = viterbi(alpha_0, trans_mat, emissions)
    viterbi_df = pd.Series(np.hstack(viterbi_path), name="viterbi")

    gamma_out = np.hstack((bins, np.vstack([g[1:] for g in gamma])))
    gamma_out = pd.DataFrame(
        gamma_out, columns=["chrom_id", "bin_pos", "bin_id"] + gamma_names
    )
    gamma_out = pd.concat(
        (gamma_out.reset_index(drop=True), viterbi_df.reset_index(drop=True)), axis=1
    )
    data_out = pd.concat(
        (
            data.reset_index(drop=True),
            gamma_out.iloc[bin_data[:, 1]].reset_index(drop=True),
        ),
        axis=1,
    )

    return gamma_out, data_out, emissions, (ll, trans_mat, cont, alpha_0)
