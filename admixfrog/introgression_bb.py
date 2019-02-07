import numpy as np
import pandas as pd
import sys
from collections import namedtuple, defaultdict
from scipy.stats import binom
from scipy.optimize import minimize
from .hmm import get_emissions_cy, split_freqs, update_contamination_cy
from .fwd_bwd import viterbi, fwd_bwd_algorithm
from .baum_welch import update_transitions, get_emissions_py, baum_welch
from .utils import data2probs, bins_from_bed, init_pars, Probs
from .distributions import dbetabinom

np.set_printoptions(suppress=True, precision=4)


def get_emissions(*args, **kwargs):
    return get_emissions_cy(*args, **kwargs)

def baum_welch_bb(
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

def p_reads_given_gt(P, n_obs, c, error, has_hap = False):
    read_emissions = np.ones((n_obs, 5)) if has_hap else np.ones((n_obs, 3))
    for g in range(3):
        p = c * P.P_cont + (1 - c) * g / 2
        p = p * (1-error) + (1-p) * error
        read_emissions[:,g] = binom.pmf(P.O, P.N, p)

    if has_hap:
        p = c * P.P_cont + (1 - c) 
        p = p * (1-error) + (1-p) * error
        read_emissions[:,3] = binom.pmf(P.O, P.N, 1-p)
        read_emissions[:,4] = binom.pmf(P.O, P.N, p)
    return read_emissions


def get_emissions_bb_py(P,
                        bins,
                        bin_data,
                        cont,
                        prior_alpha, prior_beta, fst,
                        error=2e-2, bad_bin_cutoff=1e-100,
                        has_hap = False,
):
    n_homo_states = P.alpha.shape[1]
    n_het_states = int(n_homo_states * (n_homo_states - 1) / 2)
    n_hap_states = n_homo_states * has_hap
    n_states = n_homo_states + n_het_states + n_hap_states
    n_obs = P.O.shape[0]
    n_snps = P.alpha.shape[0]
    n_bins = bins.shape[0]
    c = np.array([cont[l] for l in P.lib])

    read_emissions = p_reads_given_gt(P, n_obs, c, error, has_hap)
    #P(SNP | GT)
    gt_emissions = np.ones((n_snps, 5)) if has_hap else np.ones((n_snps, 3))
    for i, row in enumerate(bin_data['snp_id']):
        gt_emissions[row] *= read_emissions[i]


    GT = np.empty((n_snps, n_states-n_hap_states, 3))
    #P(GT | Z)
    for s in range(n_homo_states):
        for g in range(3):
            GT[:,s, g] = dbetabinom(np.zeros(n_snps, int) + g,
                            np.zeros(n_snps, int) + 2,
                            fst[s] * (P.alpha[:, s] + prior_alpha),
                            fst[s] * (P.beta[:, s] + prior_beta))

    s = n_homo_states
    fa, fb = (P.alpha + prior_alpha).T, (P.beta + prior_beta).T
    pp = fa / (fa + fb)
    qq = 1 - pp
    for s1 in range(n_homo_states):
        for s2 in range(s1+1, n_homo_states):
            for g in range(3):
                GT[:, s, 0] = qq[s1] * qq[s2]
                GT[:, s, 2] = pp[s1] * pp[s2]
                GT[:, s, 1] = 1 - GT[:, s, 0] - GT[:, s, 2]
            s += 1
    assert np.allclose( np.sum(GT, 2), 1)

    #haploid regions
    GH = np.empty((n_snps, n_hap_states, 2))
    for s in range(n_hap_states):
        GH[:, s, 0] = qq[s]
        GH[:, s, 1] = pp[s]


    snp_emissions = np.sum(GT * gt_emissions[:, np.newaxis, :3],2)   
    if has_hap:
        snp_emissions_h = np.sum(GH * gt_emissions[:, np.newaxis, 3:],2)   
        snp_emissions = np.hstack((snp_emissions, snp_emissions_h))

    bin_emissions = np.ones((n_bins, n_states))
    for bin_, snp in np.unique(bin_data[['bin_id','snp_id']],axis=0):
        bin_emissions[bin_] *= snp_emissions[snp]

    bad_bins = np.sum(bin_emissions, 1) < bad_bin_cutoff
    bin_emissions[bad_bins] += bad_bin_cutoff / n_states

    chroms = np.unique(bins['chrom'])
    bin_emissions = [bin_emissions[bins['chrom'] == chrom] for chrom in chroms]
    return bin_emissions


def p_gt_given_zo():

    
