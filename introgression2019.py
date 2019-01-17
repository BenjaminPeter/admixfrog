#!python
#cython: language_level=3
# cython: infer_types=True
# distutils: language = c++

from scipy.stats import binom
from scipy.optimize import fminbound, minimize
from collections import defaultdict, Counter
from hmm import  er_given_zoc, init_global_po
import pandas as pd
import numpy as np
from pprint import pprint
from numba import jit
import random
import pdb

np.set_printoptions(suppress=True, precision=4)          
pbinom = binom.pmf

@jit(nopython=True)
def get_po_given_zc(cont, obs, N, p_cont, p_state, bad_snp_cutoff):
    p = cont * p_cont + (1.-cont) * p_state
    return np.maximum(p ** obs * (1.-p) ** (N - obs), bad_snp_cutoff)
    #return pbinom(obs, N, p)



binned = False
n_states = 6

@jit(nopython=True)
def calc_ll(alpha0, trans_mat, emissions):
    _, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
    return np.sum([np.sum(np.log(n_)) for n_ in n])


@jit(nopython=True)
def fwd_step(alpha_prev, E, trans_mat):
    alpha_new = (alpha_prev @ trans_mat) * E
    n = np.sum(alpha_new)
    return alpha_new / n, n

@jit(nopython=True)
def fwd_algorithm(alpha0, emissions, trans_mat):
    """
    calculate P(X_t | o_[1..t], a0)
    =P(X_t , o_[1..t], a0 | o_[1..t])
    """

    alpha, n = [], []

    for em in emissions:
        print(em, em.shape)
        n_steps, n_states = em.shape
        alpha_i = np.empty((n_steps+1, n_states))
        n_i = np.empty((n_steps+1))

        alpha_i[0] = alpha0
        n_i[0] = np.sum(alpha0)

        #for i, e in enumerate(em):
        for i in range(em.shape[0]):
            e = em[i]
            alpha_i[i + 1], n_i[i + 1] = fwd_step(alpha_i[i], e, trans_mat)

        alpha.append(alpha_i)
        n.append(n_i)

    return alpha, n

@jit(nopython=True)
def bwd_step(beta_next, E, trans_mat, n):
    beta = (trans_mat * E) @ beta_next
    return beta / n

@jit(nopython=True)
def bwd_algorithm(emissions, trans_mat, n):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n])
    """

    beta = []
    for chrom, em in enumerate(emissions):
        n_steps, n_states = em.shape

        beta_i = np.empty((n_steps + 1, n_states))
        beta_i[n_steps] = 1

        #for i, e in zip(range(n_steps-1, -1, -1), reversed(em)):
        for i in range(n_steps-1, -1, -1):
            e = em[i]
            beta_i[i] = bwd_step(beta_i[i+1], e, trans_mat, n[chrom][i+1])
        beta.append(beta_i)
    return beta

def split_by_lib(var, data, libs):
    return [var[data == l] for l in libs]

def baum_welch(data, E0, C0, bin_size=1e4, bad_snp_cutoff=1e-14, max_iter=2000,
               estimate_c = True, estimate_c_em = True,
               split_lib=True, citer_likelihood=False):
    """
    baum welch algorithm using multiomial probabilities.
    at a site with allele frequency Nea:p1, Deni:p2, Afr:p3
    Delta1 = p1 - p3
    Delta2 = p2 - p3

    P(X|NN) = p1 - c Delta1
    P(X|DD) = p2 - c Delta2
    P(X|AA) = p3
    P(X|ND) = (p1+p2)/2 - c (Delta1 + Delta2) /2
    P(X|NA) = (p1+p3)/2 - c (Delta1) /2
    P(X|DA) = (p1+p3)/2 - c (Delta2) /2

    #WITH ERROR, this becomes
    P(X|p, E) = E(1-p) + (1-E)p

    """

    #init, could be made more general
    n_states = 6
    alpha0 = np.array([1/n_states] * n_states)
    trans_mat = np.zeros((n_states,n_states)) + 1e-2
    np.fill_diagonal(trans_mat, 1 - 1e-2 * (n_states -1))
    assert np.allclose(np.sum(trans_mat, 1), 1)

    ll = -np.inf
    new_trans_mat = np.zeros_like(trans_mat)
    error, cont, cont_prev = E0, C0, C0


    libs = list(pd.unique(data.lib))
    n_libs = len(libs)
    cont, cont_prev  = [C0 for _ in libs],  [C0 for _ in libs]

    if binned:
        H2E, bins = get_H2E(data, bin_size, posvar='map')
        H2E_flat = get_H2E_flat(H2E)
        H2E_flat_split = split_by_lib(H2E_flat, data.lib, libs)

    if estimate_c and estimate_c_em:
        N0 = np.array(data.ref + data.alt)
        O0 = np.array(data.alt)

        #set up numpy array to make stuff work with numba
        P0_cont = np.array(data.AFR)
        P0 = np.array([data.NEA, data.DEN, data.AFR, 
                       (data.NEA + data.AFR) / 2, 
                       (data.DEN + data.AFR)/2 , 
                       (data.DEN + data.NEA) / 2]).T
        P = P0 * (1-E0) + E0 * (1-P0)
        P_cont = np.array(P0_cont * (1-E0) + E0 * (1-P0_cont))

        N = split_by_lib(N0, data.lib, libs)
        O = split_by_lib(O0, data.lib, libs)
        P_cont = split_by_lib(P_cont, data.lib, libs)
        P = split_by_lib(P, data.lib, libs)
        P = [p.T for p in P] #P[lib][state][obs]

        # po_given_rz[L x S x O x N] gives the probability of
        # the O-th observation of Library L, given it has been
        # emitted by state S and R/N reads are contaminants
        #init_global_po(O, N, P_cont, P)


        ER = [None for _ in libs]


    if binned:
        emission_mat = get_emission_mat(data, cont, libs, error, bad_snp_cutoff)
        emissions = get_emission_prob(emission_mat, bins, H2E_flat, bad_snp_cutoff)
    else:
        emissions = get_emission_mat_p0(O0, N0, P0_cont, P0, cont, libs, error,
                                        bad_snp_cutoff)


    for it in range(max_iter):
        alpha, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
        beta = bwd_algorithm(emissions, n=n, trans_mat=trans_mat)
        gamma = [a * b for (a, b) in zip(alpha, beta)]
        for g in gamma:
            assert np.allclose(np.sum(g, 1), 1)
        ll, old_ll =  np.sum([np.sum(np.log(n_)) for n_ in n]), ll


        #update transition
        for i in range(n_states):
            for j in range(n_states):
                new_trans_mat[i, j] = 0
                for a, b, e, n_ in zip(alpha, beta, emissions, n):
                    new_trans_mat[i, j] += np.sum(a[:-1, i] *
                                                  trans_mat[i, j] *
                                                  b[1:, j] *
                                                  e[:,j] /
                                                  n_[1:])

        gamma_sum = np.sum([np.sum(g[:-1], 0) for g in gamma], 0)
        new_trans_mat /= gamma_sum[:, np.newaxis]
        assert np.allclose(np.sum(new_trans_mat,1), 1)

        #update emissions
        if estimate_c and estimate_c_em:
            cont_prev = cont.copy()
            for i, lib in enumerate(libs):
                ll_local = ll
                while True:
                    if False:
                        if binned:
                            ER[i] = np.array([er_given_zoc(cont[i], O[i], N[i],
                                                           P_cont[i], pz,
                                                           i, s)
                                for (s, pz) in enumerate(P[i])]).T
                            cont[i], cont_prev[i] = maximize_c(ER[i], gamma, N[i], H2E_flat_split[i]), cont[i]
                        else:
                            ER[i] = np.array([er_given_zoc(cont[i], O[i], N[i],
                                                           P_cont[i], pz,
                                                           i, s)
                                for (s, pz) in enumerate(P[i])]).T
                            cont[i], cont_prev[i] = np.sum(ER[i] * gamma[i][1:]) / np.sum(N[i]), cont[i]


                    #numeric minimization
                    def get_po_given_zc_all(c):
                        cont[i] = c
                        po_given_zc = np.array([get_po_given_zc(c, O[i], N[i], P_cont[i],
                                                                pz, bad_snp_cutoff) 
                                       for pz in P[i]]).T
                        prob = np.sum(np.log(np.sum(po_given_zc * gamma[i][1:],1)))
                        print("[%s]minimizing c:" % lib, c, prob)
                        return -prob
                    OO =  minimize(get_po_given_zc_all, [cont[i]], bounds=[(0., 1)])

                    print(it, "update c", lib, i, len(libs), N[i].shape[0], cont[i],
                          cont_prev[i])

                    if citer_likelihood:
                        if binned:
                            emission_mat = get_emission_mat(data, cont, libs, error)
                            emissions = get_emission_prob(emission_mat, bins, H2E_flat, bad_snp_cutoff)
                        else:
                            emissions = get_emission_mat_p0(O0, N0, P0_cont, P0, cont, 
                                                            libs, error, bad_snp_cutoff)
                        _, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
                        ll_local, prev_ll_local = np.sum([np.sum(np.log(n_)) for n_ in n]), ll_local
                        print("citer_ll: ", ll_local, ll_local - prev_ll_local)
                        if ll_local - prev_ll_local < 1: break
                    else:
                        break #single iter of R maximization
                        if abs(cont[i] - cont_prev[i]) < 0.0001: break

            if binned:
                emission_mat = get_emission_mat(data, cont, libs, error)
                emissions = get_emission_prob(emission_mat, bins, H2E_flat, bad_snp_cutoff)
            else:
                            emissions = get_emission_mat_p0(O0, N0, P0_cont, P0, cont, 
                                                            libs, error, bad_snp_cutoff)

        # numeric brute force
        elif estimate_c and not estimate_c_em: 
            for i, lib in enumerate(libs):
                def ll_cont(c):
                    cont[i] = c
                    emission_mat = get_emission_mat(data, cont, libs, error)
                    emissions = get_emission_prob(emission_mat, bins, H2E_flat, bad_snp_cutoff)
                    ll = calc_ll(alpha0, trans_mat, emissions)
                    if(np.isnan(ll)): 
                        ll = -np.inf
                    print("[%s]minimizing c:" % lib, c, -ll)
                    return -ll
                #O =  fminbound(approx_cont, 0., 1, xtol=1e-5, disp=3)
                OO =  minimize(ll_cont, [cont[i]], bounds=[(0., 1)], tol=1e-5)
                #new_ll = calc_ll(alpha0, trans_mat, emissions)
                #print(ll, new_ll)
                print(O)
                

        #update stationary probabilities
        alpha0 = np.linalg.matrix_power(new_trans_mat, 10000)[0]

        #probs = np.sum((A[has_obs]/N[has_obs])[:,np.newaxis] * gamma[1:][has_obs],0) / np.sum(gamma[1:][has_obs],0)

        n_steps = emissions[0].shape[0]
        print("iter %d [%d/%d]: %s -> %s"% (it, n_steps, n_states , ll, ll - old_ll))
#        print(probs)
        print(cont)
        print(alpha0)
        print(new_trans_mat)
        print("---")

        if ll - old_ll  < 1e-1:
            break


        trans_mat = np.copy(new_trans_mat)

    return new_trans_mat, alpha0, gamma, ll, data


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
            aux_mat = ll[i-1] + log_t + log_e[i]
        ll[i] = np.max(aux_mat, 1)
        backtrack[i] = np.argmax(aux_mat, 1)

    path = np.empty(n_steps)
    cursor = np.argmax(ll[-1])
    for i in reversed(range(n_steps)):
        cursor = backtrack[i, cursor]
        path[i] = cursor


    return path #, ll, backtrack


def run_hmm(infile, bedfile=None, 
            outfile=None, out_par=None, out_tmat=None,
                     E0 = 1e-2, C0=1e-1,
                     split_lib=False,
                     **kwargs
                     ):
    """ run baum welch to find introgressed tracts
    infile formatting:
        - chrom, map: genomic coordinates
        - tref, talt: number of ref, alt reads
        - f_nea, f_afr: frequency of african and neandertal
    """
    try:
        data=pd.read_csv(infile)
        data=data[data.ref+data.alt>0]
        data=data.dropna()
        q = np.quantile(data.ref + data.alt, .999)
        data=data[data.ref+data.alt <=q]
    except ValueError:
        data = infile


    if "lib" not in data or (not split_lib):
        data =  data.groupby(("chrom" ,"pos", "NEA", "DEN", "AFR", "PAN"),
                      as_index=False).agg({"ref" : sum, "alt" : sum})
        q = np.quantile(data.ref + data.alt, .999)
        data=data[data.ref+data.alt <=q]
        data["lib"] = "lib0"

    if bedfile is not None:
        bed = pd.read_table(bedfile, header=None)[[0, 2]]
        bed.columns = ["chrom", "pos"]
        #data = bed.merge(data, on=["chrom", "pos"], how="left")

    t,alpha,  gamma, ll, data = baum_welch(data, E0, C0, **kwargs)
    gamma_flat = np.vstack([g[1:] for g in gamma])
    

    #final_emission_mat = get_emission_mat(data, C0, E0)
    #emissions = get_emission_prob(final_emission_mat, H2E, bad_snp_cutoff)

    gamma_df = pd.DataFrame(gamma_flat, columns = ("NN", "DD", "AA", "AN", "AD", "ND"))
    df = pd.concat((data, gamma_df), 1)


    if outfile is not None:
        df.to_csv(outfile)

    if out_tmat is not None:
        pd.DataFrame(t).to_csv(out_tmat, index=False, header = False)

    if out_par is not None:
        pinf =np.linalg.matrix_power(t, 10000)[0]
        with open(out_par, "w"):
            print("state", "p_init", "p_stat", sep=",")
            for i in range(n_states):
                print("S%d" %i, alpha[i], pinf[i], sep=",")

    return t, alpha, gamma,  data, df

#@jit(nopython=True)
def get_probs_p0(c, E, P0, P0_cont):
    probs = P0_cont[:, np.newaxis] * c + (1.-c) * P0
    probs = (1.-E) * probs + E * (1.-probs)

    return probs.T


#
def get_probs_pd(data, c, E):
    pNN = data.NEA - c * data.dN
    pDD = data.DEN - c * data.dD
    pAA = data.AFR 
    #pPP = data.PAN
    pNA = (data.NEA + data.AFR)/2 - c * data.dN/2
    pDA = (data.DEN + data.AFR)/2 - c * data.dD/2
    pND = (data.DEN + data.NEA)/2 - c * (data.dD+data.dN)/2
    #pNP = (data.NEA + data.PAN)/2 - c * (data.dP+data.dN)/2
    #pDP = (data.DEN + data.PAN)/2 - c * (data.dD+data.dP)/2
    #pAP = (data.PAN + data.AFR)/2 - c * (data.dP)/2

    #probs = (pNN, pDD, pAA, pPP, pNA, pDA, pAP, pND, pNP, pDP)
    probs = (pNN, pDD, pAA, pNA, pDA, pND)
    probs = [(1-E) * p + E * (1-p) for p in probs]

    return probs

#@jit(nopython=True)
def get_emission_mat(data, cont, libs, E=1e-2, bad_snp_cutoff=1e-16):
    """ get matrix[n_obs x n_states] giving the emission
    prob for each observation given each state
    """
    c = np.empty(data.shape[0])
    for i, lib in enumerate(libs):
        c[data.lib==lib] = cont[i]

    probs = get_probs(data, c, E)
    p2 = np.array([pbinom(data.alt, data.ref+data.alt, p) for p in probs]).T
    bad_snps = np.sum(p2, 1) < bad_snp_cutoff
    p2[bad_snps] = bad_snp_cutoff
    return p2

#@jit(nopython=True)
def get_emission_mat_p0(O, N, P0_cont, P0, cont, libs, E=1e-2, bad_snp_cutoff=1e-16):
    """ get matrix[n_obs x n_states] giving the emission
    prob for each observation given each state
    """
    probs = get_probs_p0(cont[0], E, P0, P0_cont)
    #ignoring binom coeff cause numba doesnt like them and they're constants
    #p2 = np.array([pbinom(O, N, p) for p in probs]).T
    p2 = np.array([p ** O * (1-p) ** (N-O) for p in probs]).T
    p2[np.where(p2!=p2)] = 1 #missing data has emission prob 1
    bad_snps = np.sum(p2, 1) < bad_snp_cutoff
    p2[bad_snps] = bad_snp_cutoff
    return [p2]

@jit(nopython=True)
def get_emission_prob(emission_mat, bins, H2E_flat, bad_snp_cutoff=1e-4):
    """ get list[n_chrom] of matrices[n_steps x n_states] giving the emission
    prob for each chrom, step and state
    aggregates the emission mat
    """
    n_obs, n_states = emission_mat.shape
    emission_probs = []
    for chrom_bins in bins:
        emission_probs.append(np.ones((len(chrom_bins), n_states)))

    #for (chrom, bin_), e in zip(H2E_flat, emission_mat):
    for i in range(n_obs):
        chrom, bin_ = H2E_flat[i]
        emission_probs[chrom][bin_] *= emission_mat[i]

    n_bad_snps = 0.
    for e in emission_probs:
        bad_snps = np.sum(e, 1) < bad_snp_cutoff
        n_bad_snps += np.sum(bad_snps)
        e[bad_snps] = bad_snp_cutoff
    print("bad_snps ", n_bad_snps)

    return emission_probs

def get_H2E_flat(H2E):
    """
    flat version of H2E:
    - each row is an observation
    - col 0 is chrom (chain)
    - col 1 is bin of that observation (step)

    ==> we know for each observation its hidden state
    """
    x = []
    for chrom, _ in enumerate(H2E):
        for bin_, obs_list in H2E[chrom].items():
            for obs in obs_list:
                x.append((chrom, bin_))
    return np.array(x)

def get_H2E(data, bin_size, posvar='map'):
    """
    get dict[chrom][bin] : [obs1, obs2, ...]
    """
    chroms = np.unique(data.chrom)
    H2E = []
    bin_pos = []
    prev_obs = 0
    
    for i, chrom in enumerate(chroms):
        D = defaultdict(list)

        pos = np.array(data.loc[data.chrom == chrom, posvar])
        min_ = int(np.floor(pos[0] / bin_size) * bin_size)
        max_ = int(np.ceil(np.array(pos)[-1] / bin_size) * bin_size)
        bins = np.arange(min_, max_, bin_size)
        dig = np.digitize(pos, bins, right=False) - 1

        for j, d in enumerate(dig):
            D[d].append(j + prev_obs) 

        H2E.append(D)
        bin_pos.append(bins)
        prev_obs += len(dig)

    assert prev_obs == data.shape[0]
    return H2E, bin_pos

 
def maximize_c2(ER, gamma, N, H2E):                                
    S = 0                                                                              
    print(ER.shape, gamma.shape, N.shape, H2E.shape)
    n_chrom = len(gamma)
    for chrom in range(n_chrom):                                                                
        for i, g in enumerate(gamma[chrom]):                                           
            S += np.sum(  np.sum(ER[H2E[chrom][i]], 0) * gamma[chrom][i])              
    return S / np.sum(N)                                                               


@jit(nopython=True)
def maximize_c(ER, gamma, N, H2E_flat):
    S = 0.#np.zeros_like(ER)
    #for i, (chrom, bin_) in enumerate(H2E_flat):
    n_snp = H2E_flat.shape[0]
    for i in range(n_snp):
        chrom, bin_ = H2E_flat[i]
        S += np.sum(ER[i] * gamma[chrom][bin_])
    print(S, N)
    return S / np.sum(N)

