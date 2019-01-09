from scipy.stats import binom
from collections import defaultdict
import pandas as pd
import numpy as np
from pprint import pprint
from functools import lru_cache
import pdb

np.set_printoptions(suppress=True, precision=4)          
pbinom = binom.pmf

@lru_cache(maxsize=65536)
def pbinom_cache(*args, **kwargs):
    return binom.pmf(*args, **kwargs)

@lru_cache(maxsize=65536)
def pbinom_range(k_max, n, p):
    return binom.pmf(np.arange(k_max+1), n, p)
@lru_cache(maxsize=65536)
def pbinom_range_rev(k_start, k_max, n, p):
    return binom.pmf(k_start - np.arange(k_max+1), n, p)

ctr = 0
@lru_cache(maxsize=65536)
def po_given_rz_single(o, n, pc, pz):
    global ctr
    ctr += 1
    if ctr % 1000 == 0:
        print(ctr)
    s = [0] * (n+1)
    for r in range(n+1):
        s[r] = np.sum(pbinom_range(n, r, pc) * pbinom_range_rev(o, n, n-r, pz))
    return s

"""
Data: 
    either:
    S: -sites segregating in both Africans and Neandertals.
    F: -sites fixed different bewteen AFrican/Neandertals
        -These sites must be either due to gene flow between Africans / Neandertal
        - or older than Human-Neandertal split
    - # reads supporting afr / Nea allele
    A: Afr
    B: Non-Afr allele

Hidden state layer 1: 
    X: Het 
    N: Homozygous Neandertal 
    H: Homo Human /Afr allele

Hidden state layer 2:
    ILS: due to incomplete lineage sorting
    T_i due to introgression at time i


Params:
    - c: contamination rate
    - eps: error rate
    - N_anc : size of NEa-human ancestor
    - t_anc : time of split

Emission:
    P(A| A+B, X, F) ~ Binom(A+B, c [ 1- eps] + (1-c) [0.5]      |p= (1+c)/2 - c eps
    P(A| A+B, N, F) ~ Binom(A+B, c [ 1- eps] + (1-c) eps)       |p= c + eps - 2eps c
    P(A| A+B, H, F) ~ Binom(A+B, c [ 1- eps] + (1-c) (1 - eps)) |p= 1 - eps

Transitions:
    P(X | H) = r t_geneflow 
    P(X | X) = r t_geneflow
    P(N | X) = 1 - r t_geneflow - 
    P(H | H) = 1 - r t_genflow - r  
"""

def fwd_step(alpha_prev, E, trans_mat):
    alpha_new = (alpha_prev @ trans_mat) * E
    n = sum(alpha_new)
    return alpha_new / n, n

def fwd_algorithm(alpha0, emissions, trans_mat,  **kwargs):
    """
    calculate P(X_t | o_[1..t], a0)
    =P(X_t , o_[1..t], a0 | o_[1..t])
    """

    alpha, n = dict(), dict()

    for chrom, em in emissions.items():
        n_steps, n_states = em.shape
        alpha_i = np.empty((n_steps+1, n_states))
        n_i = np.empty((n_steps+1))

        alpha_i[0] = alpha0
        n_i[0] = np.sum(alpha0)

        for i, e in enumerate(em):
            alpha_i[i + 1], n_i[i + 1] = fwd_step(alpha_i[i], e, trans_mat, **kwargs)
            if n_i[i+1] < 1e-10:
                pass
                #print(chrom, i, e)

        alpha[chrom] = alpha_i
        n[chrom] = n_i

    return alpha, n

def bwd_step(beta_next, E, trans_mat, n):
    beta = (trans_mat * E) @ beta_next
    return beta / n

def bwd_algorithm(emissions, trans_mat, n, **kwargs):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n])
    """

    beta = dict()
    for chrom, em in emissions.items():
        n_steps, n_states = em.shape

        beta_i = np.empty((n_steps + 1, n_states))
        beta_i[n_steps] = 1

        for i, e in zip(range(n_steps-1, -1, -1),    reversed(em)):
            beta_i[i] = bwd_step(beta_i[i+1], e, trans_mat, n[chrom][i+1], **kwargs)
        beta[chrom] = beta_i
    return beta

def baum_welch_binom(coords, data, alpha0, trans_mat, p_init, max_iter=2000):
    n_states = len(p_init)
    n_obs = coords.shape[0]

    ll = -np.inf
    NT = np.zeros_like(trans_mat)
    probs = p_init
        
    data_dict = defaultdict(list)
    data = np.array(data)
    for chrom, a, b in zip(coords[:,0], data[:,0], data[:,1]):
        data_dict[chrom].append((a, a+b))
    for chrom in data_dict:
        data_dict[chrom] = np.array(data_dict[chrom]).T

    for it in range(max_iter):
        emissions = dict((chrom, np.array([pbinom(d[0], d[1], p) for p in probs]).T)
                        for chrom, d in data_dict.items())

        alpha, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
        beta = bwd_algorithm(emissions, n=n, trans_mat=trans_mat)

        assert alpha.keys() == beta.keys()

        gamma = dict( (chrom,alpha[chrom] * beta[chrom]) for chrom in alpha)
        for g in gamma.values():
            assert np.max(np.abs(np.sum(g, 1) - 1)) < 1e-8

#        NT = np.zeros((3, 3, 23))
        #update transition
        NT = np.zeros((3,3))
        for i in range(n_states):
            for j in range(n_states):
                NT[i, j] = 0
                for chrom in alpha.keys():
                    NT[i, j] += np.sum(alpha[chrom][:-1, i] * trans_mat[i, j] * beta[chrom][1:, j] * emissions[chrom][:,j] / n[chrom][1:])

        gamma_sum = np.sum([np.sum(g[:-1], 0) for g in gamma.values()], 0)
        NT = NT / gamma_sum[:, np.newaxis]
        assert np.max(np.abs(np.sum(NT, 1)-1) < 1e-6)

        #update emission using basic binomial
        num, denom = 0, 0
        for chrom, (A, N) in data_dict.items():
            obs = N > 0
            num += np.sum((A[obs]/N[obs])[:, np.newaxis] * gamma[chrom][1:][obs], 0)
            denom += np.sum(gamma[chrom][1:][obs], 0)

        probs = num / denom

        #update stationary probabilities
        alpha0 = np.mean([g[0] for g in gamma.values()], 0) 

        #probs = np.sum((A[has_obs]/N[has_obs])[:,np.newaxis] * gamma[1:][has_obs],0) / np.sum(gamma[1:][has_obs],0)
        ll, old_ll =  np.sum([np.sum(np.log(i)) for i in n.values()]), ll


        print("iter %d [%d/%d]: %s"% (it, n_obs, n_states , ll - old_ll))
        print(probs)
        print(alpha0)
        print(NT)
        print("---")

        if ll - old_ll  < 5e-1:
            break


        trans_mat = np.copy(NT)

    return NT, alpha0, probs, gamma


def baum_welch_multinom(data, E0, C0, bin_size=1e4, bad_snp_cutoff=1e-4, max_iter=2000,
                        estimate_c = True):
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

    - need data structure with multiple observations per hidden state
    - only thing to estimate is c, the contamination proportion
    """

    n_states = 6

    #uniform init
    alpha0 = [1/6] * 6
    trans_mat = np.zeros((6,6)) + .05
    np.fill_diagonal(trans_mat, .75)

    ll = -np.inf
    new_trans_mat = np.zeros_like(trans_mat)

    error, cont, cont_prev = E0, C0, C0

    H2E, bins = get_H2E(data, bin_size)

    if estimate_c :
        N = np.array(data.ref + data.alt)
        O = np.array(data.alt)
        P_contamination = data.AFR * (1-E0) + E0 * (1-data.AFR)
        P = np.array([data.NEA, data.DEN, data.AFR, (data.NEA + data.AFR) / 2, 
             (data.DEN + data.AFR)/2 , (data.DEN + data.NEA) / 2]).T

        P = P * (1-E0) + E0 * (1-P)
        ER = np.array([er_given_zoc(C0, O, pz, N, P_contamination) for pz in P.T]).T

    emission_mat = get_emission_mat(data, cont, error)
    emissions = get_emission_prob(emission_mat, H2E, bad_snp_cutoff)

    for it in range(max_iter):


        alpha, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
        beta = bwd_algorithm(emissions, n=n, trans_mat=trans_mat)

        assert alpha.keys() == beta.keys()

        gamma = dict( (chrom,alpha[chrom] * beta[chrom]) for chrom in alpha)
        for g in gamma.values():
            assert np.max(np.abs(np.sum(g, 1) - 1)) < 1e-8

#        NT = np.zeros((3, 3, 23))
        #update transition
        for i in range(n_states):
            for j in range(n_states):
                new_trans_mat[i, j] = 0
                for chrom in alpha.keys():
                    new_trans_mat[i, j] += np.sum(alpha[chrom][:-1, i] *
                                                  trans_mat[i, j] *
                                                  beta[chrom][1:, j] *
                                                  emissions[chrom][:,j] /
                                                  n[chrom][1:])

        gamma_sum = np.sum([np.sum(g[:-1], 0) for g in gamma.values()], 0)
        new_trans_mat /= gamma_sum[:, np.newaxis]
        assert np.max(np.abs(np.sum(new_trans_mat, 1)-1) < 1e-6)

        #update emissions
        if estimate_c :
            while True:
                cont, cont_prev = maximize_c(ER, gamma, N, H2E, bad_snp_cutoff), cont
                print("update c", cont)
                ER = np.array([er_given_zoc(cont, O, pz, N, P_contamination) for pz in P.T]).T
                if abs(cont - cont_prev) < 0.01: break
            emission_mat = get_emission_mat(data, cont, error)
            emissions = get_emission_prob(emission_mat, H2E, bad_snp_cutoff)

        #update stationary probabilities
        alpha0 = np.mean([g[0] for g in gamma.values()], 0) 

        #probs = np.sum((A[has_obs]/N[has_obs])[:,np.newaxis] * gamma[1:][has_obs],0) / np.sum(gamma[1:][has_obs],0)
        ll, old_ll =  np.sum([np.sum(np.log(i)) for i in n.values()]), ll


        n_steps = emission_mat.shape[0]
        print("iter %d [%d/%d]: %s -> %s"% (it, n_steps, n_states , ll, ll - old_ll))
#        print(probs)
        print(cont)
        print(alpha0)
        print(new_trans_mat)
        print("---")

        if ll - old_ll  < 5e-1:
            break


        trans_mat = np.copy(new_trans_mat)

    return new_trans_mat, alpha0, gamma, bins, H2E


def er_given_zoc(c, o, pz, n, pc, po_given_rz = None):
    """R * P(R = r | Z = k, o, c, n) = P(o | R, Z) P(R|Z c,n) / P(o | Z, c, n)

    return a vector of length n+1, giving the probabilty that r reads are contaminant
    """

    #pr_given_zc = pbinom(R, n, c)
    #po_given_zc = pbinom(o, n, pz + c * (pc - pz))
    pr_given_zc = [pbinom_range(n_, n_, c) for n_ in n]
    po_given_zc = [pbinom_cache(o_, n_, pz_ + c* (pc_ - pz_)) 
                for (o_, n_, pc_, pz_) in zip(o, n, pc, pz)]


    # as this does not depend on any parameters that change during EM, it should
    # be  precomputed
    if po_given_rz is None:
        po_given_rz = [po_given_rz_single(o_, n_, pc_, pz_)
                       for (pc_, pz_, n_, o_) in zip(pc, pz, n, o)]
    er= [np.sum(np.arange(n_+1) * x*y/z) 
        for (n_, x, y, z) in zip(n, po_given_rz, pr_given_zc, po_given_zc)]
    
    return er

    
def prep_data(chroms, pos, A, B, bin_size=10000):
    A_arr, B_arr = [], []
    chrom_arr, pos_arr = [], []
    for chrom in range(1, 23):
        d1 = pos[chrom==chroms]
        bins = np.arange(np.min(d1), np.max(d1), bin_size)
        dig = np.digitize(d1, bins, right=False) - 1
        A_binned = np.zeros(np.max(dig)+1, "i")
        B_binned = np.zeros(np.max(dig)+1, "i")

        for a,b, x in zip(A[chrom==chroms], B[chrom==chroms], dig): 
            A_binned[x] += a           
            B_binned[x] += b           

        assert sum(A[chrom==chroms]) == sum(A_binned)
        assert sum(B[chrom==chroms]) == sum(B_binned)

        A_arr.append(A_binned)
        B_arr.append(B_binned)
        chrom_arr.append(np.array([chrom for _ in bins]))
        pos_arr.append(bins)


    data = np.vstack((np.hstack(A_arr), np.hstack(B_arr))).T
    coords = np.vstack((np.hstack(chrom_arr), np.hstack(pos_arr))).T
    assert data.shape == data.shape

    return coords, data

    
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

    
        
        
def run_hmm_multinom(infile, outfile=None, out_par=None, out_tmat=None,
                     E0 = 1e-2, C0=1e-1,
                     bin_size=10000,
                    bad_snp_cutoff = 1e-6,
                     max_iter=200,
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
    except ValueError:
        data = infile
    data["dN"] = data.NEA - data.AFR
    data["dD"] = data.DEN - data.AFR

    t,alpha,  gamma, bins, H2E = baum_welch_multinom(data, E0, C0, bin_size, bad_snp_cutoff,
                                          max_iter)
    gamma_flat = np.vstack([g[1:] for g in gamma.values()])
    

    final_emission_mat = get_emission_mat(data, C0, E0)
    emissions = get_emission_prob(final_emission_mat, H2E, bad_snp_cutoff)

    vb_dict = dict((c,viterbi_single_obs(alpha, t, e)) for (c,e) in emissions.items())
    vb_df = pd.DataFrame(np.hstack(v for v in vb_dict.values()), columns=("viterbi",))
    vb_df.viterbi = vb_df.viterbi.astype(int)

    data["bin"] = np.nan
    for chrom, b in bins.items():
            data.loc[data.chrom==chrom, 'bin'] = np.digitize(data.pos[data.chrom==chrom], b)
    data.bin = data.bin.astype(int)


    coords =  [(chrom, *p) for (chrom, pos) in bins.items() for p in enumerate(pos)]
    coords = np.array(coords, np.int)
    coords = pd.DataFrame(coords, columns = ("chrom", "bin", "binpos"))
    gamma_df = pd.DataFrame(gamma_flat, columns = ("NN", "DD", "AA", "AN", "AD", "ND"))
    df = pd.concat((coords, gamma_df, vb_df), 1)


    data = data.merge(df)


    #reorder crap
    #o = np.argsort(probs)
    #t = t[o,:][:,o] 
    #probs = probs[o]
    #gamma_flat = gamma_float[:, o]


    if outfile is not None:
        d = pd.DataFrame(d)
        d.columns = "ref", "alt"
        df_gamma = pd.DataFrame(gamma_flat)
        df_gamma.columns = ("S%d"%i for i in range(gamma.shape[1]))
        coords = pd.DataFrame(coords)
        coords.columns = ("chrom", "map")

        q = pd.concat((coords, d, df_gamma), axis=1)
        q["obs"] = has_obs
        q.to_csv(outfile)

    if out_tmat is not None:
        pd.DataFrame(t).to_csv(out_tmat, index=False, header = False)

    if out_par is not None:
        pinf =np.linalg.matrix_power(t, 10000)[0]
        with open(out_par, "w"):
            print("state", "p_init", "p_stat", sep=",")
            for i in range(n_states):
                print("S%d" %i, alpha[i], pinf[i], sep=",")

    return t, alpha, gamma,  data, df

    

def run_hmm(infile, outfile, n_states, out_par, out_tmat,
            binsize=10000):
    """ run baum welch to find introgressed tracts
    infile formatting:
        - chrom, pos: genomic coordinates
        - tref, talt: number of ref, alt reads
        - f_nea, f_afr: frequency of african and neandertal
    """
    data=pd.read_csv(infile)
    data = data[data.comparison == "Vindija33.19"]
    data = data[np.abs(data.f_nea - data.f_afr) == 1]
    data = data[np.abs(data.f_other - data.f_afr) == 1]
    A = np.array(data.f_nea * data.tref + data.f_afr * data.talt, dtype="i")
    B = np.array(data.f_afr * data.tref + data.f_nea * data.talt, dtype="i")
    coords, d = prep_data(data.chrom, data.pos, A, B, bin_size)
    has_obs = np.sum(d,1) > 0
    t,alpha, probs, gamma = baum_welch_binom(coords, d, alpha0, trans_mat, 
                                       p_init=emission_p,
                                       max_iter=200)


    #reorder crap
    o = np.argsort(probs)
    t = t[o,:][:,o] 
    probs = probs[o]
    gamma_flat = np.vstack([g[1:] for g in gamma.values()])[:, o]


    if outfile is not None:
        d = pd.DataFrame(d)
        d.columns = "ref", "alt"
        df_gamma = pd.DataFrame(gamma_flat)
        df_gamma.columns = ("S%d"%i for i in range(gamma.shape[1]))
        coords = pd.DataFrame(coords)
        coords.columns = ("chrom", "pos")

        q = pd.concat((coords, d, df_gamma), axis=1)
        q["obs"] = has_obs
        q.to_csv(outfile)

    if out_tmat is not None:
        pd.DataFrame(t).to_csv(out_tmat, index=False, header = False)

    if out_par is not None:
        pinf =np.linalg.matrix_power(t, 10000)[0]
        with open(out_par, "w"):
            print("state", "p_init", "p_stat", sep=",")
            for i in range(n_states):
                print("S%d" %i, alpha[i], pinf[i], sep=",")

    return t, alpha, probs, gamma_flat

def get_probs(data, c, E):
    pNN = data.NEA - c * data.dN
    pDD = data.DEN - c * data.dD
    pAA = data.AFR 
    pNA = (data.NEA + data.AFR)/2 - c * data.dN/2
    pDA = (data.DEN + data.AFR)/2 - c * data.dD/2
    pND = (data.DEN + data.NEA)/2 - c * (data.dD+data.dN)/2

    probs = (pNN, pDD, pAA, pNA, pDA, pND)
    probs = [(1-E) * p + E * (1-p) for p in probs]

    return probs

def get_emission_mat(data, c, E=1e-2):
    if "lib" in data and isinstance(c, dict):
        c = [c[l] for l in data.lib]

    probs = get_probs(data, c, E)

    p2 = np.array([pbinom(data.alt, data.ref+data.alt, p) for p in probs]).T
    return p2

def get_emission_tensor(data, E=1e-2):
     [get_emission_mat(data, c/100, E0) for c in range(100)]



def get_H2E(data, bin_size, posvar='map'):
    """
    get dict[chrom, pos, hidden_state] : [obs1, obs2, ...]
    """
    chroms = np.unique(data.chrom)
    H2E = defaultdict(lambda : defaultdict(list))
    bin_pos = dict()
    prev_obs = 0
    
    for chrom in chroms:
        pos = np.array(data.loc[data.chrom == chrom, posvar])
        min_ = int(np.floor(pos[0] / bin_size) * bin_size)
        max_ = int(np.ceil(np.array(pos)[-1] / bin_size) * bin_size)
        bins = np.arange(min_, max_, bin_size)
        dig = np.digitize(pos, bins, right=False) - 1
        for i, d in enumerate(dig):
            H2E[chrom][d].append(i + prev_obs) 
        bin_pos[chrom] = bins
        prev_obs += len(dig)

    assert prev_obs == data.shape[0]
    return H2E, bin_pos

def get_emission_prob(emission_mat, H2E, bad_snp_cutoff=1e-4):
    n_obs, n_states = emission_mat.shape
    emission_probs = dict()

    for chrom, steps in H2E.items():
        d = np.ones((max(steps)+1, n_states))
        for step, snps in steps.items():
            d[step] = np.prod(emission_mat[snps], 0)  

        #set bins with very low prob to missing
        bad_snps = np.sum(d, 1) < bad_snp_cutoff
        d[bad_snps] = bad_snp_cutoff
        emission_probs[chrom] = d

    return emission_probs

def maximize_c(ER, gamma, N, H2E, bad_snp_cutoff=1e-4):
    S = 0
    for chrom in gamma:
        for i, g in enumerate(gamma[chrom]):
            S += np.sum(  np.sum(ER[H2E[chrom][i]], 0) * gamma[chrom][i])
    return S / np.sum(N)

