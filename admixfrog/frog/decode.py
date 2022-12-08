import logging
from numba import njit, uint8
from numba.typed import List
from numba.types import UniTuple
import pandas as pd
from collections import Counter
import numpy as np
from random import random, seed


@njit
def decode_runs(seq, homo_ids, het_ids, roh_ids):
    n_homo = len(homo_ids)
    n_het = len(het_ids)
    n_roh = len(roh_ids)

    n_raw_states = n_homo


    #diploid to haploid dict
    D2H = np.zeros((n_homo + n_het + n_roh, 2), dtype="i")
    for i in range(n_homo):
        D2H[i] = homo_ids[i]
    for i in range(n_het):
        D2H[i + n_homo] = het_ids[i]
        if het_ids[i][0] > n_raw_states:
            n_raw_states = het_ids[i][0] + 1
        if het_ids[i][1] > n_raw_states:
            n_raw_states = het_ids[i][1] + 1
    for i in range(n_roh):
        D2H[i + n_homo + n_het] = roh_ids[i]

    # init list of ints, numba needs typing
    runs = [[(i, i, i) for i in range(0)] for i in range(n_raw_states)]  # run lengths

    for i in range(len(seq)):

        #init for very first bin: 
        if i == 0:
            r1, r2 = D2H[seq[0]] #init ID of the two runs
            l1, l2 = 1, 1  #lengths start at 1
        else:
            c1, c2 = D2H[seq[i]] #the "continuing states"

            #match, extending both 
            if (c1 == r1 and c2 == r2) or (c1 == r2 and c2 == r1):
                l1, l2 = l1 + 1, l2 + 1

            # one strand ends
            elif (
                (c1 == r1 and c2 != r2)
                or (c1 == r2 and c2 != r1)
                or (c2 == r1 and c1 != r2)
                or (c2 == r2 and c1 != r1)
            ):
                #first case: we were homozygous before. Therefore, we are now
                #het. We end one strand randomly and write
                if r1 == r2:  
                    if random() < 0.5: #ending r1
                        runs[r1].append((l1, i - l1, i))
                        r1 = c2 if r1 == c1 else c1  # new run with non-matching
                        l1, l2 = 1, l2 + 1
                    else: #ending r2
                        runs[r2].append((l2, i - l2, i)) 
                        r2 = c2 if r2 == c1 else c1  # new run with non-matching
                        l1, l2 = l1 + 1, 1
                else:  # we were heterozygous before
                    """cases are het to homo or AB -> AC or AB -> CB """
                    if r1 == c1 and r2 != c2: #AB -> AC. End r2, keep r1
                        runs[r2].append((l2, i - l2, i))
                        r1, r2, l1, l2 = r1, c2, l1 + 1, 1
                    elif r1 == c2 and r2 != c1: #AB -> CA. End r2, switch r1, r2
                        runs[r2].append((l2, i - l2, i))
                        r1, r2, l1, l2 = c2, c1, l1 + 1, 1
                    elif r1 != c1 and r2 == c2: #AB -> CB. End r1, keep r2
                        runs[r1].append((l1, i - l1, i))
                        r1, r2, l1, l2 = c1, r2, 1, l2 + 1
                    elif r1 != c2 and r2 == c1: #AB -> BC. End r1, switch r1, r2
                        runs[r1].append((l1, i - l1, i))
                        r1, r2, l1, l2 = c2, c1, 1, l2 + 1
                    else:
                        raise ValueError("not handled decode case")

            # both strands end, write results and reset to 0
            elif c1 != r1 and c1 != r2 and c2 != r1 and c2 != r2:
                runs[r1].append((l1, i - l1, i))
                runs[r2].append((l2, i - l2, i))
                r1, r2, l1, l2 = c1, c2, 1, 1
            else:
                raise ValueError("not handled decode case")

    #finalize
    runs[r1].append((l1, i - l1, i))
    runs[r2].append((l2, i - l2, i))
    return runs


@njit
def decode_runs_single(seq, n_states, keep_loc=False):
    runs = [[(i, i, i) for i in range(0)] for i in range(n_states)]  # run lengths
    for i in range(len(seq)):
        # print(r1, r2, l1, l2)
        if i == 0:
            r, l = seq[0], 1
        else:
            if seq[i] == r:
                l += 1
            else:
                runs[r].append((l, i - l, i))
                r, l = seq[i], 1

    runs[r].append((l, i - 1, i))

    return runs


@njit
def nb_choice(n, p):
    u = random()
    for i in range(n):
        if u <= p[i]:
            return i
        u -= p[i]
    return i


@njit
def post_trans(trans, emissions, beta, beta_prev, n):
    return beta / beta_prev / n * trans * emissions


@njit
def pred_sims_rep(
    trans,
    emissions,
    beta,
    alpha0,
    n,
    homo_ids,
    het_ids,
    roh_ids,
    decode=True,
    keep_loc=False,
):
    n_steps, n_states = emissions.shape

    seq = np.zeros(n_steps, dtype=np.int64)
    for i in range(n_steps):
        if i == 0:
            state = nb_choice(n_states, p=alpha0)
        else:
            p = post_trans(
                trans=trans[state],
                emissions=emissions[i],
                beta=beta[i],
                beta_prev=beta[i - 1, state],
                n=n[i],
            )
            state = nb_choice(n_states, p)
        seq[i] = state
    if decode:
        runs = decode_runs(seq, homo_ids, het_ids, roh_ids)
    else:
        runs = decode_runs_single(seq, n_states=len(homo_ids) + len(het_ids) +
                                  len(roh_ids))
    return runs


def pred_sims_single(
    trans,
    emissions,
    beta,
    alpha0,
    n,
    states,
    n_sims=100,
    decode=True,
    keep_loc=False,
):
    sims = []

            
    homo_list = List(states.homo_ids) if len(states.homo_ids) else List.empty_list(uint8)
    het_list = List(states.het_ids) if len(states.het_ids) else List.empty_list(UniTuple(uint8, 2))
    roh_list = List(states.roh_ids) if len(states.roh_ids) else List.empty_list(uint8)
    for it in range(n_sims):
        runs = pred_sims_rep(
            trans,
            emissions,
            beta,
            alpha0,
            n,
            homo_list, het_list, roh_list,
            decode=decode,
            keep_loc=keep_loc
        )


        for i, run in enumerate(runs):
            if keep_loc:
                df = pd.DataFrame(run, columns=("len", "start", "end"))
            else:
                df = pd.DataFrame(
                    Counter(r[0] for r in run).items(), columns=("len", "n")
                )
            df["state"] = i
            df["it"] = it
            sims.append(df)
    sims = pd.concat(sims)
    return sims


def pred_sims(
    trans,
    emissions,
    beta,
    alpha0,
    n,
    states,
    n_sims=100,
    decode=True,
    keep_loc=False,
):
    """simulate runs through the model using posterior parameter.

    uses the algorithm of Nielsen, Skov et al. to generate track-length
    distribution.

    Parameters
    =====
    trans: transition matrix
    emissions: list of emission matrices
    beta: list of result of bwd-algorithm
    alpha0: initial probability
    n : list of normalizing factors from fwd-algorithm
    states : state object
    n_sims: number of reps
    decode: whether diploid states are decoded into haploid ones
    keep_loc: whether location info should be kept

    """
    output = []
    for i, (e, b, n_) in enumerate(zip(emissions, beta, n)):
        df = pred_sims_single(
            trans, e, b, alpha0, n_, states, n_sims, decode, keep_loc
        )
        df["chrom"] = i
        output.append(df)
        logging.info("Posterior Simulating chromosome %s" % i)
    return pd.concat(output)


def resampling_pars(tbl):
    in_state = tbl[["state", "it", "len"]].groupby(["state", "it"]).sum()
    tot = tbl[["it", "len"]].groupby(["it"]).sum()
    Z = in_state / tot
    sds = Z.groupby("state").std()
    means = Z.groupby("state").mean()

    sds.columns = ["sd"]
    means.columns = ["mean"]
    data = means.join(sds)
    data["lower"] = data["mean"] - 1.96 * data["sd"]
    data["upper"] = data["mean"] + 1.96 * data["sd"]

    return data
