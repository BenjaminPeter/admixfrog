from numba import njit
import pandas as pd
from collections import Counter
import numpy as np
from numpy.random import choice
from random import random


@njit  # ('i8(i8, i8)')
def get_hap_from_diploid(n_states, n_homo, est_inbreeding=False):
    n_states = n_states + n_homo if est_inbreeding else n_states
    het2homo = np.zeros((n_states, 2), np.int8)
    for s in range(n_homo):
        het2homo[s] = s
    for s1 in range(n_homo):
        for s2 in range(s1 + 1, n_homo):
            s += 1
            het2homo[s] = [s1, s2]
    if est_inbreeding:
        for i in range(n_homo):
            s += 1
            het2homo[s] = -i
        
    return het2homo


@njit
def decode_runs(seq, n_homo, n_het, est_inbreeding=False):
    D2H = get_hap_from_diploid(n_homo + n_het, n_homo, est_inbreeding)
    # print(D2H)

    # init list of ints, numba needs typing
    runs = [[i for i in range(0)] for i in range(n_homo)]
    for i in range(len(seq)):
        # print(r1, r2, l1, l2)
        if i == 0:
            r1, r2 = D2H[seq[0]]
            l1, l2 = 1, 1
        else:
            c1, c2 = D2H[seq[i]]
            # match
            if (c1 == r1 and c2 == r2) or (c1 == r2 and c2 == r1):
                # print(i, "Extending both")
                l1, l2 = l1 + 1, l2 + 1
            # one strand ends
            elif (
                (c1 == r1 and c2 != r2)
                or (c1 == r2 and c2 != r1)
                or (c2 == r1 and c1 != r2)
                or (c2 == r2 and c1 != r1)
            ):
                if r1 == r2:  # homo to het
                    if choice(2) == 0:
                        # print(i, "match r1==r2, A")
                        runs[r1].append(l1)
                        r1 = c2 if r1 == c1 else c1  # ne run with non-matching
                        l1, l2 = 1, l2 + 1
                    else:
                        # print(i, "match r1==r2, B")
                        runs[r2].append(l2)
                        r2 = c2 if r2 == c1 else c1  # ne run with non-matching
                        l1, l2 = l1 + 1, 1
                else:  # het to homo or AB -> AC or AB -> CB
                    if r1 == c1 and r2 != c2:
                        runs[r2].append(l2)
                        r1, r2, l1, l2 = c1, c2, l1 + 1, 1
                    elif r1 == c2 and r2 != c1:
                        runs[r2].append(l2)
                        r1, r2, l1, l2 = c2, c1, l1 + 1, 1
                    elif r2 == c2 and r1 != c2:
                        runs[r1].append(l1)
                        r1, r2, l1, l2 = c1, c2, 1, l2 + 1
                    elif r2 == c1 and r1 != c2:
                        runs[r1].append(l1)
                        r1, r2, l1, l2 = c2, c1, 1, l2 + 1
                    else:
                        print("not handled case", r1, r2, c1, c2, seq[i])

            # both strands end
            elif c1 != r1 and c1 != r2 and c2 != r1 and c2 != r2:
                # print(i, "change both", r1, r2, c1, c2)
                runs[r1].append(l1)
                runs[r2].append(l2)
                r1, r2, l1, l2 = c1, c2, 1, 1
            else:
                print("not handled case", r1, r2, c1, c2, seq[i])
    runs[r1].append(l1)
    runs[r2].append(l2)
    return runs


@njit
def decode_runs_single(seq, n_states):
    runs = [[i for i in range(0)] for i in range(n_states)]
    for i in range(len(seq)):
        # print(r1, r2, l1, l2)
        if i == 0:
            r, l = seq[0], 1
        else:
            if seq[i] == r:
                l += 1
            else:
                runs[r].append(l)
                r, l = seq[i], 1

    runs[r].append(l)
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
    x = beta / beta_prev / n * trans * emissions
    # if not np.isclose(np.sum(x), 1):
    #    raise ValueError("error in posterior transition")
    return x


@njit
def pred_sims_rep(trans, emissions, beta, alpha0, n, n_homo, decode=True,
                  est_inbreeding=False):
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
                n=n[i]
            )
            state = nb_choice(n_states, p)
            seq[i] = state
    if decode:
        runs = decode_runs(seq, n_homo, n_states - n_homo, est_inbreeding)
    else:
        runs = decode_runs_single(seq, n_states=np.max(seq)+1)
    return runs


def pred_sims_single(
    trans, emissions, beta, alpha0, n, n_homo, n_sims=100, decode=True,
    est_inbreeding=False
):
    sims = []
    for it in range(n_sims):
        runs = pred_sims_rep(trans, emissions, beta, alpha0, n, n_homo, decode,
                             est_inbreeding)
        for i, run in enumerate(runs):
            df = pd.DataFrame(Counter(run).items(), columns=("len", "n"))
            df["state"] = i
            df["it"] = it
            sims.append(df)
    return pd.concat(sims)


def pred_sims(trans, emissions, beta, alpha0, n, n_homo, n_sims=100,
              decode=True, est_inbreeding=False):
    """simulate runs through the model using posterior parameter.

    uses the algorithm of Nielsen, Skov et al. to generate track-length
    distribution. Could be extended to also give positions of stuff...

    """
    output = []
    for i, (e, b, n_) in enumerate(zip(emissions, beta, n)):
        df = pred_sims_single(trans, e, b, alpha0, n_, n_homo, n_sims, decode,
                              est_inbreeding)
        df["chrom"] = i
        output.append(df)
        print("Posterior Simulating chromosome %s" % i)
    return pd.concat(output)


def test():
    beta = np.ones((100000, 6))
    emissions = np.ones_like(beta)
    n = np.ones(100000)
    trans = np.zeros((6, 6)) + 1 / 6
    alpha0 = [.3, .4, .3, 0., 0., 0.]

    L = lambda runs, i: sum(runs[runs.state == i].len * runs[runs.state == i].n)
    return pred_sims_rep(trans, emissions, beta, alpha0, n, 3)


def pred_sims_slow(trans, emissions, beta, alpha0, n, n_homo, n_sims=100):
    n_steps, n_states = emissions.shape
    X = np.zeros((n_steps, n_states, n_states))
    X[:] = (beta * n[:, np.newaxis] * emissions)[:, :, np.newaxis]
    X = X * trans
    X[1:] /= beta[:-1, :, np.newaxis]

    seq = np.zeros(n_steps, dtype=int)
    sims = []
    for it in range(n_sims):
        print("sim %s" % it)
        for i in range(n_steps):
            if i == 0:
                state = choice(n_states, p=alpha0)
            else:
                state = choice(n_states, p=X[i, state])
                seq[i] = state
        runs = decode_runs(seq, n_homo, n_states - n_homo)
        for i, run in enumerate(runs):
            df = pd.DataFrame(Counter(run).items(), columns=("len", "n"))
            df["state"] = i
            df["it"] = it
            sims.append(df)

    return pd.concat(sims)
