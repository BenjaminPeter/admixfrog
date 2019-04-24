from numba import njit
import numpy as np
from .log import log_


@njit
def calc_ll(alpha0, trans_mat, emissions):
    """likelihood using forward algorithm"""
    _, n = fwd_algorithm(alpha0, emissions, trans_mat=trans_mat)
    return np.sum([np.sum(np.log(n_)) for n_ in n])


@njit
def fwd_step(alpha_prev, E, trans_mat):
    alpha_new = (alpha_prev @ trans_mat) * E
    n = np.sum(alpha_new)
    return alpha_new / n, n


@njit
def fwd_algorithm_single_obs(alpha0, emission, trans_mat):
    """
    calculate P(X_t | o_[1..t], a0)
    =P(X_t , o_[1..t], a0 | o_[1..t])
    """
    n_steps, n_states = emission.shape
    alpha = np.empty((n_steps, n_states))
    n = np.empty((n_steps))
    alpha[0] = alpha0
    n[0] = np.sum(alpha0)
    for i in range(n_steps):
        if i == 0:
            alpha[i], n[i] = fwd_step(alpha0, emission[i], trans_mat)
        else:
            alpha[i], n[i] = fwd_step(alpha[i - 1], emission[i], trans_mat)
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


@njit
def bwd_step(beta_next, E, trans_mat, n):
    beta = (trans_mat * E) @ beta_next
    return beta / n


@njit
def bwd_algorithm_single_obs(emission, trans_mat, n):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n])
    """

    n_steps, n_states = emission.shape

    beta = np.ones((n_steps, n_states))

    # for i, e in zip(range(n_steps-1, -1, -1), reversed(em)):
    for i in range(n_steps - 1, 0, -1):
        beta[i - 1] = bwd_step(beta[i], emission[i], trans_mat, n[i])
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
def fwd_bwd_algorithm(alpha0, emissions, trans_mat, gamma=None):
    alpha, n = fwd_algorithm(alpha0=alpha0, emissions=emissions, trans_mat=trans_mat)
    beta = bwd_algorithm(emissions=emissions, n=n, trans_mat=trans_mat)
    if gamma is None:
        gamma = [a * b for (a, b) in zip(alpha, beta)]
        return gamma
    for a, b, g in zip(alpha, beta, gamma):
        g[:] = a * b
    return alpha, beta, n


def viterbi(alpha0, trans_mat, emissions):
    return [viterbi_single_obs(alpha0, trans_mat, e) for e in emissions]


def viterbi_single_obs(alpha0, trans_mat, emissions):
    n_steps, n_states = emissions.shape

    ll = np.ones_like(emissions)
    backtrack = np.zeros_like(emissions, np.int8)

    log_e = np.log(emissions)
    log_t = np.log(trans_mat)

    for i in range(n_steps):
        if i == 0:
            aux_mat = np.log(alpha0) + log_t + log_e[0]
        else:
            aux_mat = ll[i - 1] + log_t + log_e[i]
        ll[i] = np.max(aux_mat, 1)
        backtrack[i] = np.argmax(aux_mat, 1)

    path = np.empty(n_steps, np.int8)
    cursor = np.argmax(ll[-1])
    for i in range(n_steps - 1, -1, -1):
        cursor = backtrack[i, cursor]
        path[i] = cursor

    return path


def update_transitions(
    old_trans_mat, alpha, beta, gamma, emissions, n, est_inbreeding=False
):
    # if no diploid / haploid data, do not update
    breakpoint()
    if len(alpha) == 0:
        return old_trans_mat

    new_trans_mat = np.zeros_like(old_trans_mat)
    n_states = old_trans_mat.shape[0]

    for i in range(n_states):
        for j in range(n_states):
            for a, b, e, n_ in zip(alpha, beta, emissions, n):
                new_trans_mat[i, j] += np.sum(
                    a[:-1, i] * old_trans_mat[i, j] * b[1:, j] * e[1:, j] / n_[1:]
                )

    gamma_sum = np.sum([np.sum(g[:-1], 0) for g in gamma], 0)
    new_trans_mat /= gamma_sum[:, np.newaxis]
    assert np.allclose(np.sum(new_trans_mat, 1), 1)

    # deal with underflow due to absence of state
    if not np.allclose(np.sum(new_trans_mat, 1), 1):
        for i in range(n_states):
            if np.any(np.isnan(new_trans_mat[i])):
                new_trans_mat[i] = 0.0
                new_trans_mat[i, i] - 1.0

    return new_trans_mat
