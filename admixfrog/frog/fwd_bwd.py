from numba import njit
import numpy as np
import logging


@njit
def fwd_step(alpha_prev, E, trans):
    alpha_new = (alpha_prev @ trans) * E
    n = np.sum(alpha_new)
    return alpha_new / n, n


@njit
def fwd_algorithm_single_obs(alpha0, emission, trans, alpha, n):
    """
    calculate P(X_t | o_[1..t], a0)
    """
    n_steps, n_states = emission.shape
    alpha[0] = alpha0
    n[0] = np.sum(alpha0)
    for i in range(n_steps):
        if i == 0:
            alpha[i], n[i] = fwd_step(alpha0, emission[i], trans)
        else:
            alpha[i], n[i] = fwd_step(alpha[i - 1], emission[i], trans)


# @jit(nopython=True)
def fwd_algorithm(alpha0, emissions, trans, alpha=None, n=None):
    """
    calculate P(X_t | o_[1..t], a0)
    """
    # alpha, n = [], []
    n_seqs = len(emissions)
    if alpha is None:
        alpha = [np.empty_like(emissions[i]) for i in range(n_seqs)]
    if n is None:
        n = [np.empty_like(emissions[i]) for i in range(n_seqs)]

    for i in range(n_seqs):
        fwd_algorithm_single_obs(alpha0, emissions[i], trans, alpha[i], n[i])
    return alpha, n


@njit
def bwd_step(beta_next, E, trans, n):
    beta = (trans * E) @ beta_next
    return beta / n


@njit
def bwd_algorithm_single_obs(emission, trans, n, beta):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n]), return to beta
    """
    n_steps, n_states = emission.shape

    for i in range(n_steps - 1, 0, -1):
        beta[i - 1] = bwd_step(beta[i], emission[i], trans, n[i])
    return beta


# @jit(nopython=True)
def bwd_algorithm(emissions, trans, n, beta=None):
    """
    calculate P(o[t+1..n] | X) / P(o[t+1..n])

    emissions : list[np.array]
        list of emission probabilities, one per observed sequence (i.e. chromosome)
    n : list[np.array]
        list of normalization constants
    trans : np.array<n_states x n_states>
        transition matrix
    """

    n_seqs = len(emissions)
    if beta is None:
        beta = [np.empty((2, 2)) for _ in range(n_seqs)]
    for i in range(n_seqs):
        n_i, em = n[i], emissions[i]
        n_steps, n_states = em.shape
        beta[i] = np.ones((n_steps, n_states))
        bwd_algorithm_single_obs(em, trans, n_i, beta[i])
    return beta


# @jit(nopython=True)
def fwd_bwd_algorithm(alpha0, emissions, trans, gamma=None):
    fwd_algorithm(alpha0=alpha0, emissions=emissions, trans_mat=trans)
    beta = bwd_algorithm(emissions=emissions, n=n, trans=trans)
    if gamma is None:
        gamma = [a * b for (a, b) in zip(alpha, beta)]
        return gamma
    for a, b, g in zip(alpha, beta, gamma):
        g[:] = a * b
    return alpha, beta, n


def fwd_bwd_algorithm2(pars, X):
    """inmemory version of forward-backward algorithm"""
    fwd_algorithm(
        alpha0=pars.alpha0,
        emissions=X.emissions,
        trans=pars.trans,
        alpha=X.alpha,
        n=X.n,
    )
    beta = bwd_algorithm(
        emissions=X.emissions,
        trans=pars.trans,
        beta=X.beta,
        n=X.n,
    )
    for a, b, g in zip(X.alpha, X.beta, X.gamma):
        g[:] = a * b


def viterbi(alpha0, trans, emissions):
    return [viterbi_single_obs(alpha0, trans, e) for e in emissions]


def viterbi_single_obs(alpha0, trans, emissions):
    n_steps, n_states = emissions.shape

    ll = np.ones_like(emissions)
    backtrack = np.zeros_like(emissions, np.uint8)

    log_e = np.log(emissions)
    log_t = np.log(trans)

    for i in range(n_steps):
        if i == 0:
            aux_mat = np.log(alpha0) + log_t + log_e[0]
        else:
            aux_mat = ll[i - 1] + log_t + log_e[i]
        ll[i] = np.max(aux_mat, 1)
        backtrack[i] = np.argmax(aux_mat, 1)

    path = np.empty(n_steps, np.uint8)
    if ll.size == 0:
        return np.zeros(0, dtype=int)
    cursor = np.argmax(ll[-1])
    for i in range(n_steps - 1, -1, -1):
        cursor = backtrack[i, cursor]
        path[i] = cursor

    return path


def update_transitions(pars, X, O):
    old_trans = pars.trans
    # if no diploid / haploid data, do not update
    if len(X.alpha) == 0:
        return pars.trans

    new_trans = np.zeros_like(old_trans)
    n_states = old_trans.shape[0]

    for i in range(n_states):
        for j in range(n_states):
            for a, b, e, n_ in zip(X.alpha, X.beta, X.emissions, X.n):
                new_trans[i, j] += np.sum(
                    a[:-1, i] * old_trans[i, j] * b[1:, j] * e[1:, j] / n_[1:]
                )

    gamma_sum = np.sum([np.sum(g[:-1], 0) for g in X.gamma], 0)
    new_trans /= gamma_sum[:, np.newaxis]
    assert np.allclose(np.sum(new_trans, 1), 1)

    # deal with underflow due to absence of state
    if not np.allclose(np.sum(new_trans, 1), 1):
        for i in range(n_states):
            if np.any(np.isnan(new_trans[i])):
                new_trans[i] = 0.0
                new_trans[i, i] - 1.0

    return new_trans
