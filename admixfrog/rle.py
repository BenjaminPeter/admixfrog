import pandas as pd
import numpy as np
from itertools import accumulate


def get_runs(targetid, penalty=0.5):
    target, id_ = targetid.target, np.array(targetid.id)
    p = np.array(np.log(target + penalty))
    p = [k for k in accumulate(p, lambda x, y: max(x + y, 0))]
    frag_score = 0
    frags = []

    for i, score in zip(reversed(id_), reversed(p)):
        if score > 0 and score > frag_score:
            end_pos, frag_score = i, score
        if score == 0 and frag_score > 0:
            if i != end_pos:
                frags.append((i, end_pos, frag_score))
            frag_score = 0
    if frag_score > 0 and i != end_pos:
        frags.append((i, end_pos, frag_score))
    return pd.DataFrame(frags, columns=['start', 'end', 'score'])


def get_rle(data, states, penalty=0.5):
    coords = data[['chrom', 'map', 'pos', 'id']]
    n_states = len(states)

    het_targets = []  # only heterozygous state
    homo_targets = []  # only homozygous state
    state_targets = []  # all homo and heterozygous states with one sample

    for i in range(n_states):
        homo_targets.append([states[i]])
        for j in range(i + 1, n_states):
            het_targets.append([states[i] + states[j]])

    for i in range(n_states):
        l = [states[i]]
        for j in range(n_states):
            if i < j:
                l.append(states[i] + states[j])
            elif i > j:
                l.append(states[j] + states[i])
        state_targets.append(l)

    targets = het_targets + state_targets + homo_targets
    types = ["het"] * len(het_targets)
    types += ["state"] * len(state_targets)
    types += ["homo"] * len(homo_targets)

    res = []

    for target, type_ in zip(targets, types):
        print(target)

        data['target'] = np.sum(data[target], 1)
        runs = data[['chrom', 'target', 'id']].groupby('chrom').apply(get_runs,penalty=penalty).reset_index()
        del data['target']
        if 'level_1' in runs:
            del runs['level_1']
        runs['target'] = target[0]
        runs['type'] = type_

        res.append(runs)

    res = pd.concat(res)
    res.score = res.score.astype(float)
    res.start = res.start.astype(int)
    res.end = res.end.astype(int)

    res = res.merge(coords,
                    left_on=['chrom', 'start'], 
                    right_on=['chrom', 'id'], 
                    how='left')
    res = res.merge(coords,
                    left_on=['chrom', 'end'], 
                    right_on=['chrom', 'id'], 
                    how='left',
                    suffixes=('', '_end'))


    res['len'] = res.end - res.start
    res['map_len'] = res.map_end - res.map
    res['pos_len'] = res.pos_end - res.pos
    res['nscore'] = res.score / res.len
    return res
