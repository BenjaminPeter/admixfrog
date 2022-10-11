import pandas as pd
import numpy as np
import logging
from itertools import accumulate


def get_runs(targetid, penalty=0.5):
    target, id_ = targetid.target, np.array(targetid.id)
    p0 = np.array(np.log(target + penalty))
    frag_score = 0
    frags = []

    while True:
        p = np.array([k for k in accumulate(p0, lambda x, y: max(x + y, 0))])
        pos_max, score_max = np.argmax(p), np.max(p)
        if score_max == 0.0:
            break
        else:
            pass
            # print(score_max)

        zeros = np.where(p[:pos_max] == 0)[0]
        if len(zeros) == 0:
            pos_min = 0
        else:
            pos_min = np.max(zeros) + 1
        if pos_max != pos_min:
            frags.append((id_[pos_min], id_[pos_max], p[pos_max] - p[pos_min]))
            # print("[%s|%s:%s] : %f" % (targetid.chrom.iloc[0], pos_min, pos_max, p[pos_max] - p[pos_min]))
        p0[pos_min : (pos_max + 1)] = 0

    # for i, score in zip(reversed(id_), reversed(p)):
    #    if score > 0 and score > frag_score:
    #        end_pos, frag_score = i, score
    #    if score == 0 and frag_score > 0:
    #        if i != end_pos:
    #            frags.append((i, end_pos, frag_score))
    #        frag_score = 0
    # if frag_score > 0 and i != end_pos:
    #    frags.append((i, end_pos, frag_score))
    return pd.DataFrame(frags, columns=["start", "end", "score"])


def get_rle(data, states, penalty=0.5):
    coords = data[["chrom", "map", "pos", "id"]]

    het_targets = [[s] for s in states.het_names]  # only heterozygous state
    homo_targets = [[s] for s in states.homo_names]  # only homozygous state
    state_targets = []  # all homo and heterozygous states with one sample
    # inbred_targets = states.roh_names  # all inbred states
    name_state = []

    for i, s in enumerate(states.states):
        name_state.append(s)
        l = []
        for (j1, j2), het_name in zip(states.het, states.het_names):
            if f"h{states[i]}" in data.columns:
                l.append(f"h{states[i]}")
            if i == j1 or i == j2:
                l.append(het_name)
        state_targets.append(l)

    targets = het_targets + state_targets + homo_targets  # + inbred_targets
    names = [*states.het_names, *name_state, *states.homo_names]
    types = ["het"] * states.n_het
    types += ["state"] * states.n_hap
    types += ["homo"] * states.n_homo
    # types += ["inbred"] * states.n_roh

    res = []

    for target, type_, name in zip(targets, types, names):
        logging.info("rle for %s", target)

        data["target"] = np.sum(data[target], 1)
        runs = (
            data[["chrom", "target", "id"]]
            .groupby("chrom")
            .apply(get_runs, penalty=penalty)
            .reset_index()
        )
        del data["target"]
        if "level_1" in runs:
            del runs["level_1"]
        runs["target"] = name
        runs["type"] = type_

        res.append(runs)

    res = pd.concat(res)
    res.score = res.score.astype(float)
    res.start = res.start.astype(int)
    res.end = res.end.astype(int)

    res = res.merge(
        coords, left_on=["chrom", "start"], right_on=["chrom", "id"], how="left"
    )
    res = res.merge(
        coords,
        left_on=["chrom", "end"],
        right_on=["chrom", "id"],
        how="left",
        suffixes=("", "_end"),
    )

    res["len"] = res.end - res.start
    res["map_len"] = res.map_end - res.map
    res["pos_len"] = res.pos_end - res.pos
    res["nscore"] = res.score / res.len
    return res
