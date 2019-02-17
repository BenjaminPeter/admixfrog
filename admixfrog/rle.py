import pandas as pd
import numpy as np


def get_rle(data, states, cutoff=0.8, lmin=3):
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

        data["target"] = np.sum(data[target], 1) > cutoff
        data["block"] = (data.target.shift(1) != data.target).astype(int).cumsum()
        x = data.reset_index().groupby(["target", "block"])["index"].apply(len)
        x = x.reset_index().sort_values('block').reset_index(drop=True)
        x = x.rename({"index": "len"}, axis="columns")
        x['gap'] = (x.len < lmin) & np.logical_not(x.target)
        x = data.merge(x, on=['block', 'target'])
        x.loc[x.gap, 'target'] = True
        print("closing gaps of size %s" % sum(x.gap))

        #second round
        del x['len']
        x["block"] = (x.target.shift(1) != x.target).astype(int).cumsum()
        y = x.reset_index().groupby(["target", "block"])["index"].apply(len)
        y = y.reset_index().sort_values('block').reset_index(drop=True)
        y = y.rename({"index": "len"}, axis="columns")
        x = x.merge(y)

        x["map_end"] = x.groupby("chrom").map.shift(-1)
        x["pos_end"] = x.groupby("chrom").pos.shift(-1)
        del data["target"]
        del data["block"]

        agg = (
            x.loc[x.target]
            .groupby(["chrom", "block"])
            .agg(
                {
                    "map": min,
                    "map_end": max,
                    "pos": min,
                    "pos_end": max,
                    "n_snps": sum,
                    "gap": sum,
                    "len": np.median,
                }
            )
            .reset_index()
        )
        agg["type"] = type_
        agg["target"] = target[0]
        res.append(agg)

    res = pd.concat(res)
    return res
