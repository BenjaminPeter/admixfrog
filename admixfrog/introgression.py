import numpy as np
import pandas as pd
import sys
from collections import namedtuple, defaultdict
from scipy.stats import binom
from scipy.optimize import minimize
from .hmm import get_emissions_cy, split_freqs, update_contamination_cy
from .fwd_bwd import viterbi, fwd_bwd_algorithm
from .baum_welch import update_transitions, get_emissions_py, baum_welch

np.set_printoptions(suppress=True, precision=4)
pbinom = binom.pmf
Freqs = namedtuple("Freqs", ("O", "N", "P_cont", "P", "lib"))

def get_emissions(*args, **kwargs):
    return get_emissions_cy(*args, **kwargs)



def data2freqs(data, state_ids, cont_id, do_hets=True):
    P_homo = [data[s] for s in state_ids]
    if do_hets:
        P_het = []
        for i, s in enumerate(state_ids):
            for s2 in state_ids[i + 1 :]:
                P_het.append((data[s] + data[s2]) / 2)
        P = np.vstack((P_homo, np.vstack(P_het))).T
    else:
        P = np.vstack((P_homo)).T

    f = Freqs(
        O=np.array(data.alt),
        N=np.array(data.ref + data.alt),
        P_cont=np.array(data[cont_id]),
        lib=np.array(data.lib),
        P=P,
    )

    return f

def bins_from_bed(bed, data, bin_size):
    chroms = np.unique(bed.chrom)
    bin_loc, data_loc = [], []
    bin0 = 0

    for i, chrom in enumerate(chroms):
        pos = bed.pos[bed.chrom == chrom]
        pos_data = data.pos[data.chrom == chrom]

        min_ = int(np.floor(pos.head(1) / bin_size) * bin_size)
        max_ = int(np.ceil(pos.tail(1) / bin_size) * bin_size)

        bins = np.arange(min_, max_, bin_size, dtype="int")
        # dig = np.digitize(pos, bins, right=False) - 1
        dig_data = np.digitize(pos_data, bins, right=False) - 1

        bin_ids = range(bin0, bin0 + len(bins))

        bin_loc.append(np.vstack((np.zeros_like(bins) + i, bins, bin_ids)).T)
        # snp_id.append(np.vstack((np.zeros_like(dig)+i, dig)).T)
        data_loc.append(
            np.vstack((np.zeros_like(dig_data) + i, dig_data + bin0, dig_data)).T
        )

        bin0 += len(bins)

    bins = np.vstack(bin_loc)
    # snps = np.hstack((bed, np.vstack(snp_id)))
    data_bin = np.vstack(data_loc)
    return bins, data_bin

def run_hmm(
    infile,
    bedfile,
    split_lib=True,
    state_ids=("AFR", "NEA", "DEN"),
    cont_id="AFR",
    bin_size=1e4,
    garbage_state=True,
    do_hets=True,
    **kwargs
):
    """ run baum welch to find introgressed tracts
    infile formatting:
        - chrom, map: genomic coordinates
        - tref, talt: number of ref, alt reads
        - f_nea, f_afr: frequency of african and neandertal
    """
    F_MAX = 0.51
    data = pd.read_csv(infile)
    data.chrom = data.chrom.astype(str)
    data = data[data.chrom != "X"]
    data = data[data.chrom != "Y"]
    print(data.shape, file=sys.stderr)
    data = data[np.logical_or(data[cont_id] <= F_MAX, data[cont_id] >= 1. - F_MAX)]
    print(data.shape, file=sys.stderr)
    print(F_MAX, file=sys.stderr)
    data.chrom = data.chrom.astype(int)

    data = data[data.ref + data.alt > 0]
    states = list(set(list(state_ids) + [cont_id]))
    cols = ["chrom", "pos", "map", "lib", "ref", "alt"] + states
    data = data[cols]
    data = data.dropna()
    q = np.quantile(data.ref + data.alt, .99)
    data = data[data.ref + data.alt <= q]

    if "lib" not in data or (not split_lib):
        data = data.groupby(
            ("chrom", "pos", "NEA", "DEN", "AFR", "PAN"), as_index=False
        ).agg({"ref": sum, "alt": sum})
        q = np.quantile(data.ref + data.alt, .999)
        data = data[data.ref + data.alt <= q]
        data["lib"] = "lib0"

    # load bed
    bed = pd.read_table(bedfile, header=None)[[0, 2]]
    bed.columns = ["chrom", "pos"]
    try:
        bed = bed[bed.chrom != "X"]
        bed = bed[bed.chrom != "Y"]
        bed.chrom = bed.chrom.astype(int)
    except TypeError:
        pass

    n_states = len(state_ids)
    if do_hets:
        n_states += int(n_states * (n_states - 1) / 2)
    n_states += garbage_state

    alpha_0 = np.array([1 / n_states] * n_states)
    trans_mat = np.zeros((n_states, n_states)) + 2e-2
    np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
    cont = defaultdict(lambda: 1e-2)

    bins, bin_data = bins_from_bed(bed, data, bin_size)
    freqs = data2freqs(data, state_ids, cont_id, do_hets=do_hets)

    # get state names
    n_homo = [s for s in state_ids]
    if do_hets:
        n_het = []
        for i, s in enumerate(state_ids):
            for s2 in state_ids[i + 1 :]:
                n_het.append(s + s2)
        gamma_names = n_homo + n_het
    else:
        gamma_names = n_homo

    if garbage_state:
        gamma_names += ["GARBAGE"]

    gamma, ll, trans_mat, alpha_0, cont, emissions = baum_welch(
        alpha_0=alpha_0,
        trans_mat=trans_mat,
        bins=bins,
        bin_data=bin_data,
        freqs=freqs,
        garbage_state=garbage_state,
        cont=cont,
        gamma_names=gamma_names,
        **kwargs
    )

    viterbi_path = viterbi(alpha_0, trans_mat, emissions)
    viterbi_df = pd.Series(np.hstack(viterbi_path), name="viterbi")

    gamma_out = np.hstack((bins, np.vstack([g[1:] for g in gamma])))
    gamma_out = pd.DataFrame(
        gamma_out, columns=["chrom_id", "bin_pos", "bin_id"] + gamma_names
    )
    gamma_out = pd.concat(
        (gamma_out.reset_index(drop=True), viterbi_df.reset_index(drop=True)), axis=1
    )
    data_out = pd.concat(
        (
            data.reset_index(drop=True),
            gamma_out.iloc[bin_data[:, 1]].reset_index(drop=True),
        ),
        axis=1,
    )

    return gamma_out, data_out, emissions, (ll, trans_mat, cont, alpha_0)
