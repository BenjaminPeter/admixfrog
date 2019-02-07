from collections import namedtuple, defaultdict
import numpy as np
import pandas as pd

Probs = namedtuple("Probs", ("O", "N", "P_cont", "alpha", "beta", "lib"))
def data2probs(data, ref, state_ids, cont_id, cont_prior=(1,1)):
    alpha_ix = ["%s_alt" % s for s in state_ids]
    beta_ix = ["%s_ref" % s for s in state_ids]
    cont = "%s_alt" % cont_id, "%s_ref" % cont_id
    pa, pb = cont_prior

    print(alpha_ix, beta_ix)

    P = Probs(
        O = np.array(data.talt),
        N = np.array(data.tref + data.talt),
        P_cont = np.array( (data[cont[0]]+pa) / (data[cont[0]] + data[cont[1]]+ pa + pb)),
        alpha = np.array(ref[alpha_ix]),
        beta = np.array(ref[beta_ix]),
        lib = np.array(data.lib)
    )
    return P

def bins_from_bed(bed, data, bin_size, pos_mode=False):
    """create a bunch of auxillary data frames for binning

    - bins: columns are chrom_id, bin_pos, bin_id
    - data_bin: chrom_id, bin_id, snp_id, local_bin_id
    """
    if pos_mode:
        data.map = data.pos
        bed.map = bed.pos
    chroms = np.unique(bed.chrom)
    snp_pos = data[['chrom', 'pos', 'map']].drop_duplicates()
    snp_pos['snp_id'] = range(snp_pos.shape[0])
    data = data.merge(snp_pos)

    bin_loc, data_loc = [], []
    bin0 = 0

    dtype_bin=  dtype=np.dtype([('chrom', 'U2'), ('pos', float), ('id',int)])
    dtype_data=  dtype=np.dtype([('chrom', 'U2'), ('bin_id', int), ('snp_id',int),
                                 ('loc_id', int)])

    for i, chrom in enumerate(chroms):
        pos = bed.map[bed.chrom == chrom]
        pos_data = data.map[data.chrom == chrom]
        snp_ids = data.snp_id[data.chrom == chrom]

        chrom_start = float(np.floor(pos.head(1) / bin_size) * bin_size)
        chrom_end = float(np.ceil(pos.tail(1) / bin_size) * bin_size)

        bins = np.arange(chrom_start, chrom_end, bin_size)
        bin_ids = range(bin0, bin0 + len(bins))
        _bin = np.empty_like(bins, dtype_bin)
        _bin['chrom'] = chrom
        _bin['pos'] = bins
        _bin['id'] = bin_ids
        bin_loc.append(_bin)

        # dig = np.digitize(pos, bins, right=False) - 1
        dig_data = np.digitize(pos_data, bins, right=False) - 1
        _data = np.empty_like(dig_data, dtype_data)
        _data['chrom'] = chrom
        _data['bin_id'] = dig_data + bin0
        _data['snp_id'] = snp_ids
        _data['loc_id'] = dig_data
        data_loc.append(_data)

        bin0 += len(bins)

    bins = np.hstack(bin_loc)
    # snps = np.hstack((bed, np.vstack(snp_id)))
    data_bin = np.hstack(data_loc)
    return bins, data_bin

def init_pars(n_states):
    alpha_0 = np.array([1 / n_states] * n_states)
    trans_mat = np.zeros((n_states, n_states)) + 2e-2
    np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
    cont = defaultdict(lambda: 1e-2)
    return alpha_0, trans_mat, cont


