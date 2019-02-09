from collections import namedtuple, defaultdict
import numpy as np
import pandas as pd

Probs = namedtuple("Probs", ("O", "N", "P_cont", "alpha", "beta", "lib"))


class _IX:
    def __init__(self):
        pass

def data2probs(data, ref, state_ids, cont_id, state_priors=(.5, .5), cont_prior=(1,1)):
    alpha_ix = ["%s_alt" % s for s in state_ids]
    beta_ix = ["%s_ref" % s for s in state_ids]
    cont = "%s_alt" % cont_id, "%s_ref" % cont_id
    pa, pb = cont_prior

    print(alpha_ix, beta_ix)

    P = Probs(
        O = np.array(data.talt),
        N = np.array(data.tref + data.talt),
        P_cont = np.array( (data[cont[0]]+pa) / (data[cont[0]] + data[cont[1]]+ pa + pb)),
        alpha = np.array(ref[alpha_ix]) + state_priors[0],
        beta = np.array(ref[beta_ix]) + state_priors[1],
        lib = np.array(data.lib)
    )
    return P

def bins_from_bed(bed, data, bin_size, pos_mode=False):
    """create a bunch of auxillary data frames for binning

    - bins: columns are chrom_id, chrom, bin_pos, bin_id, map
    - data_bin: chrom, chrom_id, bin_id, snp_id, local_bin_id
    """
    IX = _IX()
    libs =np.unique(data.lib)
    if pos_mode:
        data.map = data.pos
        bed.map = bed.pos
    chroms = np.unique(bed.chrom)
    snp_id = data[['chrom', 'pos', 'map']].drop_duplicates()
    n_snps = snp_id.shape[0]
    snp_id['snp_id'] = range(n_snps)

    IX.SNP2BIN = np.empty((n_snps,2), int)

    data = data.merge(snp_id)
    IX.OBS2SNP = np.array(data['snp_id'])

    bin_loc, data_loc = [], []
    bin0 = 0

    dtype_bin=  dtype=np.dtype([('chrom', 'U2'), ('map', float), ('pos', int), ('id',int),
                                ('chrom_id', int)])
    dtype_data=  dtype=np.dtype([('chrom', 'U2'), ('bin_id', int), ('snp_id',int),
                                 ('loc_id', int), ('chrom_id', int)])

    IX.bin_sizes = []

    for i, chrom in enumerate(chroms):
        map_ = bed.map[bed.chrom == chrom]
        pos = bed.pos[bed.chrom == chrom]
        map_data = data.map[data.chrom == chrom]

        chrom_start = float(np.floor(map_.head(1) / bin_size) * bin_size)
        chrom_end = float(np.ceil(map_.tail(1) / bin_size) * bin_size)

        bins = np.arange(chrom_start, chrom_end, bin_size)
        IX.bin_sizes.append( len(bins))
        bin_ids = range(bin0, bin0 + len(bins))
        _bin = np.empty_like(bins, dtype_bin)
        _bin['chrom'] = chrom
        _bin['pos'] = np.interp(bins, map_, pos)
        _bin['id'] = bin_ids
        _bin['map'] = bins
        _bin['chrom_id'] = i
        bin_loc.append(_bin)

        snp_ids = snp_id.snp_id[snp_id.chrom == chrom]
        dig_data = np.digitize(map_data, bins, right=False) - 1

        dig_snp = np.digitize(snp_id[snp_id.chrom==chrom].map, bins, right=False) -1
        IX.SNP2BIN[snp_ids, 0] = i
        IX.SNP2BIN[snp_ids, 1] = dig_snp

    IX.RG2OBS = dict((l, np.where(data.lib==l)[0]) for l in libs)
    IX.RG2SNP = dict((k, IX.OBS2SNP[v]) for k, v in IX.RG2OBS.items())
    IX.RG2BIN = dict((k, IX.SNP2BIN[v]) for k, v in IX.RG2SNP.items()) 
    IX.OBS2RG = np.array(data.lib)
    IX.OBS2BIN = IX.SNP2BIN[IX.OBS2SNP]

    IX.n_chroms = len(chroms)
    IX.n_bins = len(bins)
    IX.n_snps = len(IX.SNP2BIN)
    IX.n_obs = len(IX.OBS2SNP)
    # snps = np.hstack((bed, np.vstack(snp_id)))

    bins = np.hstack(bin_loc)

    return bins, IX #, data_bin

def init_pars(n_states):
    alpha_0 = np.array([1 / n_states] * n_states)
    trans_mat = np.zeros((n_states, n_states)) + 2e-2
    np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
    cont = defaultdict(lambda: 1e-2)
    return alpha_0, trans_mat, cont


if False:
        _data = np.empty_like(dig_data, dtype_data)
        _data['chrom'] = chrom
        _data['chrom_id'] = i
        _data['bin_id'] = dig_data + bin0
        _data['snp_id'] = snp_ids
        _data['loc_id'] = dig_data
        data_loc.append(_data)

        bin0 += len(bins)
        data_bin = np.hstack(data_loc)

