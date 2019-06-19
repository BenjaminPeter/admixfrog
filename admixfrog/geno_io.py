"""
read Reich geno format
"""
import numpy as np
from math import ceil
import itertools
import pandas as pd

def row_length(n_ind):
    return max(ceil(n_ind / 4 ), 48)

def read_geno(fname, pops=None, target_ind=None, guess_ploidy=True):
    with open(f'{fname}.geno', "rb") as f:
        head = f.read(48).strip(b"\x00").split()
        print(head)
        assert head[0] == b"GENO"
        n_ind, n_snp = int(head[1]), int(head[2])
        hash_ind, hash_snp = head[3:]
        rlen = row_length(n_ind)

    z = np.fromfile(f'{fname}.geno', dtype=np.uint8).reshape(n_snp+1, rlen)[1:]
    X = np.unpackbits(z, axis=1).reshape(n_snp, rlen * 4, 2)
    del z


    ind = pd.read_csv(f'{fname}.ind', sep=" ", names=['id', 'sex', 'pop'],
                      skipinitialspace=True)

    snp = pd.read_csv(f'{fname}.snp', sep=" ", skipinitialspace=True,
        names=["snp", "chrom", "map", 'pos', 'ref', 'alt'])

    #chromosome formatting
    snp.chrom = snp.chrom.astype(str)
    snp.loc[snp.chrom=='23', 'chrom'] = 'X'  
    snp.loc[snp.chrom=='24', 'chrom'] = 'Y'  
    snp.loc[snp.chrom=='90', 'chrom'] = 'mt' 
    chrom_order = pd.unique(snp.chrom)
    snp.chrom = snp.chrom.astype('category')                   
    snp.chrom.cat.reorder_categories(chrom_order, inplace=True)


    ix = pd.MultiIndex.from_frame(ind[['pop', 'sex', 'id']]) 

    Y = pd.DataFrame(X[:, :n_ind, 0] + 2 * X[:, :n_ind, 1], columns = ix)
    Y.index = pd.MultiIndex.from_frame(snp[['chrom', 'pos', 'map', 'ref', 'alt']])
    del X


    if target_ind is not None:
        sex = Y.xs(target_ind, level='id', axis=1).columns.get_level_values('sex')[0]
        Y['TARGET', sex, target_ind] = Y.xs(target_ind, level='id', axis=1)
        if pops is not None:
            pops.append("TARGET")

    if pops is not None:
        Y = Y[pops]

    if guess_ploidy:
        ploidy = Y.agg(lambda x: np.max(x * (x<3)))
        Y = Y.T.set_index(pd.Index(ploidy, name='ploidy'), append=True).T
    else:
        Y = Y.T.set_index(pd.Index([2] * Y.shape[1], name='ploidy'), append=True).T

    return Y

def ref_alt(Y):
    if 'Y' in Y.index:
        Y.loc['Y'] = Y.loc['Y'].transform(lambda x: np.where(x>1, 3, x)).values
    if 'mt' in Y.index:
        Y.loc['mt'] = Y.loc['mt'].transform(lambda x: np.where(x>1, 3, x)).values

    ref = Y.groupby(lambda x:x, level=0, axis=1).agg(ref_count)
    ref.rename(columns = lambda x: f'{x}_ref', inplace=True)
    Y.loc['1':'22'] = Y.loc['1':'22'].transform(lambda x: x.name[3] - x).values

    if 'X' in Y.index:
        Y.loc['X'] = Y.loc['X'].transform(lambda x: (x.name[3] if x.name[1] != 'M' else 1) - x).values

    if 'Y' in Y.index: #  always haploid
        Y.loc['Y'] = Y.loc['Y'].transform(lambda x: 1 - x).values

    if 'mt' in Y.index: # always haploid
        Y.loc['mt'] = Y.loc['mt'].transform(lambda x: 1 - x).values

    alt = Y.groupby(lambda x:x, level=0, axis=1).agg(ref_count)
    alt.rename(columns = lambda x: f'{x}_alt', inplace=True)

    df = ref.merge(alt, on=['chrom', 'pos', 'map', 'ref', 'alt'])
    df.rename(columns = {'TARGET_ref' : 'tref', 'TARGET_alt' : 'talt'}, 
              inplace=True)

    return df
    
def ref_alt2(Y):

    # for haploid chroms set all '2' counts to missing
    if 'Y' in Y.index:
        Y.loc['Y'] = Y.loc['Y'].transform(lambda x: np.where(x>1, 3, x)).values
    if 'mt' in Y.index:
        Y.loc['mt'] = Y.loc['mt'].transform(lambda x: np.where(x>1, 3, x)).values

    # get ref allele
    ref = Y.groupby(lambda x:x, level=0, axis=1).agg(ref_count)
    ref.rename(columns = lambda x: f'{x}_ref', inplace=True)

    
    Y.loc['1':'22'] = Y.loc['1':'22'].transform(lambda x: 2 - x).values

    if 'X' in Y.index:
        Y.loc['X'] = Y.loc['X'].transform(lambda x: - x + (x.name[1]!='M') +1).values

    if 'Y' in Y.index:
        Y.loc['Y'] = Y.loc['Y'].transform(lambda x: 1 - x).values

    if 'mt' in Y.index:
        Y.loc['mt'] = Y.loc['mt'].transform(lambda x: 1 - x).values

    alt = Y.groupby(lambda x:x, level=0, axis=1).agg(ref_count)
    alt.rename(columns = lambda x: f'{x}_alt', inplace=True)

    df = ref.merge(alt, on=['chrom', 'pos', 'map', 'ref', 'alt'])
    df.rename(columns = {'TARGET_ref' : 'tref', 'TARGET_alt' : 'talt'}, 
              inplace=True)

    return df


def ref_count(x): 
    v = np.sum(x * (x<3) , axis=1) 
    return v 

def read_geno_ref(*args, **kwargs):
    Y = read_geno(*args, **kwargs)
    return ref_alt(Y)
