import pandas as pd
import numpy as np
from admixfrog import *
AUTOSOMES = [str(i+1) for i in range(22)]
infile="hmmbb_samples/Oase_archaicadmixture.csv.xz"
ref_file="hmmbb_panels/ref_archaicadmixture.csv.xz"
infile="hmmbb_samples/altai_thirdallele.csv.xz"
ref_file="hmmbb_panels/ref_thirdallele.csv.xz"

state_ids = "AFR", "VIN", "DEN"
state_ids = "DEN", "UNIF2"
cont_id = "AFR"
split_lib = True
max_iter =1000
ll_tol = 1e-1
alpha_prior, beta_prior = .5, .5
fst = (1,1,1, 1)

do_haploid = False
optimize_cont = False
pos_mode = False
has_hap=False


states = list(set(list(state_ids) + [cont_id]))
bin_size = 10000
bin_size = bin_size  if pos_mode else bin_size * 1e-6



dtype_ = dict(chrom=str)


# initialize data
data = pd.read_csv(infile, dtype=dtype_).dropna()
if "lib" not in data or (not split_lib):
    data = data.groupby(
        ["chrom", "pos", "map"], as_index=False
    ).agg({"tref": sum, "talt": sum})
    data["lib"] = "lib0"

# rm sites with extremely high coverage
data = data[data.tref + data.talt > 0]    
q = np.quantile(data.tref + data.talt, .999)
data = data[data.tref + data.talt <= q]    


# prep data
libs = pd.unique(data.lib)
ref = pd.read_csv(ref_file, dtype=dtype_)
if "UNIF" in states:
    ref['UNIF_ref'] = 1
    ref['UNIF_alt'] = 1
if "UNIF2" in states:
    ref['UNIF2_ref'] = 0
    ref['UNIF2_alt'] = 0
if "PAN" in states:
    ref['PAN_ref'] /= 2
    ref['PAN_alt'] /= 2
ix = list(ref.columns[:5])
suffixes=["_alt", "_ref"]
cols = ix + [s + x for s in states for x in suffixes]
ref = ref[cols].dropna()
#ref = ref[ref.chrom != "X"]

if pos_mode:
    data.map = data.pos
    ref.map = ref.pos

ref = ref.merge(data.iloc[:, :3].drop_duplicates()).drop_duplicates()
data= data.merge(ref)
ref = ref.sort_values(['chrom', 'map', 'pos'])
data = data.sort_values(['chrom', 'map', 'pos'])


bins, bin_data = bins_from_bed(bed=ref.iloc[:,:5], data=data, bin_size=bin_size, pos_mode=pos_mode)
P = data2probs(data, ref, state_ids, cont_id)
assert data.shape[0] == P.O.shape[0] == bin_data.shape[0]
assert ref.shape[0] == P.alpha.shape[0]

# intialize hmm pars
n_states = len(state_ids)
n_states += int(n_states * (n_states - 1) / 2)
n_states += has_hap * len(state_ids)
alpha_0, trans_mat, cont = init_pars(n_states)
#cont_ = pd.read_csv("../admixfrog/cont.txt")
#cont = dict(("_".join([i.lib, i.deam]), i.cont) for _, i in cont_.iterrows())
#cont["SR_deam"] = 0.2
#cont["SR_nodeam"] = 0.2

emissions = get_emissions_bb_py(P, bins, bin_data, cont,
                                alpha_prior, beta_prior, fst, has_hap=has_hap)

print("done loading data")

n_seqs = len(emissions)
ll = -np.inf
for it in range(max_iter):
    alpha, beta, gamma, n = fwd_bwd_algorithm(alpha_0, emissions, trans_mat)
    ll, old_ll = np.sum([np.sum(np.log(n_i)) for n_i in n]), ll

    print("iter %d [%d/%d]: %s -> %s" % (it, n_seqs, n_states, ll, ll - old_ll))
    if ll - old_ll < ll_tol:
        break

    # update stuff
    alpha_0 = np.linalg.matrix_power(trans_mat, 10000)[0]
    trans_mat = update_transitions(trans_mat, alpha, beta, gamma, emissions, n)
    print(*["%.3f" % a for a in alpha_0], sep="\t")
