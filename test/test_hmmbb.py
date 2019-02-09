import pandas as pd
import numpy as np
from collections import defaultdict
from admixfrog import *
from admixfrog.hmm_updates import *
from admixfrog.utils import *
from admixfrog.hmmbb import *
from admixfrog.fwd_bwd import *
from admixfrog.baum_welch import *
from admixfrog.introgression_bb import *
AUTOSOMES = [str(i+1) for i in range(22)]
infile="hmmbb_samples/Chagyrskaya19_allsites.csv.xz"
ref_file="hmmbb_panels/ref_allsites.csv.xz"
#infile="hmmbb_samples/Denisova3_thirdallele.csv.xz"
#ref_file="hmmbb_panels/ref_thirdallele.csv.xz"
#infile="hmmbb_samples/UstIshim_archaicadmixture.csv.xz"
#ref_file="hmmbb_panels/ref_archaicadmixture.csv.xz"

infile="hmmbb_samples/Denisova2_thirdallele.csv.xz"
ref_file="hmmbb_panels/ref_thirdallele.csv.xz"

infile="hmmbb_samples/sample_test.csv.xz"
ref_file="hmmbb_panels/ref_test.csv.xz"

state_ids =  "X", "Y"
#state_ids = "DEN", "UNIF2"
cont_id = "EUR"
split_lib = True
max_iter =1000
ll_tol = 1e-0
prior_alpha, prior_beta = 1., 1.
tau = np.array([1.] * len(state_ids))
assert len(tau) == len(state_ids)
error = 0.

do_haploid = False
optimize_cont = False
pos_mode = False
has_hap=False


states = list(set(list(state_ids) + [cont_id]))
bin_size = 1000
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
if "X" in states:
    ref['X_ref'] = 1
    ref['X_alt'] = 0
if "Y" in states:
    ref['Y_ref'] = 0
    ref['Y_alt'] = 1
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


bins, IX = bins_from_bed(bed=ref.iloc[:,:5], data=data, bin_size=bin_size, pos_mode=pos_mode)
P = data2probs(data, ref, state_ids, cont_id, (prior_alpha, prior_beta))
#assert data.shape[0] == P.O.shape[0] == bin_data.shape[0]
assert ref.shape[0] == P.alpha.shape[0]

# intialize hmm pars
n_states = len(state_ids)
n_states += int(n_states * (n_states - 1) / 2)
n_states += has_hap * len(state_ids)
alpha_0, trans_mat, cont = init_pars(n_states)
cont = defaultdict(lambda : 0.)

emissions = get_emissions_bb_cy(P, IX, cont, tau, error)


#split_ids = split_P(libs, P, bin_data)
print("done loading data")

c_ = np.zeros_like(P.O) + 1e-2

n_seqs = len(emissions)
ll = -np.inf
for it in range(max_iter):
    G = np.zeros( (sum(IX.bin_sizes) + IX.n_chroms, n_states))
    gamma = []
    row0 = 0
    for e in emissions:
        r = e.shape[0] +1
        gamma.append(G[row0:(row0+r)])
        row0 += r
    assert False

    alpha, beta, n = fwd_bwd_algorithm(alpha_0, emissions, trans_mat, gamma)
    ll, old_ll = np.sum([np.sum(np.log(n_i)) for n_i in n]), ll

    print("iter %d [%d/%d]: %s -> %s" % (it, n_seqs, n_states, ll, ll - old_ll))
    if ll - old_ll < ll_tol:
        break

    # update stuff
    alpha_0 = np.linalg.matrix_power(trans_mat, 10000)[0]
    trans_mat = update_transitions(trans_mat, alpha, beta, gamma, emissions, n)




    print(*["%.3f" % a for a in alpha_0], sep="\t")


    n_homo = len(tau)
    n_het = int(n_homo * (n_homo - 1) / 2)

    pg = post_geno(P=P,                    
                              cont=cont,              
                              tau=tau,                
                              IX = IX,                
                              error=error)            
    assert np.all(pg >=0)                  
    assert np.all(pg <=1)                  
    assert np.allclose(np.sum(pg, 2), 1)   


#    update_contamination(cont, error, P, gamma, pg, IX, libs)
    update_contamination2(cont, error, P, gamma, np.array(tau), IX, libs)

    #update_F(tau, gamma, pg[:, :n_homo], P, IX, error, cont)
    #emissions = get_emissions_bb_cy(P, IX, cont, tau, error)
    print(np.sum(emissions[0]))


