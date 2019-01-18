from hmm import get_emissions_cy, update_contamination_cy

def test_cy_emissions():
    infile, bedfile="../tmp/Tianyuan_archaicadmixture_nofilter_C0.csv.gz",
    "archaicadmixture.bed"

    e, cont_id, state_ids, bin_size = 1e-2, "AFR", ("NEA", "AFR", "DEN"), 1e4

    n_states = len(state_ids)                                          
    n_states = int(n_states + n_states * (n_states - 1) / 2)           
    alpha_0 = np.array([1/n_states] * n_states)                        
    trans_mat = np.zeros((n_states, n_states)) + 1e-2                  
    np.fill_diagonal(trans_mat, 1 - (n_states-1) * 1e-2)               
    cont = defaultdict(lambda:1e-2)                                    
                                                                       
    bins, bin_data = bins_from_bed(bed, data, bin_size)                
    freqs = data2freqs(data, state_ids, cont_id)                       
                                                                       
    n_homo = [s for s in state_ids]                                    
    n_het = []                                                         
    for i, s in enumerate(state_ids):                                  
        for s2 in state_ids[i+1:]:                                     
            n_het.append(s + s2)                                       
    gamma_names = n_homo + n_het

    e1 = get_emissions_cy(cont, bins, bin_data, freqs, e)

    e2 = get_emissions_py(cont, bins, bin_data, freqs, e)
    all([np.allclose(ee1, ee2) for ee1, ee2 in zip(e1, e2)])

def test_cy_cont():
    infile, bedfile="../tmp/Tianyuan_archaicadmixture_nofilter_C0.csv.gz",
    "archaicadmixture.bed"

    e, cont_id, state_ids, bin_size = 1e-2, "AFR", ("NEA", "AFR", "DEN"), 1e4

    n_states = len(state_ids)                                          
    n_states = int(n_states + n_states * (n_states - 1) / 2)           
    alpha_0 = np.array([1/n_states] * n_states)                        
    trans_mat = np.zeros((n_states, n_states)) + 1e-2                  
    np.fill_diagonal(trans_mat, 1 - (n_states-1) * 1e-2)               
    cont = defaultdict(lambda:1e-2)                                    
                                                                       
    bins, bin_data = bins_from_bed(bed, data, bin_size)                
    freqs = data2freqs(data, state_ids, cont_id)                       
                                                                       
    n_homo = [s for s in state_ids]                                    
    n_het = []                                                         
    for i, s in enumerate(state_ids):                                  
        for s2 in state_ids[i+1:]:                                     
            n_het.append(s + s2)                                       
    gamma_names = n_homo + n_het

    e1 = get_emissions_cy(cont, bins, bin_data, freqs, e)

    alpha, beta, gamma, n =fwd_bwd_algorithm(alpha_0, e1, trans_mat)

infile="../tmp/Tianyuan_archaicadmixture_nofilter_C0.csv.gz" 
bedfile= "archaicadmixture.bed"

data=pd.read_csv(infile)                          
data=data[data.ref+data.alt>0]                    
data=data.dropna()                                
q = np.quantile(data.ref + data.alt, .999)        
data=data[data.ref+data.alt <=q]                  
bed = pd.read_table(bedfile, header=None)[[0, 2]] 
bed.columns = ["chrom", "pos"]


e, cont_id, state_ids, bin_size = 1e-2, "AFR", ("NEA", "AFR", "DEN"), 1e4

n_states = len(state_ids)                                          
n_states = int(n_states + n_states * (n_states - 1) / 2)           
alpha_0 = np.array([1/n_states] * n_states)                        
trans_mat = np.zeros((n_states, n_states)) + 1e-2                  
np.fill_diagonal(trans_mat, 1 - (n_states-1) * 1e-2)               
cont = defaultdict(lambda:1e-2)                                    
                                                                   
bins, bin_data = bins_from_bed(bed, data, bin_size)                
freqs = data2freqs(data, state_ids, cont_id)                       
                                                                   
n_homo = [s for s in state_ids]                                    
n_het = []                                                         
for i, s in enumerate(state_ids):                                  
    for s2 in state_ids[i+1:]:                                     
        n_het.append(s + s2)                                       
gamma_names = n_homo + n_het

e1 = get_emissions_cy(cont, bins, bin_data, freqs, e)

alpha, beta, gamma, n =fwd_bwd_algorithm(alpha_0, e1, trans_mat)
