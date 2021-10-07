from admixfrog.slug.classes import *
from admixfrog.slug.emissions_reads import *
from admixfrog.slug.em_reads import *
import numpy as np
import pytest 

def test_error_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugReads(
        READS = [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1],
        psi = [1],
        READ2RG = [0] * 6 + [1] * 5,
        READ2SNP = [0] * 11,
        FLIPPED = [False],
        SNP2SFS = [0])

    pars = SlugParsSquare(
        cont = [0, 1],
        tau = [0],
        F = [0],
        e = 0.01,
        b = 0.001
    )
    controller = SlugController(update_eb=True,  update_ftau=False, update_cont=False)
    update_pars_reads(pars, data, controller)
    print( f'e : {pars.prev_e} -> {pars.e}')
    print( f'b : {pars.prev_b} -> {pars.b}')
    print( f'll : {pars.prev_ll} -> {pars.ll}')
    assert pars.e == 0.5
    assert pars.b == 0.8

    update_pars_reads(pars, data, controller)
    assert pars.prev_ll == pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b

def test_cont_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugReads(
        READS = [0, 0, 0, 0, 
                 0, 0, 0, 1,
                 1, 1, 1, 1],
        psi = [1],
        READ2RG = [0] * 4 + [1] * 4 + [2] * 4,
        READ2SNP = [0] * 12,
        FLIPPED = [False],
        SNP2SFS = [0])

    pars = SlugPars(
        cont = [0.8, 0.5, .5],
        tau = [0],
        F = [0],
        e = 0.,
        b = 0.
    )
    controller = SlugController(update_eb=False,  update_ftau=False,
                                update_cont=True)
    pars = update_pars_reads(pars, data, controller)
    assert pars.cont[0] == 0
    assert pars.cont[2] == 1
    assert pars.cont[1] == 0.25
    print(f'C = {pars.cont}')
    ll = calc_full_ll_reads(data, pars)
    assert pars.ll  > pars.prev_ll
    print( f'll : {pars.ll} -> {ll}')


def test_ftau_est_hap():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugReads(
        READS = [0, 0, 0, 0, 
                 0, 0, 0, 1,
                 1, 1, 1, 1,
                 0, 0, 1],
        psi = [1, 1, .5, 1],
        READ2RG = [0] * 4 + [0] * 4 + [0] * 4 + [1] * 3,
        READ2SNP = [0] * 4 + [1] * 4 + [2] * 4 + [3] * 3,
        FLIPPED = [False] * 4,
        haploid_snps = [2,3],
        SNP2SFS = [0, 0, 1, 1])

    pars = SlugPars(
        cont = [0.0, .1],
        tau = [0.4, .1],
        F = [0.5, .5],
        e = 0.00,
        b = 0.00
    )
    ll0 = calc_full_ll_reads(data, pars)
    controller = SlugController(update_eb=False,  update_ftau=True, update_cont=False)
    update_pars_reads(pars, data, controller)
    print(f'eb= {pars.e}, {pars.b}')
    print(f'C = {pars.cont}')
    print(f'F = {pars.F}')
    print(f'tau = {pars.tau}')
    print( f'll : {ll0} -> {pars.ll} ')
    assert .25 < pars.tau[0] < .26
    assert .50 == pars.tau[1]
    assert pars.ll > ll0
    #assert .3 < pars.F[1] < .4
    #assert pars.F[0] == 0.5
    return pars

def test_ftau_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugReads(
        READS = [0, 0, 0, 0, 
                 0, 0, 0, 1,
                 1, 1, 1, 1,
                 0, 0, 1],
        psi = [1, 1, .5, 1],
        READ2RG = [0] * 4 + [0] * 4 + [0] * 4 + [1] * 3,
        READ2SNP = [0] * 4 + [1] * 4 + [2] * 4 + [3] * 3,
        FLIPPED = [False] * 4,
        SNP2SFS = [0, 0, 1, 1])

    pars = SlugPars(
        cont = [0.0, .1],
        tau = [0.4, .1],
        F = [0.5, .5],
        e = 0.00,
        b = 0.00
    )
    ll0 = calc_full_ll_reads(data, pars)
    controller = SlugController(update_eb=False,  update_ftau=True, update_cont=False)
    update_pars_reads(pars, data, controller)
    print(f'eb= {pars.e}, {pars.b}')
    print(f'C = {pars.cont}')
    print(f'F = {pars.F}')
    print(f'tau = {pars.tau}')
    print( f'll : {ll0} -> {pars.ll} ')
    assert .25 < pars.tau[0] < .26
    assert .50 < pars.tau[1] < .51
    assert pars.ll > ll0
    #assert .3 < pars.F[1] < .4
    #assert pars.F[0] == 0.5
    return pars

@pytest.mark.skip(reason="NYI")
def test_delta_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugData(
        REF = [10, 9, 0, 2],
        ALT = [0, 1, 5, 1],
        psi = [.5],
        READ2RG = [0, 0, 0, 0],
        READ2SNP = [0, 0, 0, 0],
        SNP2SFS = [0])

    pars = SlugPars(
        cont = [1.0],
        tau = [0.4],
        F = [0.5],
        e = 0.00,
        b = 0.00
    )
    controller = SlugController(update_eb=False,  update_ftau=False, update_cont=False)
    update_pars_reads(pars, data, controller)
    print(f'eb= {pars.e}, {pars.b}')
    print(f'C = {pars.cont}')
    print(f'F = {pars.F}')
    print(f'tau = {pars.tau}')
    update_pars_reads(pars, data, controller)
    print( f'll : {pars.prev_ll} -> {pars.ll}')
    assert .25 < pars.tau[0] < .26
    assert .64 < pars.tau[1] < .65
    assert pars.F[0] == 0
    assert .3 < pars.F[1] < .4
    return pars

@pytest.mark.skip(reason="takes very long")
def test_update_large():
    """large random test dataset for perrrrformance testing

    """
    np.random.seed(1)
    n_reads = 10_000_000
    n_snps = 100_000
    n_sfs = 40
    n_rgs = 20
    e = 0.001
    b = 0.001

    #psi = 1. / np.random.sample(n_snps)
    #psi /= np.max(psi)
    psi = np.random.sample(n_snps)
    #psi /= np.max(psi)

    true_tau = np.arange(0, 1, step = 1 / n_sfs)
    true_cont = np.arange(0, 1, step = 1 / n_rgs)

    SNP2SFS = np.random.randint(0, n_sfs, size=n_snps)
    #FLIPPED = np.zeros(n_snps, bool)


    G = np.random.binomial(2, true_tau[SNP2SFS], size = n_snps)
    A = np.random.binomial(1, psi, size = n_snps)

    #random assignments. In principle, a SNP might have repeat cats
    READ2RG = np.random.randint(0, n_rgs, size = n_reads)
    reps = np.random.multinomial(n_reads-n_snps, 
                                 pvals=np.zeros(n_snps)+ 1. / n_snps) + 1
    READ2SNP = np.repeat(np.arange(n_snps, dtype=np.int), reps)

    p = true_cont[READ2RG] * psi[READ2SNP] + (1-true_cont[READ2RG]) * G[READ2SNP] / 2. 
    p = e * (1.-p) + p * (1.-b)
    O = np.random.binomial(1, p, n_reads)

    FLIPPED = np.zeros(n_snps, bool)
    #FLIPPED = np.array(np.random.binomial(1, .4, size =n_snps), bool)
    #O[FLIPPED[READ2SNP]] = 1 - O[FLIPPED[READ2SNP]]


    np.set_printoptions(suppress=True,precision = 3)

    data = SlugReads(
        READS = O,
        psi = psi,
        READ2RG = READ2RG,
        READ2SNP = READ2SNP,
        FLIPPED = FLIPPED,
        SNP2SFS = SNP2SFS)
    pars = SlugParsSquare(
        cont = np.repeat(0.5, n_rgs),
        tau = np.repeat(.4, n_sfs),
        F = np.repeat(.5, n_sfs),
        e = 0.001,
        b = 0.001
    )
    #update_pars_reads(pars, data, data, 
    #            True, True, True)
    #print(f'eb= {pars.e}, {pars.b}')
    #print(f'C = {pars.cont}')
    #print(f'F = {pars.F}')
    #print(f'tau = {pars.tau}')
    #update_pars_reads(pars, data, data, 
    #            True, True, True)
    #print( f'll : {pars.prev_ll} -> {pars.ll}')
    controller = SlugController(update_eb=True, 
                                update_ftau=True, 
                                update_cont=True,
                                n_iter=1000,
                                ll_tol = 1e-4)
    #em(pars, data, controller)
    squarem(pars, data, controller)
    return pars


