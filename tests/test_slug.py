from admixfrog.slug.classes import *
from admixfrog.slug.emissions_reads import *
from admixfrog.slug.em_reads import *
import numpy as np
import pytest 

def test_slug_p_gt_diploid():
    tau0 = 0.4
    tau = np.zeros(4) + tau0
    F = np.array([0, 0.1, 1, 1])
    res = np.empty((4, 3))

    SNP2SFS = np.arange(4, dtype=int)

    res = p_gt_diploid(tau, F, SNP2SFS=SNP2SFS, res=res)

    pred0 = np.array([(1 - tau0) ** 2, 2 * tau0 * (1 - tau0), tau0 ** 2])
    pred2 = np.array([(1 - tau0), 0, tau0])
    pred1 = F[1] * pred2 + (1 - F[1]) * pred0

    pred = np.vstack((pred0, pred1, pred2, pred2))

    assert np.allclose(res - pred, 0)

    #return res, res - pred

def test_slug_p_gt_diploid_flipped():
    tau0 = 0.4
    tau = np.zeros(4) + tau0
    F = np.array([0, 0.1, 1, 1])
    res = np.empty((4, 3))

    SNP2SFS = np.arange(4, dtype=int)

    p_gt_diploid(tau, F, SNP2SFS=SNP2SFS, res=res)

    pred0 = np.array([(1 - tau0) ** 2, 2 * tau0 * (1 - tau0), tau0 ** 2])
    pred2 = np.array([(1 - tau0), 0, tau0])
    pred1 = F[1] * pred2 + (1 - F[1]) * pred0

    pred = np.vstack((pred0, pred1, pred2, pred2[::-1]))

    assert np.allclose(res - pred, 0)
    #return res, res - pred


def test_slug_p_gt_haploid():
    tau0 = np.arange(10) / 10.0
    res = np.empty((10, 3))

    p_gt_haploid(tau=tau0, SNP2SFS = np.arange(10, dtype=int), res=res,
                 )
    assert np.allclose(res[:, 1], 0)
    assert np.allclose(res[:, 2], tau0)
    assert np.allclose(res[:, 0], 1 - tau0)


def test_slug_fwd_p_x():
    pg = np.array([[0, .5, .5], [1, 0, 0]])
    pa = np.array([.3, 1])
    pc = np.array([0.3, 1])

    READS = np.array((0, 0, 1))
    
    n_reads = 3
    n_rgs = 2
    n_snps = 2
    READ2RG = np.array([0, 0, 1])
    READ2SNP = np.array([0, 1, 1])

    bwd_x = bwd_p_o_given_x(READS, 0, 0)
    bwd_g1 = bwd_p_one_o_given_g(bwd_x, pa, pc, READ2SNP,
                                      READ2RG, n_reads)
    bwd_g = bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps)

    px_cont = fwd_p_x_cont(pa, READ2SNP)
    px_nocont = message_fwd_p_x_nocont(pg, bwd_g, bwd_g1, READ2SNP)
    res = fwd_p_x(px_cont, px_nocont, pc, READ2RG)

    assert res[0, 1] == .3 * .3 + .75 * .7
    assert res[1, 1] == .3 
    assert res[2, 1] == 1

def test_slug_fwd_p_x2():
    fwd_g = np.array([[0.5, .5,0]])
    pa = np.array([.3, 1])
    pc = np.array([0.0])

    READS = np.array((0, 0, 0, 0, 1))
    
    n_reads = len(READS)
    n_rgs = 1
    n_snps = 1
    READ2RG = np.array([0, 0, 0, 0, 0])
    READ2SNP = np.array([0, 0, 0, 0, 0])


    bwd_x = bwd_p_o_given_x(READS, 0.01, 0.01)
    bwd_g1 = bwd_p_one_o_given_g(bwd_x, pa, pc, READ2SNP,
                                      READ2RG, n_reads)
    bwd_g = bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps)

    px_cont = fwd_p_x_cont(pa, READ2SNP)
    px_nocont = message_fwd_p_x_nocont(fwd_g, bwd_g, bwd_g1, READ2SNP)
    res = fwd_p_x(px_cont, px_nocont, pc, READ2RG)

    assert 0.4 < res[0, 1]  < 0.5
    assert 0.03 < res[4, 1]  < 0.04 

def test_data1():
    """simple test dataset for ensuring algorithm is correct

    a single SNP with one derived allele
    """
    data = SlugData(
        REF = [0, 9, 1],
        ALT = [4, 1, 0],
        psi = [0],
        OBS2RG = [0, 1, 2],
        OBS2SNP = [0, 0, 0],
        SNP2SFS = [0])

    pars0 = SlugPars(
        cont = [0, 0.5, 1],
        tau = [0],
        F = [0],
        e = 0,
        b= 0
    )


    #return pars0, data

def test_error_est():
    """simple test dataset for ensuring algorithm is correct

    """
    data = SlugReads(
        READS = [0, 0, 0, 0, 1, 1, 1, 1],
        psi = [1.],
        READ2RG = [0, 0,0,0, 0, 0, 1, 1],
        READ2SNP = np.zeros(8, dtype='i'),
        SNP2SFS = [0])

    pars = SlugPars(
        cont = [0, 0],
        tau = [0],
        F = [0],
        e = 0.01,
        b = 0.001
    )
    controller = SlugController(update_eb=True,  update_ftau=False, update_cont=False)
    update_pars_reads(pars, data, controller)
    print( f'e : {pars.prev_e} -> {pars.e}')
    print( f'b : {pars.prev_b} -> {pars.b}')
    assert pars.e == 1/2.
    assert pars.b == 1/2.

    update_pars_reads(pars, data, controller)
    assert pars.prev_ll <= pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b
    print( f'll : {pars.prev_ll} -> {pars.ll}')

def test_cont_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    reads = SlugReads(
        READS = [0,0,0,0,0,0,1,1,0,0, 1,1],
        psi = [1],
        READ2RG = [0,0,0,0,0,0,1,1,1,1, 2, 2],
        READ2SNP = [0] * 12,
        SNP2SFS = [0],
    )

    pars = SlugPars(
        cont = [0.8, 0.5, .5],
        tau = [0],
        F = [0],
        e = 0.,
        b = 0.
    )
    controller = SlugController(update_eb=False,  update_ftau=False,
                                update_cont=True)
    pars = update_pars_reads(pars, reads, controller)
    assert pars.cont[0] == 0
    assert pars.cont[2] == 1
    assert pars.cont[1] == 0.5
    print(f'C = {pars.cont}')
    ll = calc_full_ll_reads(reads, pars)
    assert pars.ll  <= ll
    print( f'll : {pars.ll} -> {ll}')

