from admixfrog.slug.slug_classes import *
from admixfrog.slug.slug_emissions_reads import *
from admixfrog.slug.slug_em_reads import *
import numpy as np
import pytest 

def test_slug_p_gt_diploid():
    tau0 = 0.4
    tau = np.array([tau0])
    F = np.array([0, 0.1, 1, 1])
    res = np.empty((4, 3))

    SNP2SFS = np.arange(4, dtype=int)

    slug_p_gt_diploid(tau, F, SNP2SFS=SNP2SFS, res=res)

    pred0 = np.array([(1 - tau0) ** 2, 2 * tau0 * (1 - tau0), tau0 ** 2])
    pred2 = np.array([(1 - tau0), 0, tau0])
    pred1 = F[1] * pred2 + (1 - F[1]) * pred0

    pred = np.vstack((pred0, pred1, pred2, pred2))

    assert np.allclose(res - pred, 0)

    return res, res - pred

def test_slug_p_gt_diploid_flipped():
    tau0 = 0.4
    tau = np.zeros(4) + tau0
    F = np.array([0, 0.1, 1, 1])
    res = np.empty((4, 3))

    SNP2SFS = np.arange(4, dtype=int)
    FLIPPED = np.zeros(4, bool)
    FLIPPED[3] = True

    slug_p_gt_diploid(tau, F, SNP2SFS=SNP2SFS, res=res, FLIPPED=FLIPPED)

    pred0 = np.array([(1 - tau0) ** 2, 2 * tau0 * (1 - tau0), tau0 ** 2])
    pred2 = np.array([(1 - tau0), 0, tau0])
    pred1 = F[1] * pred2 + (1 - F[1]) * pred0

    pred = np.vstack((pred0, pred1, pred2, pred2[::-1]))

    assert np.allclose(res - pred, 0)
    return res, res - pred


def test_slug_p_gt_haploid():
    tau0 = np.arange(10) / 10.0
    res = np.empty((10, 3))

    slug_p_gt_haploid(tau0=tau0, SNP2SFS = np.arange(10, dtype=int), res=res)
    assert np.allclose(res[:, 1], 0)
    assert np.allclose(res[:, 2], tau0)
    assert np.allclose(res[:, 0], 1 - tau0)


def test_slug_fwd_p_x():
    pg = np.array([[0, .5, .5], [1, 0, 0]])
    pa = np.array([.3, 1])
    pc = np.array([0.3, 1])

    READS = np.array((0, 0, 1))
    
    class _IX():
        n_reads = 3
        n_rgs = 2
        n_snps = 2
        READ2RG = np.array([0, 0, 1])
        READ2SNP = np.array([0, 1, 1])

    IX = _IX()

    bwd_x = slug_bwd_p_o_given_x(READS, 0, 0)
    bwd_g1 = slug_bwd_p_one_o_given_g(bwd_x, pa, pc, IX.READ2SNP,
                                      IX.READ2RG, IX.n_reads)
    bwd_g = slug_bwd_p_all_o_given_g(bwd_g1, IX.READ2SNP, IX.n_snps)

    px_cont = slug_fwd_p_x_cont(pa, IX.READ2SNP)
    px_nocont = message_fwd_p_x_nocont(pg, bwd_g, bwd_g1, IX.READ2SNP)
    res = slug_fwd_p_x(px_cont, px_nocont, pc, IX.READ2RG)

    assert res[0, 1] == .3 * .3 + .75 * .7
    assert res[1, 1] == .3 
    assert res[2, 1] == 1

def test_slug_fwd_p_x2():
    fwd_g = np.array([[0.5, .5,0]])
    pa = np.array([.3, 1])
    pc = np.array([0.0])

    READS = np.array((0, 0, 0, 0, 1))
    
    class _IX():
        n_reads = len(READS)
        n_rgs = 1
        n_snps = 1
        READ2RG = np.array([0, 0, 0, 0, 0])
        READ2SNP = np.array([0, 0, 0, 0, 0])

    IX = _IX()


    bwd_x = slug_bwd_p_o_given_x(READS, 0.01, 0.01)
    bwd_g1 = slug_bwd_p_one_o_given_g(bwd_x, pa, pc, IX.READ2SNP,
                                      IX.READ2RG, IX.n_reads)
    bwd_g = slug_bwd_p_all_o_given_g(bwd_g1, IX.READ2SNP, IX.n_snps)

    px_cont = slug_fwd_p_x_cont(pa, IX.READ2SNP)
    px_nocont = message_fwd_p_x_nocont(fwd_g, bwd_g, bwd_g1, IX.READ2SNP)
    res = slug_fwd_p_x(px_cont, px_nocont, pc, IX.READ2RG)

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

    pars0 = SlugParsSquare(
        cont = [0, 0.5, 1],
        tau = [0],
        F = [0],
        e = 0,
        b= 0
    )


    return pars0, data

def test_error_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugReads(
        REF = [6, 9],
        ALT = [6, 1],
        psi = [1],
        OBS2RG = [0, 1],
        OBS2SNP = [0, 0],
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
    assert pars.e == 0.5
    assert pars.b == 0.9

    update_pars_reads(pars, data, controller)
    assert pars.prev_ll < pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b
    assert pars.prev_ll  < pars.ll
    print( f'll : {pars.prev_ll} -> {pars.ll}')

def test_cont_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugData(
        REF = [10, 9, 0],
        ALT = [0, 1, 5],
        psi = [1],
        OBS2RG = [0, 1, 2],
        OBS2SNP = [0, 0, 0],
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
    assert pars.cont[1] == 0.1
    print(f'C = {pars.cont}')
    ll = calc_full_ll_g(data, pars)
    assert pars.ll  < ll
    print( f'll : {pars.ll} -> {ll}')


def test_ftau_est_hap():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugData(
        ALT = [0, 1, 5, 1],
        REF = [10, 9, 0, 2],
        psi = [1, 1, .5, 0],
        OBS2RG = [0, 0, 0, 1],
        OBS2SNP = [0, 1, 2, 3],
        SNP2SFS = [0, 0, 1, 1],
        haploid_snps = [2, 3]
    )

    pars = SlugPars(
        cont = [0.0, .5],
        tau = [0.4, .1],
        F = [0.5, .5],
        e = 0.00,
        b = 0.00
    )
    controller = SlugController(update_eb=False,  update_ftau=True, update_cont=False)
    update_pars_reads(pars, data, controller)
    print(f'eb= {pars.e}, {pars.b}')
    print(f'C = {pars.cont}')
    print(f'F = {pars.F}')
    print(f'tau = {pars.tau}')
    controller.update_ftau = False
    update_pars_reads(pars, data, controller)
    print( f'll : {pars.prev_ll} -> {pars.ll}')
    assert pars.prev_ll  < pars.ll
    assert .25 < pars.tau[0] < .26
    assert pars.tau[1] == 1.
    assert pars.F[0] == 0
    return pars


def test_ftau_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugData(
        REF = [10, 9, 0, 2],
        ALT = [0, 1, 5, 1],
        psi = [1, 1, .5, .1],
        OBS2RG = [0, 0, 0, 1],
        OBS2SNP = [0, 1, 2, 3],
        SNP2SFS = [0, 0, 1, 1])

    pars = SlugPars(
        cont = [0.0, .1],
        tau = [0.4, .1],
        F = [0.5, .5],
        e = 0.00,
        b = 0.00
    )
    controller = SlugController(update_eb=False,  update_ftau=True, update_cont=False)
    update_pars_reads(pars, data, controller)
    print(f'eb= {pars.e}, {pars.b}')
    print(f'C = {pars.cont}')
    print(f'F = {pars.F}')
    print(f'tau = {pars.tau}')
    ll = calc_full_ll_g(data, pars)
    print( f'll : {pars.ll} -> {ll}')
    assert .25 < pars.tau[0] < .26
    assert .64 < pars.tau[1] < .65
    assert pars.F[0] == 0
    assert pars.ll <= ll
    breakpoint()
    assert .3 < pars.F[1] < .4
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
        OBS2RG = [0, 0, 0, 0],
        OBS2SNP = [0, 0, 0, 0],
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
    n_obs = 1_000_000
    n_snps = 100_000
    n_sfs = 50
    n_rgs = 40
    avg_cov = 2
    e = 0.001
    b = 0.001

    #psi = 1. / np.random.sample(n_snps)
    #psi /= np.max(psi)
    psi = np.random.sample(n_snps)
    #psi /= np.max(psi)

    true_tau = np.arange(0, 1, step = 1 / n_sfs)
    true_cont = np.arange(0, 1, step = 1 / n_rgs)

    SNP2SFS = np.random.randint(0, n_sfs, size=n_snps)

    G = np.random.binomial(2, true_tau[SNP2SFS], size = n_snps)
    A = np.random.binomial(1, psi, size = n_snps)

    #random assignments. In principle, a SNP might have repeat cats
    OBS2RG = np.random.randint(0, n_rgs, size = n_obs)
    reps = np.random.multinomial(n_obs-n_snps, 
                                 pvals=np.zeros(n_snps)+ 1. / n_snps) + 1
    OBS2SNP = np.repeat(np.arange(n_snps, dtype=np.int), reps)

    p = true_cont[OBS2RG] * psi[OBS2SNP] + (1-true_cont[OBS2RG]) * G[OBS2SNP] / 2. 
    p = e * (1.-p) + p * (1.-b)

    N = np.random.poisson(avg_cov, size = n_obs) + 1
    O = np.random.binomial(N, p, n_obs)


    data = SlugData(
        REF = N-O,
        ALT = O,
        psi = psi,
        OBS2RG = OBS2RG,
        OBS2SNP = OBS2SNP,
        SNP2SFS = SNP2SFS)

    pars = SlugPars(
        cont = np.repeat(0.5, n_rgs),
        tau = true_tau, #np.repeat(.4, n_sfs),
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
    controller = SlugController(update_eb=False, 
                                update_ftau=True, 
                                update_cont=True,
                                n_iter=1000,
                                ll_tol = 1e-4)
    em(pars, data, controller)
    return pars


