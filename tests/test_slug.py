from admixfrog.slug.classes import *
from admixfrog.slug.emissions import *
from admixfrog.slug.em import *
import numpy as np
import pytest 

def test_slug_p_gt_diploid():
    tau = np.array([0.4, 0.2, 0, 1])
    F = np.array([0, 0.1, 1, 1])
    res = np.empty((4, 3))

    SNP2SFS = np.arange(4, dtype=int)

    res = p_gt_diploid(tau, F, SNP2SFS=SNP2SFS)

    pred0 = F * ( 1- tau) + (1-F) * (1-tau)**2
    pred1 = 2 * (1-F) * tau * (1-tau)
    pred2 = F * tau + (1-F) * tau ** 2

    pred = np.vstack((pred0, pred1, pred2)).T

    assert np.allclose(res - pred, 0)


def test_slug_p_gt_haploid():
    tau0 = np.arange(10) / 10.0
    res = np.empty((10, 3))

    res = p_gt_haploid(tau=tau0, SNP2SFS = np.arange(10, dtype=int))
    assert np.allclose(res[:, 1], 0)
    assert np.allclose(res[:, 2], tau0)
    assert np.allclose(res[:, 0], 1 - tau0)


def test_slug_fwd_p_x():
    """test example to calculate P(X|.)
    - 3 SNPs, 
    - SNP3 is flipped (alt = anc)
    - psi is 100% der for all SNPs
    - 4 Reads (2, 1, 1) with alleles (ref, alt, alt, ref)
    - 2 Read groups (0, 1, 0, 1)
    - contamination probs are (0.3, 0)


    """

    pg = np.array([[0, .5, .5], [1, 0, 0], [1/3, 1/3, 1/3]])
    pa = np.array([1, 1, 1])
    pc = np.array([0.3, 0])

    READS = np.array((0, 1, 1, 0))
    FLIPPED_READS = np.array((0, 0, 0, 1), 'bool')
    
    n_reads = 4
    n_rgs = 2
    n_snps = 3
    READ2RG = np.array([0, 1, 0, 1])
    READ2SNP = np.array([0, 0, 1, 2])

    bwd_x = bwd_p_o_given_x(READS, FLIPPED_READS, 0, 0)
    #Pr(O|X, e=b=0); 1 if X is the same as READS, except for the flipped SNP
    assert np.allclose(bwd_x[:,1] , (0,1,1,1))

    bwd_g1 = bwd_p_one_o_given_g(bwd_x, pa, pc, READ2SNP, READ2RG, n_reads)
    """Pr(O|G)
        - Read1: we observe ref=anc allele, 30% chance cont; Pr(O|G=0) =  .3 *0 + .7 * 1 = 0.7
        - Read2: we observe alt=der allele, 0% chance cont; Pr(O|G=0) =  .0 *0 + 1 * 0 = 0
        - Read3: we observe alt=der allele, 30% chance cont; Pr(O|G=0) =  .3 *1 + .7 * 0 = 0.3
        - Read4: we observe ref=der allele, 0% chance cont; Pr(O|G=0 (anc)) = 0
    """
    assert np.allclose(bwd_g1[:,0] , (0.7, 0, 0.3, 0))

    bwd_g = bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps)
    assert np.allclose(bwd_g[0] , bwd_g1[0] * bwd_g1[1])
    assert np.allclose(bwd_g[1:] , bwd_g1[2:])


    px_cont = fwd_p_x_cont(pa, READ2SNP) #all 1 in this case
    assert np.allclose(px_cont[:, 1], 1)
    px_nocont = fwd_p_x_nocont(pg, READ2SNP)
    res = fwd_p_x(px_cont, px_nocont, pc, READ2RG)

    """fwd probabilities:
     - reads1, 2: 
        - case1: both are endogenous (0.7 * 1), 
                P(X_1=0 X_2=0| .) = pg[0,1] * .5 * .5 = 1/8
                P(X_1=0 X_2=1| .) = pg[0,1] * .5 * .5 + pg[1,1] * 0  = 1/8
                P(X_1=1 X_2=0| .) = pg[0,1] * .5 * .5 + pg[1,1] * 0  = 1/8
                P(X_1=1 X_2=1| .) = pg[0,1] * .5 * .5 + pg[1,1] * 1  = 5/8
        - case2: only 2 endogenous (0.3 * 1)
                P(X_1=0 X_2=0| .) = 0
                P(X_1=0 X_2=1| .) = 0
                P(X_1=1 X_2=0| .) = pg[0,1] * 0.5 = 0.25
                P(X_1=1 X_2=1| .) = pg[0,1] * 0.5 + pg[1,1] * 1  = 3/4
        total:
                P(X_1=0 X_2=0| .) = 1/8 * 7/10 = 7/80
                P(X_1=0 X_2=1| .) = 7/80
                P(X_1=1 X_2=0| .) = 13/80
                P(X_1=1 X_2=1| .) = 53/80 
        ergo
                p(X_1 = 0) = 14/80 = 0.175
                p(X_2 = 0) = 20/80 = 0.25

      - read3: endo is always 0, cont always 1
            P(X_3=0) = 1 * .7 = 0.7
      - read4: p(cont) = 0, pg is symmetric, P(X_4 = 0.5)
    """
    assert np.allclose(res[:,0], [14/80, 20/80, 0.7, 0.5])


def test_slug_fwd_p_x2():
    fwd_g = np.array([[0.5, .5,0]])
    pa = np.array([.3, 1])
    pc = np.array([0.0])

    READS = np.array((0, 0, 0, 0, 1))
    FLIPPED_READS = np.array([False] * 5)
    
    n_reads = len(READS)
    n_rgs = 1
    n_snps = 1
    READ2RG = np.array([0, 0, 0, 0, 0])
    READ2SNP = np.array([0, 0, 0, 0, 0])



    bwd_x = bwd_p_o_given_x(READS, FLIPPED_READS, 0.01, 0.01)
    bwd_g1 = bwd_p_one_o_given_g(bwd_x, pa, pc, READ2SNP,
                                      READ2RG, n_reads)
    bwd_g = bwd_p_all_o_given_g(bwd_g1, READ2SNP, n_snps)

    px_cont = fwd_p_x_cont(pa, READ2SNP)
    px_nocont = fwd_p_x_nocont(fwd_g,  READ2SNP)
    res = fwd_p_x(px_cont, px_nocont, pc, READ2RG)


    assert np.allclose(res[:,0], 3/4)

def test_error_est():
    """simple test dataset for ensuring algorithm is correct

    """
    data = SlugReads(
        READS = [0, 0, 0, 0, 1, 1, 1, 1],
        psi = [1.],
        READ2RG = [0, 0,0,0, 0, 0, 1, 1],
        READ2SNP = np.zeros(8, dtype='i'),
        FLIPPED = np.zeros(1, dtype='bool'),
        SNP2SFS = [0])

    pars = SlugPars(
        cont = [0, 0],
        tau = [0],
        F = [0],
        e = 0.01,
        b = 0.001
    )
    controller = SlugController(update_eb=True,  update_ftau=False, update_cont=False)
    update_pars(pars, data, controller)
    print( f'e : {pars.prev_e} -> {pars.e}')
    print( f'b : {pars.prev_b} -> {pars.b}')
    assert pars.e == 1/2.
    assert pars.b == 1/2.

    update_pars(pars, data, controller)
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
        FLIPPED = [False],
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
    pars = update_pars(pars, reads, controller)
    assert pars.cont[0] == 0
    assert pars.cont[2] == 1
    assert pars.cont[1] == 0.5
    print(f'C = {pars.cont}')
    ll = calc_full_ll(reads, pars)
    assert pars.ll  <= ll
    print( f'll : {pars.ll} -> {ll}')

