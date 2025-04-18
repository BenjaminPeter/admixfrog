from admixfrog.slug.classes import *
from admixfrog.slug.emissions import *
from admixfrog.slug.em import *
import numpy as np
import pytest


def test_error_bias_est():
    """simple test dataset for ensuring algorithm is correct
    SNP 1 has ref=anc, alt = derived (FLIPPED = False)
    SNP 2 has ref=der, alt = anc (FLIPPED = TRUE)

    - error is the probability a read has the alt allele when it really should
        have reference
    - bias is the prob a read has ref allele when it should have alt
    - both SNPs have truth = homozygous anc, i.e no error
    - SNP1 has 6 reads, half of which are derived
    - SNP2 has 4 reads, one of which is derived

    """
    data = SlugReads(
        READS=[0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
        psi=[1, 1],
        READ2RG=[0] * 10,
        READ2SNP=[0] * 6 + [1] * 4,
        FLIPPED=[False, True],
        SNP2SFS=[0, 0],
    )

    pars = SlugPars(cont=[0.0], tau=[0], F=[0], e=0.50, b=0.55)
    controller = SlugController(
        update_eb=True, update_ftau=False, update_cont=False, update_bias=True
    )
    update_pars(pars, data, controller)
    print(f"e : {pars.prev_e} -> {pars.e}")
    print(f"b : {pars.prev_b} -> {pars.b}")
    print(f"ll : {pars.prev_ll} -> {pars.ll}")
    assert pars.e == 0.5
    assert pars.b == 0.25

    update_pars(pars, data, controller)
    assert pars.prev_ll == pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b


def test_error_est():
    """simple test dataset for ensuring algorithm is correct
    SNP 1 has ref=anc, alt = derived (FLIPPED = False)
    SNP 2 has ref=der, alt = anc (FLIPPED = TRUE)

    - error is the probability a read has the alt allele when it really should
        have reference
     - bias is the prob a read has ref allele when it should have alt

     both SNPs have truth = homozygous anc, i.e no error

    """
    data = SlugReads(
        READS=[0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
        psi=[1, 1],
        READ2RG=[0] * 10,
        READ2SNP=[0] * 6 + [1] * 4,
        FLIPPED=[False, True],
        SNP2SFS=[0, 0],
    )

    pars = SlugPars(cont=[0.0], tau=[0], F=[0], e=0.50, b=0.25)
    controller = SlugController(
        update_eb=True, update_ftau=False, update_cont=False, update_bias=False
    )
    update_pars(pars, data, controller)
    print(f"e : {pars.prev_e} -> {pars.e}")
    print(f"b : {pars.prev_b} -> {pars.b}")
    print(f"ll : {pars.prev_ll} -> {pars.ll}")
    assert pars.e == 0.4
    assert pars.b == pars.e

    update_pars(pars, data, controller)
    assert pars.prev_ll == pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b


def test_error_bias_est2():
    """simple test dataset for ensuring algorithm is correct
    SNP 1 has ref=anc, alt = derived (FLIPPED = False)
         thus if tau=0 (only expect anc) error is the proportion of der
    SNP 2 has ref=der, alt = anc (FLIPPED = TRUE)
         thus if tau=1 (only expect der) error is the proportion of alt)

    - error is the probability a read has the alt allele when it really should
        have reference
     - bias is the prob a read has ref allele when it should have alt

     both SNPs have truth = ref

    """
    data = SlugReads(
        READS=[0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
        psi=[1, 1],
        READ2RG=[0] * 10,
        READ2SNP=[0] * 6 + [1] * 4,
        FLIPPED=[False, True],
        SNP2SFS=[0, 1],
    )

    pars = SlugPars(cont=[0.0], tau=[0, 1], F=[0, 0], e=0.70, b=0.35)
    controller = SlugController(
        update_eb=True, update_ftau=False, update_cont=False, update_bias=True
    )
    update_pars(pars, data, controller)

    print(f"e : {pars.prev_e} -> {pars.e}")
    print(f"b : {pars.prev_b} -> {pars.b}")
    print(f"ll : {pars.prev_ll} -> {pars.ll}")
    assert pars.e == 6 / 10
    assert pars.b == 0

    update_pars(pars, data, controller)
    assert pars.prev_ll == pars.ll
    assert pars.e == pars.prev_e
    assert pars.b == pars.prev_b


def test_cont_est():
    """simple test dataset for ensuring algorithm is correct

    we have 2 SNP (1 flipped), with 3 RGs with 4 reads each;
    tau=0, psi=1, i.e. all contaminant reads are derived and endo are anc
    """
    data = SlugReads(
        READS=[0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        psi=[1, 1],
        READ2RG=[0] * 4 + [1] * 4 + [2] * 4,
        READ2SNP=[0] * 10 + [1] * 2,
        FLIPPED=[False, True],
        SNP2SFS=[0, 0],
    )

    pars = SlugPars(cont=[0.8, 0.5, 0.5], tau=[0, 0], F=[0, 0], e=0.0, b=0.0)
    controller = SlugController(update_eb=False, update_ftau=False, update_cont=True)
    pars = update_pars(pars, data, controller)
    assert pars.cont[0] == 0
    assert pars.cont[1] == 0.25
    assert pars.cont[2] == 0.5
    print(f"C = {pars.cont}")
    ll = calc_full_ll(data, pars)
    assert pars.ll > pars.prev_ll
    print(f"ll : {pars.ll} -> {ll}")


def test_ftau_est_hap():
    """simple test dataset for ensuring algorithm is correct

    2 SFS cats, 2  RG, 4 SNPs, no error;
    RGs have cont = 0, .0

    """
    data = SlugReads(
        READS=[
            0,
            0,
            0,
            0,  # SNP0, RG0, SFS0
            0,
            0,
            0,
            1,  # SNP1, RG0, SFS1
            1,
            1,
            1,
            1,  # SNP2, RG1, haploid, SFS0
            0,
            0,
            0,
        ],  # SNP3< from RG1, haploid, SFS1
        psi=[1, 1, 1, 1],
        READ2RG=[0] * 4 + [0] * 4 + [1] * 4 + [1] * 3,
        READ2SNP=[0] * 4 + [1] * 4 + [2] * 4 + [3] * 3,
        FLIPPED=[False] * 4,
        haploid_snps=[2, 3],
        SNP2SFS=[0, 1, 0, 1],
    )

    pars = SlugPars(cont=[0.0, 0.0], tau=[0.4, 0.1], F=[0.5, 0.5], e=0.00, b=0.00)
    ll0 = calc_full_ll(data, pars)
    controller = SlugController(update_eb=False, update_ftau=True, update_cont=False)
    update_pars(pars, data, controller)
    assert pars.ll > ll0
    update_pars(pars, data, controller)
    assert pars.ll > ll0
    update_pars(pars, data, controller)
    assert pars.ll > ll0

    print(f"eb= {pars.e}, {pars.b}")
    print(f"C = {pars.cont}")
    print(f"F = {pars.F}")
    print(f"tau = {pars.tau}")
    print(f"ll : {ll0} -> {pars.ll} ")
    assert np.allclose(pars.tau, [0.5, 0.25], atol=1e-4)
    assert pars.ll > ll0


def test_ftau_est():
    """simple test dataset for ensuring algorithm is correct"""
    data = SlugReads(
        READS=[
            0,
            0,
            0,
            0,  # SNP1 only has ref, hence is likely homoz
            0,
            0,
            0,
            1,  # SNP2 has both, hence is het
            0,
            0,
            0,
            0,  # SNP3 only has ref, hence is homoz
            1,
            1,
            1,
            1,
        ],
        psi=[1, 1, 0.5, 1],
        READ2RG=[0] * 16,
        READ2SNP=[0] * 4 + [1] * 4 + [2] * 4 + [3] * 4,
        FLIPPED=[False] * 2 + [True] * 2,
        SNP2SFS=[0, 0, 1, 1],
    )

    pars = SlugPars(cont=[0.0], tau=[0.4, 0.1], F=[0.5, 0.5], e=0.00, b=0.00)
    ll0 = calc_full_ll(data, pars)
    controller = SlugController(update_eb=False, update_ftau=True, update_cont=False)
    update_pars(pars, data, controller)
    print(f"eb= {pars.e}, {pars.b}")
    print(f"C = {pars.cont}")
    print(f"F = {pars.F}")
    print(f"tau = {pars.tau}")
    print(f"ll : {ll0} -> {pars.ll} ")
    assert 0.25 < pars.tau[0] < 0.26
    assert 0.45 < pars.tau[1] < 0.5
    assert pars.ll > ll0
    assert pars.F[1] > 0.8
    # assert pars.F[0] == 0.5


@pytest.mark.skip(reason="NYI")
def test_delta_est():
    """simple test dataset for ensuring algorithm is correct

    one RG with only endo (all ref) and one RG with only cont ( all 1)
    """
    data = SlugData(
        REF=[10, 9, 0, 2],
        ALT=[0, 1, 5, 1],
        psi=[0.5],
        READ2RG=[0, 0, 0, 0],
        READ2SNP=[0, 0, 0, 0],
        SNP2SFS=[0],
    )

    pars = SlugPars(cont=[1.0], tau=[0.4], F=[0.5], e=0.00, b=0.00)
    controller = SlugController(update_eb=False, update_ftau=False, update_cont=False)
    update_pars(pars, data, controller)
    print(f"eb= {pars.e}, {pars.b}")
    print(f"C = {pars.cont}")
    print(f"F = {pars.F}")
    print(f"tau = {pars.tau}")
    update_pars(pars, data, controller)
    print(f"ll : {pars.prev_ll} -> {pars.ll}")
    assert 0.25 < pars.tau[0] < 0.26
    assert 0.64 < pars.tau[1] < 0.65
    assert pars.F[0] == 0
    assert 0.3 < pars.F[1] < 0.4


@pytest.mark.skip(reason="takes very long")
def test_update_large():
    """large random test dataset for perrrrformance testing"""
    np.random.seed(1)
    n_reads = 10_000_000
    n_snps = 100_000
    n_sfs = 40
    n_rgs = 20
    e = 0.001
    b = 0.001

    # psi = 1. / np.random.sample(n_snps)
    # psi /= np.max(psi)
    psi = np.random.sample(n_snps)
    # psi /= np.max(psi)

    true_tau = np.arange(0, 1, step=1 / n_sfs)
    true_cont = np.arange(0, 1, step=1 / n_rgs)

    SNP2SFS = np.random.randint(0, n_sfs, size=n_snps)
    # FLIPPED = np.zeros(n_snps, bool)

    G = np.random.binomial(2, true_tau[SNP2SFS], size=n_snps)
    A = np.random.binomial(1, psi, size=n_snps)

    # random assignments. In principle, a SNP might have repeat cats
    READ2RG = np.random.randint(0, n_rgs, size=n_reads)
    reps = (
        np.random.multinomial(n_reads - n_snps, pvals=np.zeros(n_snps) + 1.0 / n_snps)
        + 1
    )
    READ2SNP = np.repeat(np.arange(n_snps, dtype=np.int64), reps)

    p = (
        true_cont[READ2RG] * psi[READ2SNP]
        + (1 - true_cont[READ2RG]) * G[READ2SNP] / 2.0
    )
    p = e * (1.0 - p) + p * (1.0 - b)
    O = np.random.binomial(1, p, n_reads)

    FLIPPED = np.zeros(n_snps, bool)
    # FLIPPED = np.array(np.random.binomial(1, .4, size =n_snps), bool)
    # O[FLIPPED[READ2SNP]] = 1 - O[FLIPPED[READ2SNP]]

    np.set_printoptions(suppress=True, precision=3)

    data = SlugReads(
        READS=O,
        psi=psi,
        READ2RG=READ2RG,
        READ2SNP=READ2SNP,
        FLIPPED=FLIPPED,
        SNP2SFS=SNP2SFS,
    )
    pars = SlugPars(
        cont=np.repeat(0.5, n_rgs),
        tau=np.repeat(0.4, n_sfs),
        F=np.repeat(0.5, n_sfs),
        e=0.001,
        b=0.001,
    )
    # update_pars_reads(pars, data, data,
    #            True, True, True)
    # print(f'eb= {pars.e}, {pars.b}')
    # print(f'C = {pars.cont}')
    # print(f'F = {pars.F}')
    # print(f'tau = {pars.tau}')
    # update_pars_reads(pars, data, data,
    #            True, True, True)
    # print( f'll : {pars.prev_ll} -> {pars.ll}')
    controller = SlugController(
        update_eb=True, update_ftau=True, update_cont=True, n_iter=1000, ll_tol=1e-4
    )
    squarem(pars, data, controller)
    assert True
