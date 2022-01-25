from dataclasses import dataclass
import numpy as np
from typing import Any
from scipy.special import logit, expit


@dataclass
class SlugPars(object):
    """Parameters inferred by admixslug
    """

    cont: np.ndarray  # contamination estimate per rg
    tau: np.ndarray  # tau est per SFS entry
    F: np.ndarray  # F est per SFS entry
    e: float   # error  0-> 1
    b: float   # refernce bias;  error 1 -> 0
    ll: float = np.NINF

    # self.delta : float = 0 #accuracy of contamination; currently unused

    def __post_init__(self):
        assert len(self.F) == len(self.tau)
        super().__setattr__("cont", np.array(self.cont))
        super().__setattr__("F", np.array(self.F))
        super().__setattr__("tau", np.array(self.tau))
        self.prev_tau = np.zeros_like(self.tau)
        self.prev_F = np.zeros_like(self.F)
        self.prev_cont = np.zeros_like(self.cont)
        self.prev_ll = np.NINF

    @property
    def n_sfs(self):
        return len(self.F)

    @property
    def n_rgs(self):
        return len(self.cont)

    @property
    def delta_cont(self):
        return np.sum(np.abs(self.prev_cont - self.cont))

    @property
    def delta_F(self):
        return np.sum(np.abs(self.prev_F - self.F))

    @property
    def delta_tau(self):
        return np.sum(np.abs(self.prev_tau - self.tau))

    @property
    def delta_ll(self):
        return self.ll - self.prev_ll

    @property
    def delta_e(self):
        return self.e - self.prev_e

    @property
    def delta_b(self):
        return self.b - self.prev_b

class SlugParsSquare(object):
    """Parameters inferred by admixslug
    """
    def __init__(self, cont, tau, F, e, b):
        n_rgs, n_sfs = len(cont), len(tau)
        assert len(tau) == len(F)


        k = n_rgs + 2 * n_sfs
        self.cont_slice = slice(n_rgs)
        self.tau_slice = slice(n_rgs, n_rgs + n_sfs)
        self.F_slice = slice(n_rgs + n_sfs, k)
        self.e_slice = slice(k, k+1)
        self.b_slice = slice(k+1, k+2)

        self._pars = np.empty(k + 2)

        self._pars[self.cont_slice] = cont
        self._pars[self.tau_slice] = tau 
        self._pars[self.F_slice] = F
        self._pars[self.b_slice] = b 
        self._pars[self.e_slice] = e
    
        self.ll = np.NINF


        self.prev_tau = np.zeros_like(self.tau)
        self.prev_F = np.zeros_like(self.F)
        self.prev_cont = np.zeros_like(self.cont)
        self.prev_ll = np.NINF


    def __sub__(self, other):
        return self.pars - other.pars
    def not__sub__(self, other): #replaced
        return SlugParsSquare(cont = self.cont - other.cont,
                              tau = self.tau - other.tau,
                              F = self.F - other.F,
                              e = self.e - other.e,
                              b = self.b - other.b
                              )

    @property
    def tau(self):
        return self._pars[self.tau_slice]

    @tau.setter
    def tau(self, tau):
        self._pars[self.tau_slice] = tau

    @property
    def F(self):
        return self._pars[self.F_slice]

    @F.setter
    def F(self, F):
        self._pars[self.F_slice] = F

    @property
    def b(self):
        return self._pars[self.b_slice]

    @b.setter
    def b(self, b):
        self._pars[self.b_slice] = b

    @property
    def e(self):
        return self._pars[self.e_slice]

    @e.setter
    def e(self, e):
        self._pars[self.e_slice] = e

    @property
    def cont(self):
        return self._pars[self.cont_slice]

    @cont.setter
    def cont(self, cont):
        self._pars[self.cont_slice] = cont

    @property
    def norm(self):
        return np.sqrt(np.sum(np.power(self._pars, 2)))

    @property
    def pars(self):
        return self._pars

    @property
    def n_sfs(self):
        return len(self.F)

    @property
    def n_rgs(self):
        return len(self.cont)

    @property
    def n_pars(self):
        return len(self._pars)

    @property
    def delta_cont(self):
        return np.sum(np.abs(self.prev_cont - self.cont))

    @property
    def delta_F(self):
        return np.sum(np.abs(self.prev_F - self.F))

    @property
    def delta_tau(self):
        return np.sum(np.abs(self.prev_tau - self.tau))

    @property
    def delta_ll(self):
        return self.ll - self.prev_ll

    @property
    def delta_e(self):
        return self.e - self.prev_e

    @property
    def delta_b(self):
        return self.b - self.prev_b

@dataclass(frozen=True)
class SlugData:
    """container for immutable data used for admixslug
    in total requires O x 4 + L x 1 integers of memory, and Lx1 float memory
    """

    REF: np.ndarray  # number of ref allele obs [O x 1]
    ALT: np.ndarray  # number of alt allele obs [O x 1]

    psi: np.ndarray  # freq of contaminant [L x 1]

    OBS2RG: np.ndarray  # which obs is in which rg [O x 1]
    OBS2SNP: np.ndarray  # which obs belongs to which locus [O x 1]
    SNP2SFS: np.ndarray  # which snp belongs to which sfs entry [L x 1]


    states : Any = None  #the SFS references used to generate the data
    rgs : Any = None    # the read group names used
    chroms : Any = None  #the chromosome names used
    haplo_chroms : Any = None #names of haploid chromosomes
    haploid_snps : Any = None #which snp are haploid
    sex : str = 'm'

    def __post_init__(self):
        assert len(self.REF) == len(self.ALT)
        assert len(self.OBS2RG) == len(self.OBS2SNP)
        assert len(self.REF) == len(self.OBS2RG)
        assert len(self.psi) == len(self.SNP2SFS)
        super().__setattr__("REF", np.array(self.REF))
        super().__setattr__("ALT", np.array(self.ALT))
        super().__setattr__("psi", np.array(self.psi))
        super().__setattr__("OBS2RG", np.array(self.OBS2RG))
        super().__setattr__("OBS2SNP", np.array(self.OBS2SNP))
        super().__setattr__("SNP2SFS", np.array(self.SNP2SFS))
        super().__setattr__("n_sfs", np.max(self.SNP2SFS) + 1)
        super().__setattr__("n_rgs", np.max(self.OBS2RG) + 1)
        assert 1 + np.max(self.OBS2SNP) == self.n_snps
        assert 1 + np.max(self.OBS2RG) == self.n_rgs
        self.REF.flags.writeable = False
        self.ALT.flags.writeable = False
        self.psi.flags.writeable = False
        self.OBS2RG.flags.writeable = False
        self.OBS2SNP.flags.writeable = False
        self.SNP2SFS.flags.writeable = False

    @property
    def O(self):
        return self.ALT

    @property
    def N(self):
        return self.ALT + self.REF

    @property
    def n_snps(self):
        return len(self.psi)

    @property
    def n_obs(self):
        return len(self.O)

    @property
    def n_reads(self):
        return np.sum(self.N)

    @property
    def OBS2SFS(self):
        return self.SNP2SFS[self.OBS2SNP]

@dataclass(frozen=True)
class SlugReads:
    """container for immutable data
    """

    READS: np.ndarray  # 0 for REF, 1 for ALT call  [R x 1]

    psi: np.ndarray  # freq of contaminant [L x 1]

    READ2RG: np.ndarray  # which read is in which rg [R x 1]
    READ2SNP: np.ndarray  # which read belongs to which locus [R x 1]
    SNP2SFS: np.ndarray  # which snp belongs to which sfs entry [L x 1]
    FLIPPED : np.ndarray # which SNP are flipped around [L x 1]


    states : Any = None  #the SFS references used to generate the data
    rgs : Any = None    # the read group names used
    chroms : Any = None  #the chromosome names used
    haplo_chroms : Any = None #names of haploid chromosomes
    haploid_snps : Any = None #which snp are haploid
    sex : str = 'm'

    def __post_init__(self):
        assert len(self.READ2RG) == len(self.READ2SNP)
        assert len(self.READS) == len(self.READ2RG)
        assert len(self.psi) == len(self.SNP2SFS)
        assert len(self.FLIPPED) == len(self.SNP2SFS)
        super().__setattr__("READS", np.array(self.READS))
        super().__setattr__("psi", np.array(self.psi))
        super().__setattr__("FLIPPED", np.array(self.FLIPPED, 'bool'))
        super().__setattr__("READ2RG", np.array(self.READ2RG))
        super().__setattr__("READ2SNP", np.array(self.READ2SNP))
        super().__setattr__("SNP2SFS", np.array(self.SNP2SFS))
        super().__setattr__("n_sfs", np.max(self.SNP2SFS) + 1)
        super().__setattr__("n_rgs", np.max(self.READ2RG) + 1)
        assert 1 + np.max(self.READ2SNP) <= self.n_snps
        assert 1 + np.max(self.READ2RG) <= self.n_rgs
        self.READS.flags.writeable = False
        self.FLIPPED.flags.writeable = False
        self.psi.flags.writeable = False
        self.READ2RG.flags.writeable = False
        self.READ2SNP.flags.writeable = False
        self.SNP2SFS.flags.writeable = False

    def jackknife_sample(self, i, n_samples=20):
        """return a jackknife-subsample of the data
        """
        RSLICE = np.ones(self.n_reads, bool)
        lsp = np.linspace(0, self.n_reads, n_samples+1, dtype=int)
        RSLICE[lsp[i]:lsp[i+1]] = 0

        return SlugReads(READS = self.READS[RSLICE],
                         READ2SNP = self.READ2SNP[RSLICE],
                         READ2RG = self.READ2RG[RSLICE],
                         SNP2SFS = self.SNP2SFS, #[SSLICE],
                         FLIPPED = self.FLIPPED, #[SSLICE],
                         haploid_snps=self.haploid_snps,
                         states = self.states,
                         rgs = self.rgs,
                         psi = self.psi,
                         chroms = self.chroms,
                         haplo_chroms = self.haplo_chroms,
                         sex = self.sex)
                         

    @property
    def n_snps(self):
        return len(self.psi)

    @property
    def n_reads(self):
        return len(self.READS)

    @property
    def READ2SFS(self):
        return self.SNP2SFS[self.READ2SNP]


@dataclass
class SlugController:
    """class for options for the admixslug em
    """

    do_ll: bool = True
    update_eb: bool = True
    update_ftau: bool = True
    update_cont: bool = True
    update_delta: bool = False
    update_bias: bool = True
    update_F : bool = True
    n_iter: int = 200
    ll_tol: float = 0.01
    param_tol : float = 1e-4
    copy_pars : bool = False
    squarem_min : float = 1.
    squarem_max : float = 1.
    squarem_mstep : float = 2.
    n_resamples : int = 10
