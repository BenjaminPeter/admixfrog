from dataclasses import dataclass
import numpy as np
from typing import Any
from scipy.special import logit, expit


@dataclass(frozen=True)
class SlugData:
    """container for immutable data"""

    READS: np.ndarray  # 0 for REF, 1 for ALT call  [R x 1]

    psi: np.ndarray  # freq of contaminant [L x 1]

    READ2RG: np.ndarray  # which read is in which rg [R x 1]
    READ2SNP: np.ndarray  # which read belongs to which locus [R x 1]
    SNP2SFS: np.ndarray  # which snp belongs to which sfs entry [L x 1]
    FLIPPED: np.ndarray  # which SNP are flipped around [L x 1]

    states: Any = None  # the SFS references used to generate the data
    rgs: Any = None  # the read group names used
    chroms: Any = None  # the chromosome names used
    haplo_chroms: Any = None  # names of haploid chromosomes
    haploid_snps: Any = None  # which snp are haploid
    sex: str = "m"

    n_sfs: Any = None
    n_rgs: Any = None

    def __post_init__(self):
        assert len(self.READ2RG) == len(self.READ2SNP)
        assert len(self.READS) == len(self.READ2RG)
        assert len(self.psi) == len(self.SNP2SFS)
        assert len(self.FLIPPED) == len(self.SNP2SFS)
        super().__setattr__("READS", np.array(self.READS))
        super().__setattr__("psi", np.array(self.psi))
        super().__setattr__("FLIPPED", np.array(self.FLIPPED, "bool"))
        super().__setattr__("READ2RG", np.array(self.READ2RG))
        super().__setattr__("READ2SNP", np.array(self.READ2SNP))
        super().__setattr__("SNP2SFS", np.array(self.SNP2SFS))
        # give option to set n_sfs and n_rgs manually for resampling which may toss out some..
        if self.n_sfs is None:
            super().__setattr__("n_sfs", np.max(self.SNP2SFS) + 1)
        if self.n_rgs is None:
            super().__setattr__("n_rgs", np.max(self.READ2RG) + 1)
        self.READS.flags.writeable = False
        self.FLIPPED.flags.writeable = False
        self.psi.flags.writeable = False
        self.READ2RG.flags.writeable = False
        self.READ2SNP.flags.writeable = False
        self.SNP2SFS.flags.writeable = False

    def jackknife_sample(self, i, n_samples=20):
        """return a jackknife-subsample of the data"""
        RSLICE = np.ones(self.n_reads, bool)
        lsp = np.linspace(0, self.n_reads, n_samples + 1, dtype=int)
        RSLICE[lsp[i] : lsp[i + 1]] = 0
        old_snp_ids = np.unique(self.READ2SNP[RSLICE])  # old snp ids
        udict = dict(zip(old_snp_ids, range(len(old_snp_ids))))
        new_snp_ids = np.array(list(udict[snp] for snp in self.READ2SNP[RSLICE]))

        sr = SlugData(
            READS=self.READS[RSLICE],
            READ2SNP=new_snp_ids,
            READ2RG=self.READ2RG[RSLICE],
            SNP2SFS=self.SNP2SFS[old_snp_ids],  # [SSLICE],
            FLIPPED=self.FLIPPED[old_snp_ids],  # [SSLICE],
            haploid_snps=self.haploid_snps,
            states=self.states,
            rgs=self.rgs,
            psi=self.psi[old_snp_ids],
            chroms=self.chroms,
            haplo_chroms=self.haplo_chroms,
            n_sfs=self.n_sfs,
            n_rgs=self.n_rgs,
            sex=self.sex,
        )
        return sr

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
    """class for options for the admixslug em"""

    do_ll: bool = True
    update_eb: bool = True
    update_ftau: bool = True
    update_cont: bool = True
    update_delta: bool = False
    update_bias: bool = True
    update_F: bool = True
    n_iter: int = 200
    ll_tol: float = 0.01
    param_tol: float = 1e-4
    copy_pars: bool = False
    squarem_min: float = 1.0
    squarem_max: float = 1.0
    squarem_mstep: float = 2.0
    n_resamples: int = 10


@dataclass(frozen=True)
class FrogData:
    """container for immutable data"""

    O: np.ndarray  # number of non-reference alleles [R x 1]
    N: np.ndarray  # number of alleles [R x 1]

    psi: np.ndarray  # freq of contaminant [L x 1]

    alpha: np.ndarray  # alpha parameter for each SNP [L x 1]
    beta: np.ndarray  # beta parameter for each SNP [L x 1]
    alpha_hap: np.ndarray  # alpha parameter for each haploid SNP [L_h x 1]
    beta_hap: np.ndarray  # beta parameter for each haploid SNP [L_h x 1]

    OBS2RG: np.ndarray  # which read is in which rg [R x 1]
    OBS2SNP: np.ndarray  # which read belongs to which locus [R x 1]
    SNP2BIN: np.ndarray  # which snp belongs to which sfs entry [L x 1]
    n_bins: int

    bin_sizes: list

    states: Any = None  # the bins
    rgs: Any = None  # the read group names used
    chroms: Any = None  # the chromosome names used
    haplo_chroms: Any = None  # names of haploid chromosomes
    diplo_chroms: Any = None  # names of diploid chromosomes
    haploid_snps: Any = None  # which snp are haploid
    diploid_snps: Any = None  # which snp are diploid
    sex: str = "m"

    def __post_init__(self):
        assert len(self.OBS2RG) == len(self.OBS2SNP)
        # assert len(self.psi) == len(self.SNP2BIN)
        super().__setattr__("psi", np.array(self.psi))
        super().__setattr__("OBS2RG", np.array(self.OBS2RG))
        super().__setattr__("OBS2SNP", np.array(self.OBS2SNP))
        super().__setattr__("SNP2BIN", np.array(self.SNP2BIN))
        # give option to set n_sfs and n_rgs manually for resampling which may toss out some..
        self.psi.flags.writeable = False
        self.OBS2RG.flags.writeable = False
        self.OBS2SNP.flags.writeable = False
        self.SNP2BIN.flags.writeable = False

    @property
    def n_rgs(self):
        return len(self.rgs)

    @property
    def n_snps(self):
        return len(self.SNP2BIN)

    @property
    def n_obs(self):
        return len(self.OBS2SNP)

    @property
    def n_reads(self):
        return sum(self.N)

    @property
    def OBS2BIN(self):
        return self.SNP2BIN[self.OBS2SNP]

    pass


@dataclass
class FrogOptions:
    est_contamination: bool = True
    est_F: bool = True
    est_tau: bool = True
    est_trans: bool = True
    est_error: bool = True
    gt_mode: bool = False
    est_inbreeding: bool = True
    scale_probs: bool = True
    freq_contamination: int = 1
    freq_F: int = 1
    ll_tol: float = 0.01
    param_tol: float = 0.01
    copy_pars: bool = False
    squarem_min: float = 1.0
    squarem_max: float = 1.0
    squarem_mstep: float = 2.0
    n_iter: int = 100


class FrogX:
    def __init__(self, P):
        n_bins, n_states = P.n_bins, P.states.n_states
        n_snps = P.n_snps
        n_hap_states = P.states.n_hap

        self.Z = np.zeros((n_bins, n_states))  # P(Z | O)
        self.A = np.zeros((n_bins, n_states))  # fwd/bwd alpha
        self.B = np.zeros((n_bins, n_states))  # fwd/bwd beta
        self.N = np.zeros((n_bins))  # fwd/bwd beta
        self.E = np.ones((n_bins, n_states))  # P(O| Z)
        self.SNP = np.ones(
            (n_snps, n_states, 3)
        )  # P(O, G | Z), scales such that the max is 1
        self.PG = np.zeros((n_snps, n_states, 3))  # P(G, Z | O)

        # pointers to the same data, but split by chromosome
        self.H = FrogXHap()
        self.gamma, self.emissions = [], []
        self.alpha, self.beta, self.n = [], [], []
        self.H.gamma, self.H.emissions = [], []
        self.H.alpha, self.H.beta, self.H.n = [], [], []

        row0 = 0
        for r, chrom in zip(P.bin_sizes, P.chroms):
            if chrom in P.haplo_chroms:
                self.H.gamma.append(self.Z[row0 : (row0 + r), :n_hap_states])
                self.Z[row0 : (row0 + r), :n_hap_states] = 1 / n_hap_states
                self.H.emissions.append(self.E[row0 : (row0 + r), :n_hap_states])
                self.H.alpha.append(self.A[row0 : (row0 + r), :n_hap_states])
                self.H.beta.append(self.B[row0 : (row0 + r), :n_hap_states])
                self.H.n.append(self.N[row0 : (row0 + r)])
            else:
                self.gamma.append(self.Z[row0 : (row0 + r)])
                self.Z[row0 : (row0 + r)] = 1 / n_states
                self.emissions.append(self.E[row0 : (row0 + r)])
                self.alpha.append(self.A[row0 : (row0 + r)])
                self.beta.append(self.B[row0 : (row0 + r)])
                self.n.append(self.N[row0 : (row0 + r)])
            row0 += r

        self.scaling = 0


class FrogXHap:
    pass
