from dataclasses import dataclass
import numpy as np
from typing import Any
from scipy.special import logit, expit


class SlugPars(object):
    """Parameters inferred by admixslug"""

    def __init__(self, cont, tau, F, e, b):
        self.par_names = ['cont', 'tau', 'F', 'e', 'b']

        lcs = locals()
        self.slices = dict()
        i = 0
        for par in self.par_names:
            try:
                new_i = i + len(lcs[par])
            except TypeError:
                new_i = i + 1
            self.slices[par] = slice(i, new_i)
            i = new_i

        self._pars = np.empty(i)
        for par in self.par_names:
            self._pars[ self.slices[par] ] = lcs[par]

        self.ll = np.NINF

        self.ll = np.NINF

        self.prev = np.zeros_like(self._pars)
        self.prev_ll = np.NINF

    def __setattr__(self, name, value):
        if name == 'par_names': 
            object.__setattr__(self, name, value)
        elif name in self.par_names:
            self._pars[self.slices[name]] = value
        elif name.startswith('prev_') and name != 'prev_ll':
            self.prev[self.slices[name[5:]]] = value
        else: 
            object.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name =='par_names':
            return super().__getattr__(name)
        elif name in self.par_names:
            return self._pars[self.slices[name]]
        elif name.startswith('prev_') and name != 'prev_ll':
            return self.prev[self.slices[name[5:]]] 
        elif name.startswith('delta_'):
            name = name[6:]
            if name == 'll':
                return self.ll - self.prev_ll
            ix = self.slices[name]
            return np.sum(np.abs(self.prev[ix] - self._pars[ix]))
        else:
            return super().__getattribute__(name)


    def __sub__(self, other):
        return self.pars - other.pars

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
    def delta_ll(self):
        return self.ll - self.prev_ll

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
