from dataclasses import dataclass
import numpy as np

class FrogPars(object):
    """Parameters inferred by admixfrog"""

    def __init__(self, alpha0, alpha0_hap, trans, trans_hap, cont, error, F, tau):
        par_names = ['alpha0', 'alpha0_hap', 'trans', 'trans_hap',
                     'cont', 'error', 'F', 'tau']

        self.par_names = par_names

        lcs = locals()
        self.slices = dict()
        i = 0
        for par in par_names:
            self.slices[par] = slice(i, i + len(lcs[par]))
            i += len(lcs[par])

        self._pars = np.empty(i)
        for par in par_names:
            self._pars[ self.slices[par] ] = lcs[par]

        self.ll = np.NINF



    def __sub__(self, other):
        return self.pars - other.pars

    def __setattr__(self, name, value):
        if name == 'par_names': 
            self.__dict__['par_names'] = value
        elif name in self.par_names:
            self._pars[self.slices[name]] = value
        else: 
            self.__dict__[name] = value

    def __getattr__(self, name):
        if name in self.par_names:
            return self._pars[self.slices[name]]
        else:
            return self.__dict__[name]


    @property
    def norm(self):
        return np.sqrt(np.sum(np.power(self._pars, 2)))

    @property
    def pars(self):
        return self._pars

    @property
    def n_states(self):
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
class FrogData:
    pass


@dataclass
class FrogController:
    pass

