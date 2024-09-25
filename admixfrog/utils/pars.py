import numpy as np


class Pars(object):
    def __init__(self, **kwargs):
        lcs = kwargs
        self.slices = dict()
        self.shapes = dict()
        i = 0
        for par in self.par_names:
            arr = np.array(lcs[par])
            self.shapes[par] = arr.shape
            new_i = i + len(arr.flatten())
            self.slices[par] = slice(i, new_i)
            i = new_i

        self._pars = np.empty(i)
        for par in self.par_names:
            arr = np.array(lcs[par])
            self._pars[self.slices[par]] = arr.flatten()

        self.ll = -np.inf

        self.prev = np.zeros_like(self._pars)
        self.prev_ll = -np.inf

    def __setattr__(self, name, value):
        """get attributes by names.
        I use this function as a more readable (but slightly slower) version of previosu extended one.
        Key idea is that all parameters are stored in a flat numpy array for squarem/em
        """
        if name == "par_names":
            object.__setattr__(self, name, value)
        elif name in self.par_names:
            self._pars[self.slices[name]] = value.flatten()
        elif name.startswith("prev_") and name != "prev_ll":
            self.prev[self.slices[name[5:]]] = value.flatten()
        else:
            object.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name == "par_names":
            return super().__getattr__(name)
        elif name in self.par_names:
            return self._pars[self.slices[name]].reshape(self.shapes[name])
        elif name.startswith("prev_") and name != "prev_ll":
            return self.prev[self.slices[name[5:]]].reshape(self.shapes[name[5:]])
        elif name.startswith("delta_"):
            name = name[6:]
            if name == "ll":
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
    def n_pars(self):
        return len(self._pars)

    @property
    def delta_ll(self):
        return self.ll - self.prev_ll


class SlugPars(Pars):
    """Parameters inferred by admixslug"""

    par_names = ["cont", "tau", "F", "e", "b"]

    @property
    def n_sfs(self):
        return len(self.F)

    @property
    def n_rgs(self):
        return len(self.cont)

    @property
    def n_pars(self):
        return len(self._pars)


class FrogPars(Pars):
    """Parameters inferred by admixfrog"""

    par_names = [
        "alpha0",
        "alpha0_hap",
        "trans",
        "trans_hap",
        "cont",
        "error",
        "F",
        "tau",
    ]

    class _H:
        pass

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.H = self._H()
        self.H.alpha0 = self.alpha0_hap
        self.H.trans = self.trans_hap

    @property
    def n_states(self):
        return len(self.F)

    @property
    def n_rgs(self):
        return len(self.cont)
