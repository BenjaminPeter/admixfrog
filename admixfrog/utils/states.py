from itertools import chain


class States(object):
    def __init__(
        self, state_ids, homo_states=None, het_states=None, 
        roh_states=None, est_inbreeding=False
    ):
        """get info about which states to include"""

        self.states = state_ids

        if homo_states is None:
            self.homo_ids = range(len(self.states))
        else:
            self.homo_ids = [i for i, s in enumerate(self.states) if s in homo_states]

        self.het_ids = []
        for i, s in enumerate(state_ids):
            for j, s2 in enumerate(state_ids[i + 1 :]):
                if het_states is None or s + s2 in het_states or s2 + s in het_states:
                    self.het_ids.append((i, j + i + 1))

        if est_inbreeding:
            if roh_states is None:
                self.roh_ids = self.homo_ids
            else:
                self.roh_ids = [i for i, s in enumerate(self.states) if s in roh_states]
        else:
            self.roh_ids = list()

        self.hap_ids = self.homo_ids

        self.est_inbreeding = est_inbreeding


    @property
    def homo(self):
        for i in self.homo_ids:
            yield i

    @property
    def het(self):
        for i in self.het_ids:
            yield i

    @property
    def roh(self):
        if self.est_inbreeding:
            for i in self.roh_ids:
                yield i

    @property
    def hap(self):
        for i in self.hap_ids:
            yield i

    @property
    def n_homo(self):
        return len(self.homo_ids)

    @property
    def n_het(self):
        return len(self.het_ids)

    @property
    def n_roh(self):
        return len(self.roh_ids)

    @property
    def n_hap(self):
        return len(self.hap_ids)


    @property
    def homo_names(self):
        for i in self.homo:
            yield self.states[i]

    @property
    def het_names(self):
        for i, j in self.het:
            yield self.states[i] + self.states[j]

    @property
    def roh_names(self):
        for i in self.roh:
            yield f'h{self.states[i]}'

    @property
    def hap_names(self):
        for i in self.hap:
            yield self.states[i]

    #n state and state names do not contian haploid states, as I reuse the
    #homozygous names
    @property
    def n_states(self):
        return self.n_homo + self.n_het + self.n_roh

    @property
    def state_names(self):
        return chain(self.homo_names, self.het_names, self.roh_names)

    @property
    def n_raw_states(self):
        """number of input/possible states"""
        return len(self.states)

    def __iter__(self):
        return self.states.__iter__()

    def __getitem__(self, i):
        return self.states[i]
