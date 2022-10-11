from itertools import chain

class States(object):
    def __init__(self, 
        state_ids,
        homo_states=None,
        het_states=None,
        est_inbreeding=False
        ):
        """get info about which states to include
        """

        self.states = state_ids


        if homo_states is None:
            self.homo_ids = range(len(self.states))
        else:
            self.homo_ids = [i for i, s in enumerate(self.states) if s in homo_states]


        self.het_ids = [] 
        for i, s in enumerate(state_ids):
            for j,s2 in enumerate(state_ids[i + 1 :]):
                if het_states is None or s+s2 in het_states or s2+s in het_states:
                    self.het_ids.append((i, j+i+1))


        if est_inbreeding:
            #self.state_names.extend((f'h{s}' for s in self.all_homo_states[homo_ids]))
            pass

        self.all_roh_states = [f'h{s}' for s in state_ids]
        self.all_hap_states = [s for s in state_ids] 
        self.hap_ids = range(len(self.all_hap_states))

    @property
    def homo(self):
        for i in self.homo_ids:
            yield i

    @property
    def het(self):
        for i in self.het_ids:
            yield i

    @property
    def n_homo(self):
        return len(self.homo_ids)

    @property
    def n_het(self):
        return len(self.het_ids)

    @property
    def n_hap(self):
        return len(self.states)

    @property
    def homo_names(self):
        for i in self.homo:
            yield self.states[i]

    @property
    def het_names(self):
        for i,j in self.het:
            yield self.states[i] + self.states[j]

    @property
    def n_states(self):
        return self.n_homo + self.n_het 

    @property
    def state_names(self):
        return chain(self.homo_names, self.het_names)

    def __iter__(self):
        return self.states.__iter__()

    def __getitem__(self, i):
        return self.states[i]
