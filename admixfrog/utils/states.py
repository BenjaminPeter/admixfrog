from itertools import chain
import yaml


class States(object):
    """class to track different states"""

    def __init__(
        self,
        state_ids,
        homo_states=None,
        het_states=None,
        roh_states=None,
        est_inbreeding=False,
        ancestral = None,
        contamination = None,
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

        self.ancestral = ancestral
        self.contamination = contamination

    @classmethod
    def from_commandline(
        cls, raw_states, state_file, ancestral, cont_id=None, *args, **kwargs
    ):

        state_dict = cls.parse_state_string(
            raw_states + [ancestral, cont_id], state_file=state_file
        )
        state_dict2 = cls.parse_state_string(raw_states, state_file=state_file)
        states = list(dict(((x, None) for x in state_dict2.values())))
        new_states = cls(states, ancestral=ancestral, 
                         contamination=cont_id, *args, **kwargs)
        new_states.state_dict = state_dict
        return new_states

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
            yield f"h{self.states[i]}"

    @property
    def hap_names(self):
        for i in self.hap:
            yield self.states[i]

    # n state and state names do not contian haploid states, as I reuse the
    # homozygous names
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

    @staticmethod
    def parse_state_string(states, ext=[""], operator="+", state_file=None):
        """parse shortcut parse strings

        the basic way states are defined is as
            --states AFR NEA DEN

        this assmues the columns AFR, NEA and DEN are known and present in the reference
        We can rename and merge stuff here using the following syntax

        -- states AFR=YRI+SAN NEA=VIN DEN=Denisova2

        return rename dict for reference

        """
        d1 = [s.split("=") for s in states if s is not None]
        d2 = [(s if len(s) > 1 else (s[0], s[0])) for s in d1]
        state_dict = dict(
            (f"{k}{ext_}", f"{i}{ext_}")
            for i, j in d2
            for k in j.split(operator)
            for ext_ in ext
        )
        if state_file is not None:
            """logic here is:
            - yaml file contains some groupings (i.e. assign all
            Yorubans to YRI, all San to SAN
            - state_dict contains further groupings / renaming / filtering
                i.e. AFR=YRI+SAN
            - stuff not in state_dict will not be required

            therefore:
            1. load yaml
            2. get all required target pops from state_dict
            3. add all expansions replacements
            """
            Y = yaml.load(open(state_file), Loader=yaml.BaseLoader)
            # D = dict((v_, k) for (k,v) in Y.items() for v_ in v)
            s2 = dict()

            for state, label in state_dict.items():
                if state in Y:
                    for ind in Y[state]:
                        s2[ind] = label
                else:
                    s2[state] = label
            state_dict = s2
        return state_dict

    #@property
    #def ancestral(self):
        pass
