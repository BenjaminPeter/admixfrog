import random, string
from .gt import ReadGT, DiploidGT, HaploidGT, GTNoData


class Sample(object):
    """
    represents a single individual
    this class has two main purposes. First, it allows setting up populations.
    Second, it defines how stuff is sampled
    """

    def __init__(self, file_, sample_name, populations=[]):
        if type(file_) is str:
            from .io import BamFile, VCFFile, FastaFile

            file_name = file_
            if file_name.endswith("bam"):
                self._file = BamFile(file_name, sample_name)
            elif file_name.endswith("vcf") or file_name.endswith("vcf.gz"):
                self._file = VCFFile(file_name, sample_name)
            elif file_name.endswith("fa") or file_name.endswith("fa.gz"):
                self._file = FastaFile(file_name, sample_name)
            else:
                raise ValueError("cant initiate sample with file name %s" % file_name)
        else:
            self._file = file_
        self._name = sample_name
        self._outname = sample_name
        self.populations = populations
        self._subsamples = []
        for pop in populations:
            pop.add(self)

    @property
    def subsamples(self):
        return self._subsamples

    @staticmethod
    def from_file(file_name):
        from .io import BamFile, VCFFile, FastaFile

        if file_name.endswith("bam"):
            file_ = BamFile(file_name)
        elif file_name.endswith("vcf") or file_name.endswith("vcf.gz"):
            file_ = VCFFile(file_name)
        elif file_name.endswith("fa") or file_name.endswith("fa.gz"):
            file_ = FastaFile(file_name)
        else:
            raise ValueError("cant _handle filetype")
        assert len(file_.samples) == 1
        s = next(iter(file_.samples))
        return s

    @property
    def name(self):
        """name is the internal name of the sample, that is immutable.
        in a vcf file, this corresponds to the identifier;
        """
        return self._name

    @property
    def outname(self):
        """name is the output name of the sample, that is mutable."""
        return self._outname

    @outname.setter
    def name(self, name):
        self._outname = name

    def __str__(self):
        return "<%s %s>" % (type(self), self._outname)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self._name)

    def __eq__(self, other):
        return hash(self) == hash(other)

    @property
    def file(self):
        return self._file

    @property
    def has_subsamples(self):
        return len(self._subsamples) > 0

    def __iter__(self):
        if self.has_subsamples:
            for ss in self.subsamples:
                yield ss
        else:
            yield self
        return


class SampleWithSubsamples(Sample):
    """ """

    def __init__(
        self,
        sample,
    ):
        self.sample = sample

    @property
    def file(self):
        return self.sample._file

    @property
    def name(self):
        return "|".join((self.sample.name, str(self.rg), str(self.deam)))


class Population:
    """group individuals in a population for calculations"""

    def __init__(self, samples=[], subpopulations=[], name=None, **kwargs):
        self._samples = samples
        self.subpopulations = subpopulations
        self.superpopulations = []
        if name is None:
            rstr = "".join(random.choices(string.ascii_uppercase, k=6))
            self.name = "Pop" + rstr
        else:
            self.name = name

        for sample in self._samples:
            self.add(sample)

        for subpop in subpopulations:
            self.add_subpop(subpop)

        # additional args can be stored, unused for now
        self.__dict__.update(kwargs)

    def add(self, sample):
        if sample not in self._samples:
            self._samples.append(sample)
        if self not in sample.populations:
            sample.populations.append(self)

    def add_subpop(self, subpop):
        if subpop not in self.subpopulations:
            self.subpopulations.append(subpop)
        if self not in subpop.superpopulations:
            subpop.superpopulations.add(self)

    @property
    def samples(self):
        """samples
        returns all samples from self and all subpopulations
        """
        s = [subpop.samples for subpop in self.subpopulations]
        return self._samples + s

    def __iter__(self):
        for s in self.samples:
            yield s


__all__ = ["Sample", "SampleWithSubsamples", "Population"]
