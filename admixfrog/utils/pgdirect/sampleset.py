import itertools
from pybedtools import BedTool
from .snp import SNP
from .io import FastaFile, VCFFile, BamFile
from .sample import *
from collections import defaultdict
import pdb


class SampleSet(object):
    """
    Represents a set of samples that are processed concurrently 
    using common set of variants saved in pos

    This class is used to jointly iterate over several files.
    This can be done either using the ::meth:`~gt`-method or by directly
    iterating over the object, i.e.

    .. testsetup::

       from pgdirect import *
       from pybedtools import BedTool

    .. testcode:: 

       bed = BedTool("data/sample.bed")
       files = ["data/sample_{CHROM}.vcf.gz", "data/sample_spy.bam", 
           "data/chimp.fa.gz"]
       s = SampleSet.from_file_names(files, bed=bed, max_snps=10)

       for snp in s:
           if snp.is_biallelic:
               print(snp.chrom, snp.pos, "".join(sorted(snp.alleles)))

    .. testoutput::

       1 777121 AT
       11 300363 CT
       11 306790 CT
       11 306919 AT
       11 307243 CT


    Arguments
    ----------
    samples : list<File>
        A list of File object this FileSet contains
    bed : bed file handle
        An iterator over genomic positions to read. Can be constructed from a 
    ref : FastaFile
        Fasta File of desired reference allele
    mode : int
        Define the format of the alleles returned
        0 : bases(ACGT).
        1 : int(4)
        2 : int(2) (known alleles)
        3 : bialleiic (ref)
    """

    def __init__(
        self,
        samples,
        bed,
        files=None,
        ref=None,
        anc=None,
        max_snps=None,
        chrom=None,
        start=0,
        end=100000000000000,
        mode=0,
    ):
        #if type(bed) is BedTool:
        self.bed = bed
        #else:
        #    self.bed = BedTool(bed)
        self.samples = samples
        if files is None:
            self.files = set(s.file for s in self.samples)
        else:
            self.files = files
        self.ref = ref
        self.max_snps = max_snps

        self.chrom = chrom
        self.start = start
        self.end = end

        self.anc = anc
        assert self.anc is None or self.anc in self.samples

        self.reset_iter()

        # only get needed inds from files
        for file_ in self.files:
            file_.samples = file_.samples.intersection(self.samples)

    def reset_iter(self):
        self.snp_iterator = self.pos_from_bed(
            self.bed, self.max_snps, self.chrom, self.start, self.end
        )
        # set up joint iterator for all files
        file_iter = (file_.gt(self.pos_from_bed(self.bed)) for file_ in self.files)
        self.file_iter = zip(*file_iter)

    @classmethod
    def from_file_names(cls, file_names, chroms=None, *args, **kwargs):
        """from_file_names
        This is one of the ways a SampleSet object can be created. It takes in file names
        and returns all samples contained in this file. It's propbably the easiest to use,
        as no preliminary steps are necessary
        
        Arguments
        ----------
        file_names : list<File>
            the names of all the files being loaded. By default, all samples are
            loaded
        
        Returns
        -------
        sample_set : cls
            the SampleSet (or inherited) object being generated
        """
        files = set()
        for file_name in file_names:
            if file_name.endswith("bam"):
                files.add(BamFile(file_name, chroms))
            elif file_name.endswith("vcf") or file_name.endswith("vcf.gz"):
                files.add(VCFFile(file_name, chroms))
            elif file_name.endswith("fa") or file_name.endswith("fa.gz"):
                files.add(FastaFile(file_name, chroms))
        samples = set((sample for f in files for sample in f.samples))
        return cls(samples=samples, files=files, *args, **kwargs)

    @classmethod
    def from_samples(cls, samples, *args, **kwargs):
        """generate SampleSet object containing specified samples
        
        Args:
            samples (iter[Sample]) : samples to be contained in analysis
            *args : Further 
            **kwargs : type
        
        Returns
        -------
        sample_set : cls
            the SampleSet (or inherited) object being generated
        """
        files = set()
        samples = set(samples)
        for sample in samples:
            files.add(sample.file)
        return cls(samples=samples, files=files, *args, **kwargs)

    @staticmethod
    def pos_from_bed(bed, max_snps=None, chrom=None, start=0, end=1000000000):
        """create iterator over all positions in a bed file

        Parameters
        ----------
        bed : pybedtools.bed
            The bed file to be parsed
        max_snps : int or None
            The maximum number of SNPs. None reads entire file
        chrom : str or None
            Chromosome to be read. None reads all chromosomes
        start : int (default 0)
            starting position
        end : int (default 10^10)
            ending position
        
        Returns
        -------
        it : iterator
            iterator yielding tuple (chrom, pos) for each position to be parsed
        """
        if chrom is None:
            if 'pos' in bed: # custom file with chrom pos
                snp_iterator = zip(bed.chrom.astype(str), bed.pos-1)
            else:
                pos = (zip(itertools.cycle([b.chrom]), range(b.start, b.end)) for b in bed)
                snp_iterator = itertools.chain.from_iterable(pos)
        else:
            chrom = str(chrom)

            def bed_ints(bed):
                for b in bed:
                    if b.chrom == chrom:
                        if start <= b.start <= end or start <= b.end <= end:
                            s, e = max(start, b.start), min(end, b.end)
                            pos = zip(itertools.cycle([b.chrom]), range(s, e))
                            yield pos
                return

            snp_iterator = itertools.chain.from_iterable(bed_ints(bed))

        if max_snps is None:
            return snp_iterator
        else:
            return itertools.islice(snp_iterator, max_snps)

    @property
    def sample_names(self):
        l1 = [s.name for s in self.samples]
        return l1

    @property
    def n_samples(self):
        return len(self.samples)

    def __contains__(self, sample):
        return sample in self.samples

    def gt(self, verbose=False):
        """iterates over all genotypes in all samples
        """
        for i, ((chrom, pos), file_) in enumerate(
            zip(self.snp_iterator, self.file_iter)
        ):
            if i % 1000 == 0 and verbose:
                print(i)
                # print("samplesetiter", i, chrom, pos)
            yield SNP((chrom, pos), file_)

    def __iter__(self):
        return self.gt()


class BlockedSampleSet(SampleSet):
    """A SampleSet with output created in Blocks, i.e. sets of SNP or chromosome
    
    .. testsetup::

       from pgdirect import *
       from pybedtools import BedTool

    .. testcode:: 

       bed = BedTool("data/sample.bed")
       files = ["data/sample_{CHROM}.vcf.gz", "data/sample_spy.bam", 
           "data/chimp.fa.gz"]
       sblock = SingleBlocks(bed, max_snps=10)
       s = BlockedSampleSet.from_file_names(files, blocks=sblock)

       for block, snp in s:
           if snp.is_biallelic:
               print(block, snp.chrom, snp.pos, "".join(sorted(snp.alleles)))

    .. testoutput::

       4 1 777121 AT
       5 11 300363 CT
       6 11 306790 CT
       7 11 306919 AT
       8 11 307243 CT


    Parameters
    ----------
    blocks (:class: `pgdirect.Block`) :
        The :class: `Block`-setting to be used
    *args : further args passed to Sampleset
    """

    def __init__(self, blocks, *args, **kwargs):
        self.blocks = blocks
        super().__init__(bed=blocks.bed, *args, **kwargs)

    def gt(self):
        for i, ((block, snp), file_) in enumerate(zip(self.blocks, self.file_iter)):
            # for (block, snp), file_ in zip(self.blocks, self.file_iter):
            # if i % 1000 == 0 : print(i)
            # print("samplesetiter", i)
            yield block, SNP(snp, file_)


class CallBackSampleSet(BlockedSampleSet):
    """a sample set object that allows multiple things be done on the fly

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.callbacks = []

    def add_callback(self, cb):
        self.callbacks.append(cb)

    def run_callbacks(self):
        for cb in self.callbacks:
            cb.preprocess(self)

        for block, snp in self:
            for cb in self.callbacks:
                cb.process_snp(block, snp)

        for cb in self.callbacks:
            cb.postprocess(self)


__all__ = ["SampleSet", "BlockedSampleSet", "CallBackSampleSet"]
