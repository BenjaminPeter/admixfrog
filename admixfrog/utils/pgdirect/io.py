from abc import ABC, abstractmethod
from collections import defaultdict, namedtuple
import pysam
from . import DEAMINATION, SAMPLING
from .gt import ReadGT, DiploidGT, HaploidGT, GTNoData
import os.path
from .sample import Sample
from math import ceil
from .utils import rev_complement, Read



class File(object):
    """Abstract File Handling Class

    Base class for handling streaming various genomic data sets.
    Most likely, one of the derived :class:`VCFFile`, :class:`BamFile` or 
    :class:`FastaFile` will
    be used.

    Args:
        fname (str) : name of file to be read

    Attributes:
        fname (str) : name of file to be read

    """

    def __init__(self, fname, chroms=None, sample_name=None):
        self.fname = fname
        self.chroms = chroms
        if sample_name is None:
            sample_name = os.path.basename(fname)
        self._samples = {Sample(self, sample_name)}
        #        self.subsamples = [("GT", None) ]

    def __str__(self):
        return "<File: " + self.fname + ">"

    def __repr__(self):
        return "<File: " + self.fname + ">"

    @property
    def sample(self):
        """the sample object if file contains exactly one sample

        Returns
        -------
        <Sample> object if it is unique
        """
        if len(self.samples) == 1:
            return next(iter(self.samples))
        raise ValueError("File contains more than one sample")

    @property
    def samples(self):
        """a list of samples contained in this file"""
        return self._samples

    @samples.setter
    def samples(self, new_samples):
        """a list of samples contained in this file"""
        self._samples = new_samples


class ComplexFile(File, ABC):
    """Abstract class for jointly iterating bed file and data file


    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._handle = None

    @property
    @abstractmethod
    def _offset(self):
        """_offset
        the offset for accessing positions compared to zero-based indexing
        (O for bam, 1 for vcf)
        """
        return 0

    @abstractmethod
    def _open_file(self, chrom):
        pass

    @abstractmethod
    def process_snp(self, snp_data):
        """process one snp of data
        """
        pass

    def _close_file(self):
        self._handle.close()

    def no_data_default(self, chrom, pos):
        return GTNoData(self.sample, chrom, pos)

    def gt(self, target_positions):
        cur_chrom = None
        site = None

        for i, (target_chrom, target_pos) in enumerate(target_positions):

            # we found last target position
            if target_chrom != cur_chrom:
                print(target_chrom)
                if cur_chrom is not None:
                    print("closing", target_chrom)
                    self._close_file()

                sites = self._open_file(target_chrom)

                cur_chrom = target_chrom
                try:
                    site = next(sites)
                except StopIteration:  # no sites in bam
                    site = None

            # case 0: no data for chromosome
            if site is None:
                yield self.no_data_default(target_chrom, target_pos)
                #print("no chrom data", target_chrom, target_pos)
                continue

            # case 1: pos in bam not in target, advance bam
            while (
                target_pos > site.pos - self._offset
            ):  # and target_chrom == site.reference_name:
                try:
                    site = next(sites)
                except StopIteration:  # no more sites in bam,
                    # print("stopping bam", target_pos, site.pos, target_chrom)
                    yield self.no_data_default(
                        target_chrom, target_pos
                    )  # taret snp after last snp in bam file
                    break

            # case 2: target SNP not in the bam, we return NoData and advance target list
            if target_pos < site.pos - self._offset:
                yield self.no_data_default(target_chrom, target_pos)
                # print("no data", (chrom, pos), (site.reference_name, site.pos), cur_chrom)
                continue

            if target_pos == site.pos - self._offset:  # and chrom == col.contig:
                # print("gt_iter", i, target_pos)
                yield self.process_snp(site)

        return


class VCFFile(ComplexFile):
    """VCFFile
    use for VCF with multiple individuals to be parsed
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.chroms is None:
            fname = self.fname.format(CHROM="1")
        else:
            fname = self.fname.format(CHROM=chroms[0])

        with pysam.VariantFile(fname) as h:
            self._samples = set()
            for sample_name in h.header.samples:
                self._samples.add(Sample(self, sample_name))

    def no_data_default(self, chrom, pos):
        f = lambda s: GTNoData(s, chrom, pos)
        return ((s, f(s)) for s in self.samples)

    def _open_file(self, chrom):
        self._handle = pysam.VariantFile(self.fname.format(CHROM=chrom))
        sample_list = [s.name for s in self.samples]
        self._handle.subset_samples(sample_list)
        try:
            vcf = self._handle.fetch(contig=chrom)
        except ValueError:  #chrom not found
            return iter([])
        return vcf

    def process_snp(self, snp):
        def gt_iter():
            for sample in self.samples:
                # pdb.set_trace()
                gt = DiploidGT(
                    snp.samples[sample.name].alleles,
                    sample,
                    snp.chrom,
                    snp.pos - self._offset,
                )
                yield sample, gt

        return gt_iter()

    @property
    def _offset(self):
        return 1

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, new_samples):
        """we can only subset vcfs, not add new sample
        better implementation would separate samples and all possible samples
        """
        for new_sample in new_samples:
            assert new_sample in self._samples
        self._samples = new_samples


class BamFile(ComplexFile):
    """read and process BAM-files

    Note
    ------
    note that pysam is ZERO-indexed; i.e. the positions are one less than
    what would be expected from e.g. running samtools pileup or checking
    ucsc. 

    i.e. the first column of a bed file is what we need


    samples : ["single", "rg", "none"]
        how samples are defined. If single, a single sample is created.
        if rg, a sample per read group is created
        for other values, no samples are created and this needs to be done manually
    """

    def __init__(self, *args, samples="single", **kwargs):
        super().__init__(*args, **kwargs)
        if self.chroms is None:
            bam = pysam.AlignmentFile(self.fname.format(CHROM="1"))
        else:
            bam = pysam.AlignmentFile(self.fname.format(CHROM=self.chroms[0]))

        if "RG" in bam.header:
            self.rgs = set(i["ID"] for i in bam.header["RG"])
        else:
            self.rgs = []

        self.cleanup_rgs()
        bam.close()

        self.n_indel_rm = 0

        if samples == "single":
            """handled by default constructor"""
            print(self.sample, len(self.samples))
            self.rg2sample = defaultdict(lambda: self.sample)
        elif samples == "rg":
            self.samples = set()
            self.rg2sample = dict()

            for rg in self.rgs:
                sample = Sample(self, rg)
                self.samples.add(sample)
                self.rg2sample[rg] = sample

    def cleanup_rgs(self):
        """get simplified read group information

        some heuristic to clean up read groups learned from file.
        currently, all readgroups starting with Phi (spike-in?)
        are removed, and rgs starting of form "P???-12345" are
        shortened to P???
        """
        self.rgs = list(
            set(self.cleanup_rg(rg) for rg in self.rgs if not rg.startswith("Phi"))
        )

    @staticmethod
    def cleanup_rg(rg):
        """simplified read-group
        
        Parameters
        ----------
        rg : str
            read-group name to be simplified
            
        
        Returns
        -------
        simplified rg (underscore replaced by dash)
        """
        return rg.replace("_", "-")

    def _open_file(self, chrom):
        self._handle = pysam.AlignmentFile(self.fname.format(CHROM=chrom))
        pileup = self._handle.pileup(reference=chrom, multiple_iterators=False)
        print(f"opening {chrom}")
        return pileup

    def process_snp(self, col):
        """process a single pileup column"""
        # pileup_dict is a dict[sample : read]
        pileup_dict = self.process_pileup(col)
        if len(pileup_dict) == 0:
            # print("returned no data")
            return ((sample, GTNoData()) for sample in self.samples)
        else:
            return ((sample, ReadGT(data)) for sample, data in pileup_dict.items())

    def no_data_default(self, chrom, pos):
        return ((s, GTNoData()) for s in self.samples)

    @property
    def _offset(self):
        return 0

    # from Martin, added R
    # def process_pileup(column, ref_base):
    def process_pileup(self, column):
        """Process pileup column

        This function extracts some basic info from aligned reads at a
        particular position.
        For each read, it returns the tuple
        (rg, deam), base, len
        where 

            - rg is the read group
            - deam is the deaination pattern
            - base is the particular base at this position
            - len is the total length of the read
        

        Arguments
        ---------
        column : pysam.Pileup.column
            all reads overlapping a single position


        Returns
        -------
        list(tuple)
        """
        pileup = defaultdict(list)
        # walk through all reads overlapping the current column
        # and accumulate bases at that position
        for pileup_read in column.pileups:
            r = self.process_read(pileup_read)
            if r is not None:
                pileup[self.rg2sample[r.RG]].append(r)

        return pileup

    def process_read(self, pileup_read):
        if pileup_read.is_del:
            return None

        A = pileup_read.alignment

        if "I" in A.cigarstring or "D" in A.cigarstring or "N" in A.cigarstring:
            self.n_indel_rm += 1
            return None
        if "S" in A.cigarstring:
            return None

        mapq = A.mapping_quality
        ref_seq = A.get_reference_sequence()
        pos_in_read = pileup_read.query_position
        read_len = A.query_length
        read_base = A.query_sequence[pos_in_read]
        try:
            ref_base = ref_seq[pos_in_read].upper()
        except IndexError:
            #print("no ref", ref_seq, pos_in_read, A.cigarstring)
            #print("no re2", A.query_sequence, pos_in_read, A.cigarstring)
            return None

        try:
            read_group = BamFile.cleanup_rg(A.get_tag("RG"))
        except KeyError:
            try:
                #guess RG from XI XJ tag, which in some bam files are p7/p5 ix
                read_group = f'{A.get_tag("XI")}-{A.get_tag("XJ")}'
            except KeyError:
                read_group = "NONE"
        try:
            baseq = ord(A.qual[pos_in_read]) - 33
        except UnicodeDecodeError:
            baseq = 0

        # deam_type = BamFile.get_deamination_type(A, ref_seq, A.is_reverse)
        d5, d3, deam_full = BamFile.get_deamination_pattern(A, ref_seq)
        deam_pattern = d5, d3

        if read_base in "ACGT":
            tpl = Read(
                read_group,
                deam_pattern,
                read_base,
                read_len,
                baseq,
                mapq,
                A.is_reverse,
                pos_in_read,
                deam_full
            )
            return tpl
        else:
            return None

    @staticmethod
    def get_deamination_pattern(A, ref_seq):
        fwd_iter = zip(ref_seq, A.query_sequence)
        bwd_iter = zip(ref_seq[::-1], A.query_sequence[::-1])
        if A.is_reverse:
            deam5 = [i for i, (r, q) in enumerate(fwd_iter) if r == "g" and q == "A"]
            deam3 = [i for i, (r, q) in enumerate(bwd_iter) if r == "g" and q == "A"]
        else:
            deam5 = [i for i, (r, q) in enumerate(fwd_iter) if r == "c" and q == "T"]
            deam3 = [i for i, (r, q) in enumerate(bwd_iter) if r == "c" and q == "T"]

        if len(deam5) == 0 or len(deam3) == 0:
            return -1, -1, []

        return deam5[0], deam3[0], deam5


class BamFileFilter(BamFile):
    def __init__(self, outfile, min_length=0, minq=0, pos_in_read_cutoff=0,
                 chroms=None,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if chroms is None:
            bam = pysam.AlignmentFile(self.fname.format(CHROM="1"))
        else:
            bam = pysam.AlignmentFile(self.fname.format(CHROM=chroms[0]))
        self.outfile = pysam.AlignmentFile(outfile, 'wb', template=bam)
        self.min_length = min_length
        self.minq = minq
        self.pos_in_read_cutoff = pos_in_read_cutoff
        bam.close()

    def process(self, target_positions):
        cur_chrom = None
        self.handle = pysam.AlignmentFile(self.fname)
        for target_chrom, target_pos in target_positions:
            if target_chrom != cur_chrom:
                reads = self.handle.fetch(target_chrom)
                cur_chrom = target_chrom
                try:
                    read = next(reads)
                except StopIteration: #no sites in bam
                    break

            # case 1: pos in bam not in target, advance bam
            while (
                target_pos > read.pos - self._offset + read.qlen
            ):  # and target_chrom == site.reference_name:
                try:
                    read = next(reads)
                except StopIteration:  # no more sites in bam,
                    # print("stopping bam", target_pos, site.pos, target_chrom)
                    yield self.no_data_default(
                        target_chrom, target_pos
                    )  # taret snp after last snp in bam file
                    break

            if target_pos + self._offset in read.positions:  # and chrom == col.contig:
                # print("gt_iter", i, target_pos)
                yield self.process_snp(read)
                read = next(reads)
                

    def process_snp(self, A):
        if "I" in A.cigarstring or "D" in A.cigarstring or "N" in A.cigarstring:
            return -4
        if A.alen < self.min_length:
            return -1
        if A.mapping_quality < self.minq:
            return -2
        self.outfile.write(A)
        return 1 

    def no_data_default(self, *args, **kwargs):
        return -3

    def close(self):
        self.outfile.close()


class FastaFile(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def gt(self, target_positions):
        cur_chrom = None
        fa = None
        for chrom, pos in target_positions:
            if cur_chrom != chrom:
                if cur_chrom is not None:
                    fa.close()
                fa = pysam.FastaFile(self.fname.format(CHROM=chrom))
                gt = HaploidGT(fa.fetch(chrom, pos, pos + 1), self.sample, chrom, pos)
            yield (i for i in [(self.sample, gt)])


__all__ = ["File", "ComplexFile", "VCFFile", "BamFileFilter", "BamFile", "FastaFile"]
