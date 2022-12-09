from abc import ABC, abstractmethod
import lzma
from collections import defaultdict, namedtuple


class Callback(ABC):
    """
    callbacks are objects with three methods

        -preprocess: stuff that is done before iterations starts
        -postprocess: stuff that is done after we're done iterating

        -process_snp: functions that accept a tuple
            - block
            - SNP
    """

    @abstractmethod
    def preprocess(self, sampleset):
        pass

    @abstractmethod
    def postprocess(self, sampleset):
        pass

    @abstractmethod
    def process_snp(self, block, snp):
        pass


class Coverage(Callback):
    """get coverage with reference/ancestral allele from bam file"""

    def __init__(self, outfile, **kwargs):
        self.outfile = outfile
        self.kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = lzma.open(self.outfile, "wt")
        print(
            "chrom",
            "pos",
            "ref",
            "alt",
            "tref",
            "talt",
            "tdeam",
            "tother",
            sep=",",
            file=self.f,
        )

    def process_snp(self, block, snp):
        reads = snp.reads(**self.kwargs)
        n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        for r in reads:
            if r.base == snp.ref:
                n_ref += 1
            elif r.base == snp.alt:
                n_alt += 1
            elif r.base == "T" and not r.is_reverse and "C" in (snp.ref, snp.alt):
                n_deam += 1
            elif r.base == "A" and r.is_reverse and "G" in (snp.ref, snp.alt):
                n_deam += 1
            elif r.base != "N":
                n_other += 1
        print(
            snp.chrom,
            snp.pos,
            snp.ref,
            snp.alt,
            n_ref,
            n_alt,
            n_deam,
            n_other,
            file=self.f,
            sep=",",
        )

    def postprocess(self, sampleset):
        self.f.close()


class ExtCoverage(Callback):
    """get coverage with reference/ancestral allele from bam file"""

    class Obs:
        def __init__(self):
            self.n_ref = 0
            self.n_alt = 0
            self.n_deam = 0
            self.n_other = 0

    def __init__(self, outfile, deam_cutoff, **kwargs):
        self.outfile = outfile
        self.deam_cutoff = deam_cutoff
        self.kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = open(self.outfile, "w")
        print(
            "chrom",
            "pos",
            "rg",
            "deam",
            "ref",
            "alt",
            "deam",
            "other",
            sep=",",
            file=self.f,
        )

    def process_snp(self, block, snp):
        reads = snp.reads(**self.kwargs)
        D = defaultdict(lambda: self.Obs())
        # n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        for r in reads:
            DEAM = (
                r.deam[0] < self.deam_cutoff or r.deam[1] < self.deam_cutoff
            ) and r.deam[0] >= 0
            if r.base == snp.ref:
                D[r.RG, DEAM].n_ref += 1
            elif r.base == snp.alt:
                D[r.RG, DEAM].n_alt += 1
            elif r.base == "T" and not r.is_reverse and "C" in (snp.ref, snp.alt):
                D[r.RG, DEAM].n_deam += 1
            elif r.base == "A" and r.is_reverse and "G" in (snp.ref, snp.alt):
                D[r.RG, DEAM].n_deam += 1
            elif r.base != "N":
                D[r.RG, DEAM].n_other += 1
        for (rg, deam), r in D.items():
            print(
                snp.chrom,
                snp.pos,
                rg,
                deam,
                r.n_ref,
                r.n_alt,
                r.n_deam,
                r.n_other,
                file=self.f,
                sep=",",
            )

    def postprocess(self, sampleset):
        self.f.close()


class VCFWriter(Callback):
    def __init__(self, outfile, **kwargs):
        self.outfile = outfile
        self.write_kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = open(self.outfile, "w")
        self.print_vcf_header(self.f, sampleset.samples)

    def process_snp(self, block, snp):
        rec = snp.to_vcfrecord(self.samples, **self.write_kwargs)
        # print("writing", rec)
        self.f.write(rec)

    def postprocess(self, sampleset):
        self.f.close()

    def print_vcf_header(self, handle, samples):
        """Print VCF header.; from Martin"""
        s = "##fileformat=VCFv4.1\n"
        s += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        s += '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">\n'
        s += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        self.samples = [s for s in samples]  # ensure ordering
        sample_names = [s.name for s in samples]
        s += "\t".join(sample_names)
        print(s, file=handle)

    def write(self, **kwargs):
        with open(self.outfile, "w") as f:
            samples = self.print_vcf_header(f)
            for _, snp in self:
                rec = snp.to_vcfrecord(samples, **kwargs)
                f.write(rec)


__all__ = ["VCFWriter", "Callback", "Coverage", "ExtCoverage"]
