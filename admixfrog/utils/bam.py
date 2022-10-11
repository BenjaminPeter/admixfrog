from . import pgdirect as pg
import logging
import pandas as pd
from collections import defaultdict
import lzma

default_filter = {
    "deam_only": False,
    "pos_in_read_cutoff": 2,
    "min_length": 35,
    "max_length": 1000,
    "minq": 25,
}


class AdmixfrogInput(pg.ExtCoverage):
    def __init__(
        self,
        outfile,
        deam_cutoff,
        length_bin_size=None,
        random_read_sample=False,
        report_alleles=False,
        max_reads=100,
        **kwargs,
    ):
        self.outfile = outfile
        self.deam_cutoff = deam_cutoff
        self.length_bin_size = length_bin_size
        self.random_read_sample = random_read_sample
        self.report_alleles = report_alleles
        self.max_reads = max_reads
        if random_read_sample:
            raise NotImplementedError
        try:
            self.min_length = kwargs["min_length"]
        except KeyError:
            self.min_length = 35
        self.kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = lzma.open(self.outfile, "wt")
        print(
            "chrom",
            "pos",
            "lib",
            "tref",
            "talt",
            "tdeam",
            "tother",
            sep=",",
            file=self.f,
        )

    def process_snp(self, block, snp):
        reads = snp.reads(**self.kwargs)
        D = defaultdict(lambda: self.Obs())
        # n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        i = 0
        for r in reads:
            i += 1
            if i > self.max_reads:
                print(
                    f"Warning at {snp.chrom}:{snp.pos+1} more than {i-1} reads found ({len(reads)}). Skipping the rest..."
                )
                break

            DEAM = (
                "deam"
                if (r.deam[0] < self.deam_cutoff or r.deam[1] < self.deam_cutoff)
                and r.deam[0] >= 0
                else "nodeam"
            )
            if self.length_bin_size is None:
                LEN = 0
            else:
                LEN = (r.len - self.min_length) // self.length_bin_size
            if r.base == snp.ref:
                D[r.RG, DEAM, LEN].n_ref += 1
            elif r.base == snp.alt:
                D[r.RG, DEAM, LEN].n_alt += 1
            elif r.base == "T" and not r.is_reverse and "C" in (snp.ref, snp.alt):
                D[r.RG, DEAM, LEN].n_deam += 1
            elif r.base == "A" and r.is_reverse and "G" in (snp.ref, snp.alt):
                D[r.RG, DEAM, LEN].n_deam += 1
            elif r.base != "N":
                D[r.RG, DEAM, LEN].n_other += 1
        for (rg, deam, len_), r in D.items():
            if self.length_bin_size is None:
                # lib = f"{rg}_{deam}"
                lib = f"{rg}_{deam}"
            else:
                # lib = f"{rg}_{len_}_{deam}"
                lib = f"{rg}_{len_}_{deam}"
            if self.report_alleles:
                # breakpoint()
                alleles = "".join(sorted(snp.ref + snp.alt))
                lib = f"{lib}_{alleles}"
            print(
                snp.chrom,
                snp.pos + 1,
                lib,
                r.n_ref,
                r.n_alt,
                r.n_deam,
                r.n_other,
                file=self.f,
                sep=",",
            )


class AdmixfrogInput2(pg.ExtCoverage):
    def __init__(
        self,
        outfile,
        deam_cutoff,
        length_bin_size=None,
        random_read_sample=False,
        report_alleles=False,
        **kwargs,
    ):
        self.outfile = outfile
        self.deam_cutoff = deam_cutoff
        self.length_bin_size = 1
        self.random_read_sample = random_read_sample
        self.report_alleles = report_alleles
        if random_read_sample:
            raise NotImplementedError
        try:
            self.min_length = kwargs["min_length"]
        except KeyError:
            self.min_length = 35
        self.kwargs = kwargs

    def preprocess(self, sampleset):
        self.f = lzma.open(self.outfile, "wt")
        print(
            "chrom",
            "pos",
            "tref",
            "talt",
            "lib",
            "len",
            "deam",
            "dmgsite",
            sep=",",
            file=self.f,
        )

    def process_snp(self, block, snp):
        reads = snp.reads(**self.kwargs)
        D = defaultdict(lambda: self.Obs())
        # n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        for r in reads:
            DEAM = (
                "deam"
                if (r.deam[0] < self.deam_cutoff or r.deam[1] < self.deam_cutoff)
                and r.deam[0] >= 0
                else "nodeam"
            )
            if self.length_bin_size is None:
                LEN = 0
            else:
                LEN = r.len // self.length_bin_size

            DEAM = min((r.deam[0], r.deam[1]))

            if r.base == snp.ref:
                if r.base == "T" and snp.alt == "C" and not r.is_reverse:
                    dmgsite = 1
                elif r.base == "A" and snp.alt == "G" and r.is_reverse:
                    dmgsite = 1
                else:
                    dmgsite = 0
                D[r.RG, DEAM, LEN, dmgsite].n_ref += 1
            elif r.base == snp.alt:
                if r.base == "T" and snp.ref == "C" and not r.is_reverse:
                    dmgsite = 1
                elif r.base == "A" and snp.ref == "G" and r.is_reverse:
                    dmgsite = 1
                else:
                    dmgsite = 0
                D[r.RG, DEAM, LEN, dmgsite].n_alt += 1

        for (rg, deam, len_, dmgsite_), r in D.items():
            print(
                snp.chrom,
                snp.pos + 1,
                r.n_ref,
                r.n_alt,
                rg,
                len_,
                deam,
                dmgsite_,
                file=self.f,
                sep=",",
            )


class RefIter:
    def __init__(self, ref):
        self.ref = pd.read_csv(ref)
        self.bed = self.ref

    def __iter__(self):
        for ix, row in self.ref.iterrows():
            yield 0, (row.chrom, row.pos - 1, row.ref, row.alt)


def process_bam(
    outfile,
    bamfile,
    ref,
    deam_cutoff,
    length_bin_size,
    random_read_sample=False,
    max_reads=100,
    **kwargs,
):
    """generate input file from bam-file"""
    blocks = RefIter(ref)
    sampleset = pg.CallBackSampleSet.from_file_names([bamfile], blocks=blocks)

    default_filter.update(kwargs)
    logging.info("Filter is %s", default_filter)
    cov = AdmixfrogInput(
        **default_filter,
        random_read_sample=random_read_sample,
        length_bin_size=length_bin_size,
        deam_cutoff=deam_cutoff,
        outfile=outfile,
        max_reads=max_reads,
    )
    sampleset.add_callback(cov)
    sampleset.run_callbacks()


def process_bam2(
    outfile,
    bamfile,
    ref,
    deam_cutoff,
    length_bin_size,
    random_read_sample=False,
    **kwargs,
):
    """generate 2nd generation input file from bam-file
    file will have output cols
    chrom, pos, tref, talt, lib, len, deam, score
    """
    blocks = RefIter(ref)
    sampleset = pg.CallBackSampleSet.from_file_names([bamfile], blocks=blocks)

    default_filter.update(kwargs)
    logging.info("Filter is %s", default_filter)
    cov = AdmixfrogInput2(
        **default_filter,
        random_read_sample=random_read_sample,
        length_bin_size=length_bin_size,
        deam_cutoff=deam_cutoff,
        outfile=outfile,
    )
    sampleset.add_callback(cov)
    sampleset.run_callbacks()
