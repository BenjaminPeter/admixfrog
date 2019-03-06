import admixfrog.pgdirect as pg
import pandas as pd
from collections import defaultdict
import lzma
import logging

default_filter = {
    "deam_only": False,
    "pos_in_read_cutoff": 2,
    "min_length": 35,
    "max_length": 1000,
    "minq": 25,
}


class AdmixfrogInput(pg.ExtCoverage):
    def __init__(self, outfile, deam_cutoff, length_bin_size=None, **kwargs):
        self.outfile = outfile
        self.deam_cutoff = deam_cutoff
        self.length_bin_size = length_bin_size
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
                lib = "{rg}_{deam}".format(rg=rg, deam=deam)
            else:
                # lib = f"{rg}_{len_}_{deam}"
                lib = "{rg}_{len_}_{deam}".format(rg=rg, len_=len_, deam=deam)
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


class RefIter:
    def __init__(self, ref):
        self.ref = pd.read_csv(ref)
        self.bed = self.ref

    def __iter__(self):
        for ix, row in self.ref.iterrows():
            yield 0, (row.chrom, row.pos - 1, row.ref, row.alt)


def process_bam(outfile, bamfile, ref, deam_cutoff, length_bin_size, **kwargs):
    blocks = RefIter(ref)
    sampleset = pg.CallBackSampleSet.from_file_names([bamfile], blocks=blocks)

    default_filter.update(kwargs)
    logging.info("Filter is %s", default_filter)
    cov = AdmixfrogInput(
        **default_filter,
        length_bin_size=length_bin_size,
        deam_cutoff=deam_cutoff,
        outfile=outfile
    )
    sampleset.add_callback(cov)
    sampleset.run_callbacks()
