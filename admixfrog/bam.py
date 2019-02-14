import admixfrog.pgdirect as pg
from collections import defaultdict
import lzma

default_filter = {
    "deam_only" : False,
    "pos_in_read_cutoff" : 0,
    "min_length" : 35,
    "max_length" : 1000,
    "minq" : 25
}

class AdmixfrogInput(pg.ExtCoverage):
    def preprocess(self, sampleset):
        self.f = lzma.open(self.outfile, "wt")
        print("chrom", "pos", "lib", "tref", "talt", "tdeam", "tother",
              sep = ",", file=self.f)

    def process_snp(self, block, snp):
        reads = snp.reads(**self.kwargs)
        D = defaultdict(lambda : self.Obs())
        #n_ref, n_alt, n_deam, n_other = 0, 0, 0, 0
        for r in reads:
            DEAM = 'deam' if (r.deam[0] < 3 or r.deam[1] < 3) and r.deam[0] >=0 else 'nodeam'
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
            print(snp.chrom, snp.pos + 1, f"{rg}_{deam}", r.n_ref, r.n_alt, r.n_deam, r.n_other,
              file=self.f, sep =",")

def process_bam(outfile, bamfile, bedfile, deam_cutoff, length_bin_size, **kwargs):
    print(kwargs)
    blocks = pg.NoBlocks(bedfile)
    sampleset = pg.CallBackSampleSet.from_file_names([bamfile], blocks=blocks)

    default_filter.update(kwargs)
    print(default_filter)
    cov = AdmixfrogInput(**default_filter, deam_cutoff=deam_cutoff, outfile=outfile)
    sampleset.add_callback(cov)
    sampleset.run_callbacks()
