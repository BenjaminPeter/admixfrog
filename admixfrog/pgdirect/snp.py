from collections import Counter
from functools import lru_cache as cache
from .utils import NA
from collections import defaultdict
from pysam import VariantRecord
from .sample import Population

import pdb


VCF_ROW = "{chrom}\t{pos}\t{chrom}_{pos}\t{ref}\t{alt}\t.\t.\t.\tGT"
VCF_ROWDP = "{chrom}\t{pos}\t{chrom}_{pos}\t{ref}\t{alt}\t.\t.\t.\tGT:DP"


class SNP(object):
    def __init__(self, coords, data, anc=None):
        """__init__
            Parameters
            ----------
            dataset : [removed] SampleSet | maybe just pass Structure object
                Sample set that stores all possible individual samples
            coords : str , int, [ref] [alt]
                Chromosome and position of SNP
            anc : :Sample:
                the sample reflecting the ancestral allele
            -------
        """
        self.has_ref = len(coords) == 4
        self.coords = coords
        self.chrom, self.pos = coords[:2]

        self.ref = coords[2] if self.has_ref else None
        self.alt = coords[3] if self.has_ref else None

        self._data = dict()

        for gtset in data:
            for sample, gt in gtset:
                self._data[sample] = gt

        self.anc_sample = anc if anc in self._data else None
        self._freq, self._freq2 = None, None

    @property
    @cache(maxsize=1)
    def freq_tot(self):
        return  Counter(self._data[s] for s in self._data)

    @property
    @cache(maxsize=1)
    def freq(self):
        return Counter(
                self._data[s].gt() for s in self._data if self._data[s].gt() != NA
            )

    @property
    def n_missing(self):
        return self.freq_tot[NA]

    @property
    def n_called(self):
        return sum(self.freq.values())

    @property
    @cache(maxsize=1)
    def f(self):
        return dict((al, cnt / self.n_called) for al, cnt in self.freq.items())

    @property
    def maf(self):
        """minor allele frequency"""
        if self.is_monoallelic:
            return 0
        elif self.is_biallelic:
            return min(self.f.values())
        else:
            raise ValueError("maf undefined for more than two alleles")

    @property
    def minor_allele(self):
        f = self.freq
        if len(f) > 1:
            return f.most_common(2)[1][0]
        return None

    @property
    def major_allele(self):
        f = self.freq
        if len(f) > 1:
            return f.most_common(1)[1][0]
        return None

    def __iter__(self):
        """iterate over (sample, gt)
        """
        return iter(self._data.items())

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            if type(item) is Population:
                return [self[s] for s in item.samples]

    def reads(self, **kwargs):
        reads = []
        for k, v in self._data.items():
            r = v.reads(**kwargs)
            reads.extend(r)

        return reads
            

    def to_vcfrecord(self, samples, depth=True, **kwargs):
        # VCF_ROW = '{chrom}\t{pos}\t{chrom}_{pos}\t{ref}\t{alt}\t.\t.\t.\tGT:DP\t'

        has_data = False

        # if self.ref is None or self.alt is None:
        #    ref, alt = self.major_allele, self.minor_allele
        #    if ref is None:
        #        ref = "."
        #    if alt is None:
        #        alt = "."
        # else:

        ref, alt = self.ref, self.alt
        if depth:
            row = VCF_ROWDP.format(chrom=self.chrom, pos=self.pos + 1, ref=ref, alt=alt)
        else:
            row = VCF_ROW.format(chrom=self.chrom, pos=self.pos + 1, ref=ref, alt=alt)

        for s in samples:
            gt = self._data[s].vcfgt(
                ref, alt, depth=depth, coords=self.coords, **kwargs
            )
            if not gt.startswith("."):
                has_data = True
            row = row + "\t" + gt

        if has_data:
            return row + "\n"
        else:
            return ""



__all__ = ["SNP"]
