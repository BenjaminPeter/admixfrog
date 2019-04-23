import random
from functools import lru_cache
from .utils import unambiguous, NA, NA_READ, Read
from .utils import NA as NA_
from abc import ABC, abstractmethod
import pdb

VCF_GTDP = "{gt}:{dp}"
VCF_GT = "{gt}"


class GT(ABC):
    def __init__(self, gt, sample=None, chrom=None, pos=None):
        self.chrom = chrom
        self._gt = gt
        self.pos = pos
        self.sample = sample

        # make properties
        #        self.sampling = sample.sampling if sample is not None else None
        self.name = None if sample is None else sample.name

    def __str__(self):
        tpl = (self.name, self.chrom, self.pos, "".join(self.gt()))
        return "%s[%s:%s]: %s" % tpl

    def __repr__(self):
        tpl = (self.name, self.chrom, self.pos, "".join(self.gt()))
        return "%s[%s:%s]: %s" % tpl

    @property
    def alleles(self):
        s = set(self.gt())
        s.discard(NA)
        return s

    @abstractmethod
    def gt():
        """
        the main method of obtianing a genotype; will return either
        ACGT/ 0,1,2, NA depending on what sampling scheme is asked
        """
        pass

    @abstractmethod
    def vcfgt(self, depth=False, **kwargs):
        """
        get gt in vcf format
        """
        pass


class GTNoData(GT):
    def __init__(self, *args, **kwargs):
        super().__init__(NA, *args, **kwargs)

    def gt(self, report_depth=False, return_read=False, **kwargs):
        NA = NA_READ if return_read else NA_
        if report_depth:
            return NA, 0
        return NA

    def vcfgt(self, ref=".", alt=".", depth=True, *args, **kwargs):
        if depth:
            return "./.:0"
        else:
            return "./."

    def reads(self, **kwargs):
        return []


class DiploidGT(GT):
    def __init__(self, gt, sample, chrom=None, pos=None, record=None):
        super().__init__(gt, sample, chrom, pos)
        if gt is None or gt[0] is None:
            self._gt = NA
        else:
            self._gt = "".join(gt)

        # self.record = None if record is None else record

    @lru_cache(maxsize=64)
    def gt(self,  return_read=False,**kwargs):
        if return_read:
            return Read(base=random.choice(self._gt))
        return self._gt

    def vcfgt(self, ref, alt, depth=True, **kwargs):
        gt = self.gt()
        if depth:
            if gt == NA:
                return "./.:0"
            a1, a2 = gt
            if a1 != ref and a1 != alt:
                return "./.:0"
            if a2 != ref and a2 != alt:
                return "./.:0"

            GT = "%d/%d:-1" % (a1 != ref, a2 != ref)
            return GT
        else:
            if gt == NA:
                return "./."
            a1, a2 = gt
            if a1 != ref and a1 != alt:
                return "./."
            if a2 != ref and a2 != alt:
                return "./."

            GT = "%d/%d:-1" % (a1 != ref, a2 != ref)
            return GT


class HaploidGT(GT):
    def __init__(self, gt, sample, chrom=None, pos=None):
        super().__init__(gt, sample, chrom, pos)
        self._gt = unambiguous[self._gt]

    def gt(self, return_read=False,**kwargs):
        if return_read:
            return Read(base=self._gt)
        return self._gt

    def vcfgt(self, ref, alt, depth=True, **kwargs):
        gt = self.gt()
        if gt == NA:
            return ".:0"
        if depth:
            if gt != ref and gt != alt:
                return ".:0"
            GT = "%d:1" % (gt != ref)
            return GT
        else:
            if gt != ref and gt != alt:
                return "."
            GT = "%d" % (gt != ref)
            return GT


class ReadGT(GT):
    # def __init__(self, gt, *args, **kwargs):
    #    super().__init__(self, *args, **kwargs)

    def __str__(self):
        tpl = (self.name, self.chrom, self.pos)
        return "%s[%s:%s]:" % tpl

    @staticmethod
    def qc(
        read,
        deam_only=False,
        deam_bases=3,
        pos_in_read_cutoff=3,
        min_length=30,
        max_length=1000,
        minq=25,
        minmapq=25
           ):
        if read.len < min_length or read.len > max_length:
            return False
        if read.bq < minq:
            return False
        if read.mq < minmapq:
            return False
        if read.pos_in_read < pos_in_read_cutoff:
            return False
        if read.len - read.pos_in_read < pos_in_read_cutoff:
            return False
        if deam_only:
            deam5, deam3 = read.deam
            if deam5 < 0 and deam3 < 0:
                return False
            if deam5 > deam_bases and deam3 > deam_bases:
                return False
        return True

    def __repr__(self):
        return self.__str__()

    def reads(self, **kwargs):
        return [r for r in self._gt if self.qc(r, **kwargs)]

    @lru_cache(maxsize=64)
    def gt(
        self,
        max_depth=400,
        report_depth=False,
        return_read=False,
        **kwargs
    ):
        """assume random sampling for now
        self._gt is a list of Reads
        """
        #print("READGT:: report_length", report_length)

        NA = NA_READ if return_read else NA_
        r_iter = self.reads(**kwargs)

        depth = len(r_iter)
        if depth > max_depth:
            if report_depth:
                return NA, depth
            return NA

        if report_depth:
            if depth == 0:
                return NA, 0
            return random.choice(r_iter), depth
        elif return_read:
            try:
                rand_read = random.choice(r_iter)
                return rand_read  #.base, rand_read.len, (*rand_read.deam), rand_read.pos_in_read
            except IndexError:
                return NA #-1, None, None, None
        else:
            try:
                return random.choice(r_iter).base
            except IndexError:
                return NA

    def vcfgt(self, ref, alt, depth=True, coords=None,  **kwargs):
        if depth:
            gt, d = self.gt(report_depth=depth,  **kwargs)
        else:
            gt = self.gt(report_depth=depth, **kwargs)

            
        if depth:
            if gt == NA:
                return ".:0"
            if gt.base != ref and gt.base != alt:
                print("third allele", gt, ref, alt, coords)
                return ".:0"

            GT = "%d:%d" % (gt.base != ref, d)
            return GT
        else:
            if gt == NA:
                return "."
            if gt.base != ref and gt.base != alt:
                return "."

            GT = "%d" % (gt.base != ref)
            return GT



__all__ = ["GT", "GTNoData", "DiploidGT", "HaploidGT", "ReadGT"]
