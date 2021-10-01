from .sampleset import *
from collections import defaultdict, Counter
from .callback import Callback
from .utils import NA
import itertools


class Stat(Callback):
    """Stat
    generic Statistic

    the intended use-case is the following. We would like to calculate a 
    statistic over multiple genomes, from a variety of source files.
    
    The following features are desirable:
    1. handle plink, vcf, bam and fasta files
    1. group multiple individuals into populations
    3. group multiple populations into equivalent ``slots'', calculating
        statistics for each slot
    4. restrict to a subset of genome
    4. use (blocked) resampling to estimate error
    5. handle ancestral/derived state
    6. allow different seq qualities, including on-the-fly read sampling
    7. subset samples according to read group, deamination level
    8. streaming processing

    what we want to run is ideally

    as an use-case, suppose we want to run the following analysis:

    1. analyze only sites where we know Neandertals are polymorphic
    2. use blocks of 5000 SNP each for resampling
    3. split bam-files such that each read group is analyzed separately
    3. analyze terminally deaminated and non-deaminated fragments separately
    4. Merge modern humans from the same population
    5. treat high-coverage neandertals individually

    bed = Bedfile("nea_sites.bed")  # read bed with desired regions
    blocks = BlockSNP(bed, block_size = 5000)  #make blocks of 5000 SNP each

    giantVCF = VCFFile("giantvcf_chr{CHROM}.gz")
    Mbuti = Population(s for s in giantVCF.samples if "Mbuti" in s.name)
    Khomani = Population(s for s in giantVCF.samples if "Khomani" in s.name)
    French = Population(s for s in giantVCF.samples if "French" in s.name)

    slot1 = [Sample.from_file(bam, sampling=SAMPLING.RANDOM_RGDEAM) for bam in bam_files]
    slot2 = [s for s in giantVCF.samples if s.name in NEA]
    slot3 = [Mbuti, Khomani, French]
    slot4 = [Sample.from_file(fa) for fa in fasta_files]

    anc = [Chimp_chr{CHROM}.fa]

    DStat([slot1, slot2, slot3, slot4], blocks=blocks) #all D-stats
    F3Stat([slot1, slot2, slot3], blocks=blocks) #all F3-stats
    PiStat(slot1 + slot2 + slot3 + slot4, blocks=blocks) #pw differences matrix
    NDerived(slot1 + slot2 + slot3, anc=anc) #get avg number of derived alleles


    """

    def __init__(self, slots, anc=None, blocks=None, *args, **kwargs):
        self.slots = slots  # [Slot1, Slot2, Slot3, Slot4]
        self.analysis_units = None
        """
        [slot1.Pop1, slot1.Pop2, slot1.Sample3, slot2.Sample1, 
        slot3.sample1.Subsample1, slot3.sample1.subsample2, slot4.pop1]
        """
        self.samples = None
        """[slot1.pop1.sample1, slot1.pop2.sample2, slot1.pop2.sample1, ...
        slot3.sample1, slot4.pop1.sample1
        ]
        """


class Dstat(Callback):
    def __init__(self, slots, outname="d.test", verbose=False, **kwargs):
        self.slots = slots
        for slot in self.slots:
            assert len(slot) > 0
        self.verbose = verbose
        self.outname = outname

    def preprocess(self, sampleset):
        # dicts of form ABBA[block][(s1, s2, s3, s4)]
        self.patterns = Counter()  # key is (quartet, block, pattern)

    def postprocess(self, sampleset):
        with open(self.outname, "w") as f:
            line = ["S1", "S2", "S3", "S4"]
            line += ["block", "pattern" "n"]
            print(*line, sep=",", file=f)

            for (quartet, block, pattern), n in self.patterns.most_common():
                print(",".join(q.name for q in quartet), file=f, end=",")
                print(block, pattern, n, sep=",", file=f)

    def process_snp(self, block, snp):
        for quartet in itertools.product(*self.slots):
            n = len(list(itertools.product(*(snp[ind].gt() for ind in
                                             quartet))))

            pattern = itertools.product(*(snp[ind].gt() for ind in quartet))
            for p in pattern:
                # p = "".join(p)
                if NA not in p:
                    # print(p)
                    # if self.is_abba(p) or self.is_baba(p) or self.is_bbaa(p):
                    #    self.patterns[quartet, block, p] += 1
                    if self.is_abba(p):
                        self.patterns[quartet, block, "ABBA"] += 1. / n
                    elif self.is_baba(p):
                        self.patterns[quartet, block, "BABA"] += 1. / n
                    elif self.is_bbaa(p):
                        self.patterns[quartet, block, "BBAA"] += 1. / n

    def has_data(snp, q):
        return all(snp[i] != NA for i in q)

    @staticmethod
    def is_abba(p):
        return p[0] == p[3] and p[1] == p[2] and p[0] != p[1]

    @staticmethod
    def is_baba(p):
        return p[0] == p[2] and p[1] == p[3] and p[0] != p[1]

    @staticmethod
    def is_bbaa(p):
        return p[0] == p[1] and p[2] == p[3] and p[0] != p[2]

class Dstatl(Dstat):
    def postprocess(self, sampleset):
        with open(self.outname, "w") as f:
            line = ["S1", "S2", "S3", "S4"]
            line += ["L3"]
            line += ["d5", "d3", "n_deam"]
            line += ["block", "pattern", "n"]
            print(*line, sep=",", file=f)

            for (quartet, block, length, deam, pattern), n in self.patterns.most_common():
                #print(*quartet, *length,  *deam, block, pattern, n, sep=",")
                print(*quartet, length,  *deam, block, pattern, n, sep=",", file=f)

    def process_snp(self, block, snp):
        for quartet in itertools.product(*self.slots):
            n = 1 #len(list(itertools.product(*(snp[ind].gt() for ind in
                 #                            quartet))))

            #tmp_list = []
            #for i, ind in enumerate(quartet):
            #    gt = snp[ind].gt(report_length=True)
            #    tmp_list.append(gt)
            #    if i == 2:
            #        print(snp[ind].reads())


            tmp_list = [snp[ind].gt(return_read=True,
                                    pos_in_read_cutoff=0,
                                    min_length=0) for ind in quartet]


            try:
                p = tuple((gt.base for gt in tmp_list  ))
            except Exception as e:
                print(tmp_list)
                raise
            length = tuple((gt.len for gt in tmp_list  ))
            deam = tuple(((*gt.deam, len(gt.deam_full)) for gt in
                      tmp_list  ))

            #for p, length, deam in zip(pattern, lengths, other):
                # p = "".join(p)
            if NA not in p:
                # if self.is_abba(p) or self.is_baba(p) or self.is_bbaa(p):
                #    self.patterns[quartet, block, p] += 1
                qn = tuple(q.name for q in quartet)
                if self.is_abba(p):
                    self.patterns[qn, block, length[2], deam[2], "ABBA"] += 1. / n
                elif self.is_baba(p):
                    self.patterns[qn, block, length[2], deam[2], "BABA"] += 1. / n
                elif self.is_bbaa(p):
                    self.patterns[qn, block, length[2], deam[2], "BBAA"] += 1. / n
                else :
                    self.patterns[qn, block, length[2], deam[2], "other"] += 1. / n


__all__ = ["Stat", "Dstat", "Dstatl"]
