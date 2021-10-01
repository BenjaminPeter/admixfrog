import yaml

from abc import ABC, abstractmethod
from pybedtools import BedTool
from collections import Counter
import itertools

import pgdirect as pg

from math import ceil

"""
Brief overview:

We define a bunch of classes for a basic random read sampling api

The basic unit are objects of the File class. That reflect single sampled 
individuals and the associated files. Currently, 3 subclasses are implemented:
 - BamFile
 - VCFFile
 - FastaFile

 They represent samples in one of those formats. Files can have subsamples, which for now
 are defined by read-group and deamination levels.

 The most important function is gt(), an iterator over all SNPs in the sample

The most basic is the SampleSet, which collects all the samples, and allows
iterating over all files

The following hierarchy of samples is used
1. Population: a set of samples assumed to be exchangeable (e.g. CEU, YRI)
2. File: data from a single individual, will be procured from a single file
3. Subsample: a subset of data from an individual, e.g. technical filtering


Orthogonal, we have a set of FILES the data are stored in
Files have a one-to-many relationship with files, i.e. a file may contain
many samples, but a sample must come from a single file

thus, set up requires the following

-  defining populations
-  defining all samples
-  defining all files to be read
-  assigning files to samples
-  assigning samples to populations

before doing any computation

this is cumbersome, so the lazy set-up should be
-  assume each sample is a population, unless specified otherwise
-  assume all samples in all files are desired


"""
"""
Indexing note:

samtools pileup : 1 -indexed
pysam.pileup : 0 indexed

vcf : 1-indexed
pysam.vcf: 1-indexed !!!!

pysam.fasta: 0-indexed

"""


if __name__ == "__main__":

    with open("config/paths.yaml", "r") as f:
        y = yaml.load(f)

    def run_(C):
        global CHROM
        CHROM = C
        print(CHROM)
        bam_files = [
            y["paths"]["bam"]["VindijaG1"],
            y["paths"]["bam"]["Spy"],
            y["paths"]["bam"]["Mezmaiskaya2"],
            y["paths"]["bam"]["Les_Cottes"],
            y["paths"]["bam"]["Goyet"],
            y["paths"]["bam"]["Vindija"].format(CHROM=CHROM),
        ]

        bam_names = [
            "VindijaG1",
            "Spy",
            "Mezmaiskaya2",
            "Les_Cottes",
            "Goyet",
            "Vindija",
        ]
        bam_files += ["data/bam/salkhit_BigYoA_390k_840k_L35MQ25_deam53_merged.bam"]
        bam_names += ["Salkhit"]

        vcf_files = [
            y["paths"]["vcf"]["Vindija"],
            y["paths"]["vcf"]["Altai"],
            y["paths"]["vcf"]["Denisova"],
            y["paths"]["vcf"]["Mbuti"],
            y["paths"]["vcf"]["Khomani"],
        ]

        vcf_names = ["Vindija33.19", "Altai", "Denisova", "Mbuti", "Khomani"]

        fa_files = [y["paths"]["fa"]["Chimp"], y["paths"]["ref"]]

        fa_names = ["Chimp", "ref"]

        """now handiling this through fa_files"""
        # ref_file = y["paths"]["ref"]
        # ref = pysam.FastaFile(ref_file)
        # outgroup = pysam.FastaFile(fa_files[0].format(CHROM=CHROM))

        # file with snps to read
        bed_file = "/mnt/expressions/mateja/Early_modern_humans/Upper_Palaeolithic_project_Mateja/nuclear_captures/targets/2.2_million.snp.as.bed"
        bed = BedTool(bed_file)

        bams = [BamFile(bam_file, name) for bam_file, name in zip(bam_files, bam_names)]
        fastas = [
            FastaFile(fa.format(CHROM=CHROM), name)
            for fa, name in zip(fa_files, fa_names)
        ]
        # vcfs = [VCFFile(v.format(CHROM=CHROM), name) for v, name in zip(vcf_files, vcf_names)]

        snp_positions = [int(b[1]) for b in bed if b.chrom == CHROM]

        pos = (zip(itertools.cycle(b.chrom), range(b.start, b.end)) for b in bed)
        snp_positions = itertools.chain.from_iterable(pos)

        data = SampleSet(bams + vcfs + fastas, pos=snp_positions, ref=None)
        # pi = data.pi()
        Pi(data).write_csv("pi/pi_chr%s.txt" % CHROM)

        # snps = []
        # for pos in snp_positions:
        #    snps.append(SNP(data, None, CHROM, pos))

    CHROM = 1
    bam_files = [
        y["paths"]["bam"]["VindijaG1"],
        y["paths"]["bam"]["Spy"],
        y["paths"]["bam"]["Mezmaiskaya2"],
        y["paths"]["bam"]["Les_Cottes"],
        y["paths"]["bam"]["Goyet"],
        y["paths"]["bam"]["Vindija"],
    ]

    bam_names = ["VindijaG1", "Spy", "Mezmaiskaya2", "Les_Cottes", "Goyet", "Vindija"]
    bam_files += ["data/bam/salkhit_BigYoA_390k_840k_L35MQ25_deam53_merged.bam"]
    bam_names += ["Salkhit"]

    vcf_files = [
        y["paths"]["vcf"]["Vindija"],
        y["paths"]["vcf"]["Altai"],
        y["paths"]["vcf"]["Denisova"],
        y["paths"]["vcf"]["Mbuti"],
        y["paths"]["vcf"]["Khomani"],
    ]

    vcf_names = ["Vindija33.19", "Altai", "Denisova", "Mbuti", "Khomani"]

    fa_files = [y["paths"]["fa"]["Chimp"], y["paths"]["ref"]]

    fa_names = ["Chimp", "ref"]

    """now handiling this through fa_files"""
    # ref_file = y["paths"]["ref"]
    # ref = pysam.FastaFile(ref_file)
    # outgroup = pysam.FastaFile(fa_files[0].format(CHROM=CHROM))

    # file with snps to read
    bed_file = (
        "/mnt/expressions/mateja/Early_modern_humans/"
        + "Upper_Palaeolithic_project_Mateja/"
        + "nuclear_captures/targets/2.2_million.snp.as.bed"
    )
    bed = BedTool(bed_file)

    bams = [pg.BamFile(bam_file, name) for bam_file, name in zip(bam_files, bam_names)]
    fastas = [pg.FastaFile(fa, name) for fa, name in zip(fa_files, fa_names)]
    # vcfs = [VCFFile(v, name) for v, name in zip(vcf_files, vcf_names)]

    vcf_multi_test = [
        "/mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/"
        "merged_all_sites_arch_apes_sgdp1_g1000_chr1.vcf.gz"
    ]

    files = bam_files + fa_files + vcf_multi_test

    # data = SampleSet(bams +  vcfs + fastas, bed=bed, ref=None, max_snps=100, chrom='1')
    data = pg.SampleSet.from_file_names(
        files, bed=bed, ref=None, max_snps=100, chrom="1"
    )
    # pi = data.pi()
    # Pi(data).write_csv("pi/pi_chr%s.txt" % CHROM)

    def test_basic():
        with open("config/paths.yaml", "r") as f:
            y = yaml.load(f)
        bed_file = (
            "/mnt/expressions/mateja/Early_modern_humans/"
            + "Upper_Palaeolithic_project_Mateja/"
            + "nuclear_captures/targets/2.2_million.snp.as.bed"
        )
        NEA = ["Vindija33.19", "AltaiNeandertal", "Chagyrskaya-Phalanx"]
        bam_files = [
            y["paths"]["bam"]["VindijaG1"],
            y["paths"]["bam"]["Spy"],
            y["paths"]["bam"]["Mezmaiskaya2"],
            y["paths"]["bam"]["Les_Cottes"],
            y["paths"]["bam"]["Goyet"],
            y["paths"]["bam"]["Vindija"].format(CHROM=CHROM),
        ]
        vcf_multi_test = [
            "/mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/"
            "merged_all_sites_arch_apes_sgdp1_g1000_chr1.vcf.gz"
        ]

        p1 = [
            pg.Sample.from_file(b, sampling=pg.SAMPLING.RANDOM_RG)
            for b in bam_files[:3]
        ]
        # p1 = [ss for sample in p0 for ss in sample.subsamples]
        vcf_multi = pg.VCFFile(vcf_multi_test[0])
        p2 = [s for s in vcf_multi.samples if s.name in NEA]
        p3 = {s for s in vcf_multi.samples if "Mbuti" in s.name}
        p4 = {pg.Sample.from_file(f) for f in fa_files[:1]}
        fa_files = [y["paths"]["fa"]["Chimp"], y["paths"]["ref"]]

        pops = [Population(p1), Population(p2), p3, p4]
        blocks = NoBlocks(bed, max_snps=1005)
        self = Dstat(pops, blocks)
        block, snp = next(self.sampleset.gt())
        block, snp = next(self.sampleset.gt())
        quartet = next(itertools.product(*self.popssub))
        [(snp[s], s.name, s.sampling) for s in self.subsamples]
