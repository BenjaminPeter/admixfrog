import logging
from collections import Counter
import numpy as np
import lzma
from pysam import VariantFile
from collections import defaultdict
import pandas as pd
import yaml
from pprint import pprint, pformat

EXT = "ref", "alt"


def parse_chroms(arg):
    if arg is None: return None
    chroms = []
    for s in arg.split(","):
        if "-" in s:
            a, b = s.split("-")
            chroms.extend([str(s) for s in range(int(a), int(b) + 1)])
        else:
            chroms.append(s)
    return chroms


def load_pop_file(pop_file=None, pops=None):
    P = dict()
    P2 = dict()
    if pop_file is not None:
        with open(pop_file) as f:
            y = yaml.load(f, Loader=yaml.FullLoader)
            y = y["sampleset"] if "sampleset" in y else y
            P.update(y)
    print(pops)
    if pops is None:
        return P
    else:
        for p in pops:
            if p in P:
                P2[p] = P[p]
            elif "=" in p:
                k, v = p.split("=")
                P2[k] = v.split(",")
                tmp_samples = []
                for set_ in P2[k]:
                    if set_ in P:
                        tmp_samples.extend(P[set_])
                    else:
                        tmp_samples.append(set_)
                P2[k] = tmp_samples
            else:
                P2[p] = [p]
        return P2


def vcf_to_ref(
    outfile,
    vcf_file,
    rec_file,
    pop2sample,
    random_read_samples=[],
    pos_id="Physical_Pos",
    map_ids=["AA_Map"],
    default_map="AA_Map",
    rec_rate=1e-8,
    chroms=None,
    bed=None,
    lax_alleles=False,
):

    pprint(pop2sample)

    if chroms is None:
        vcf_first = vcf_file
    else:
        vcf_first = vcf_file.format(CHROM=chroms[0])

    #  get chromosomes
    with VariantFile(vcf_first) as vcf:
        if chroms is None:
            chroms = [i for i in vcf.header.contigs]
        else:
            chroms = parse_chroms(chroms)
        logging.info("chroms found: %s", chroms)

        sample2pop = defaultdict(list)
        for pop, v in pop2sample.items():
            for sample in v:
                if sample in vcf.header.samples:
                    sample2pop[sample].append(pop)

    samples = sample2pop.keys()
    pops = sorted(set(pop for s, v in sample2pop.items() for pop in v))
    pprint(sample2pop)
    pprint(pops)

    map_ids = ["map"] + map_ids

    data_cols = [f"{p}_{e}" for p in pops for e in EXT]

    with lzma.open(outfile, "wt") as ref:
        ref.write("chrom,pos,ref,alt,")
        if rec_file is None:
            ref.write("map,")
        else:
            ref.write(",".join(map_ids))
            ref.write(",")
        ref.write(",".join(data_cols))
        ref.write("\n")
        for chrom in chroms:

            # set up rec file
            if rec_file is not None:
                rec = pd.read_csv(rec_file.format(CHROM=chrom), sep=" ")
                if "chrom" in rec:
                    rec = rec[rec.chrom == chrom]

                rec["map"] = rec[default_map]
                rec_file_cols = list((pos_id, *map_ids))
                rec = rec[rec_file_cols]

                rec_iter = rec.iterrows()
                R0 = next(rec_iter)[1]
                R1 = next(rec_iter)[1]

            # skip chrom if empty
            with VariantFile(vcf_file.format(CHROM=chrom)) as vcf:
                try:
                    V = next(vcf)
                except StopIteration:
                    continue
            with VariantFile(vcf_file.format(CHROM=chrom)) as vcf:
                vcf.subset_samples(samples)
                for row in vcf.fetch(chrom):

                    alt_ix = 0

                    if len(row.alleles) <= 1 or len(row.alleles) > 3:
                        continue

                    if len(row.alleles) == 3:
                        alleles = [i for v in row.samples.values() for i in v["GT"]]
                        if 3 in alleles:
                            continue
                        elif 1 in alleles and 2 in alleles:
                            continue
                        elif 1 not in alleles and 2 not in alleles:
                            continue
                        elif 1 in alleles:
                            alt_ix = 0
                        elif 2 in alleles:
                            alt_ix = 1
                        else:
                            raise ValueError(f"weird alleles {row.alleles}")
                        logging.debug(
                            f"{row.chrom}, {row.pos}, {row.alleles}, {Counter(alleles)}"
                        )

                    if row.alts[alt_ix] not in "ACGT" or lax_alleles:
                        continue

                    D = defaultdict(int)
                    # rec stuff
                    if rec_file is None:
                        map_ = row.pos * rec_rate
                        ref.write(
                            f"{row.chrom},{row.pos},{row.ref},{row.alts[alt_ix]},{map_},"
                        )
                    else:
                        if R1 is None:
                            map_ = R0[map_ids]
                        elif row.pos <= R0[pos_id]:
                            map_ = R0[map_ids]
                        elif R0[pos_id] < row.pos <= R1[pos_id]:
                            slope = (R1[map_ids] - R0[map_ids]) / (
                                R1[pos_id] - R0[pos_id]
                            )
                            map_ = R0[map_ids] + slope * (row.pos - R0[pos_id]) / (
                                R1[pos_id] - R0[pos_id]
                            )
                        elif row.pos > R1[pos_id]:
                            try:
                                while row.pos > R1[pos_id]:
                                    R0, R1 = R1, next(rec_iter)[1]
                            except StopIteration:
                                R0, R1 = R1, None
                            if R1 is None:
                                map_ = R0[map_ids]
                            else:
                                slope = (R1[map_ids] - R0[map_ids]) / (
                                    R1[pos_id] - R0[pos_id]
                                )
                                map_ = R0[map_ids] + slope * (row.pos - R0[pos_id]) / (
                                    R1[pos_id] - R0[pos_id]
                                )

                        ref.write(
                            f"{row.chrom},{row.pos},{row.ref},{row.alts[alt_ix]},"
                        )
                        map_str = ",".join((str(m) for m in map_))
                        ref.write(f"{map_str},")

                    sample_data = row.samples
                    for s in sample_data:
                        if s in random_read_samples:
                            allele = sample_data[s]["GT"][0]
                            if allele is not None:
                                for pop in sample2pop[s]:
                                    D[f"{pop}_{EXT[allele > 0]}"] += 1
                        else:
                            for allele in sample_data[s]["GT"]:
                                if allele is not None:
                                    for pop in sample2pop[s]:
                                        D[f"{pop}_{EXT[allele > 0]}"] += 1

                    ref.write(",".join((str(D[c]) for c in data_cols)))
                    ref.write("\n")


def vcf_to_sample(
    outfile, vcf_file, ref_file, sample_id, random_read=False, chroms=None
):
    if chroms is None:
        with VariantFile(vcf_file.format(CHROM="1")) as vcf:
            chroms = [i for i in vcf.header.contigs]
    else:
        chroms = parse_chroms(chroms)

    logging.debug("chroms found: %s", chroms)

    ref = pd.read_csv(ref_file)
    ref.chrom = ref.chrom.astype(str)

    with lzma.open(outfile, "wt") as infile:
        infile.write(f"chrom,pos,tref,talt\n")

        for chrom in chroms:

            ref_local = ref[ref.chrom == chrom]
            with VariantFile(vcf_file.format(CHROM=chrom)) as vcf:
                vcf.subset_samples([sample_id])
                for row in vcf.fetch(chrom):
                    ALT_INDEX = 1

                    if row.pos in ref_local.pos.values:
                        if len(row.alleles) == 1:
                            ALT_INDEX = -1
                        elif len(row.alleles) != 2:
                            ref_row = ref_local[ref_local.pos == row.pos]
                            ref_alt = ref_row.alt.values
                            ALT_INDEX = np.where(ref_alt == row.alleles)[0].item()
                        # logging.debug(f"{row.chrom}:{row.pos} in ref")
                    else:
                        # logging.debug(f"{row.chrom}:{row.pos} not in ref, skipping")
                        continue

                    infile.write(f"{row.chrom},{row.pos},")

                    sample_data = row.samples
                    for s in sample_data:
                        if random_read:
                            allele = sample_data[s]["GT"][0]
                            if allele == 0:
                                infile.write("1,0\n")
                            elif allele == ALT_INDEX:
                                infile.write("0,1\n")
                            else:
                                infile.write("0,0\n")
                        else:
                            alleles = Counter(sample_data[s]["GT"])
                            # breakpoint()
                            infile.write(f"{alleles[0]},{alleles[ALT_INDEX]}\n")
            logging.debug(f"done processing {chrom}")


def load_random_read_samples(pop_file=None, random_read_samples=[]):
    pseudo_haps = []
    if pop_file is not None:
        with open(pop_file) as f:
            y = yaml.load(f, Loader=yaml.FullLoader)
            if "pseudo_haploid" in y:
                pseudo_haps.extend(y["pseudo_haploid"])
    pseudo_haps.extend(random_read_samples)
    print(pseudo_haps)
    return pseudo_haps
