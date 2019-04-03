import pdb
import lzma
from pysam import VariantFile
from collections import defaultdict
import pandas as pd
import yaml
from pprint import pprint, pformat

EXT = "ref", "alt"

x = defaultdict(lambda s: s)
x['AFR'] = 'S_Mbuti-1', 'S_Mbuti-2', 'S_Mbuti-3'
x['PAN'] = 'panTro4',
x['OCE'] = 'S_Papuan-1', 'S_Papuan-2', 'S_Papuan-3', 'S_Papuan-4' 
pop2sample = x

rec_file = '/home/benjamin_peter/proj/100a/basic_processing/recs/maps_b37/maps_chr.{CHROM}'
vcf_file = "/mnt/expressions/bpeter/100a/basic_processing/vcfs/merged/tinyvcf_afr.vcf.gz"

def load_pop_file(pop_file=None, pops=None):
    P = dict()
    P2 = dict()
    if pop_file is not None:
        with open(pop_file) as f:
            y = yaml.load(f, Loader=yaml.FullLoader)
            y = y['sampleset'] if 'sampleset' in y else y
            P.update(y)
    if pops is None:
        return P
    else:
        for p in pops:
            if p in P:
                P2[p] = P[p]
            else:
                P2[p] = [p]
        return P2

def load_random_read_samples(pop_file=None):
    if pop_file is None:
        return []
    with open(pop_file) as f:
        y = yaml.load(f, Loader=yaml.FullLoader)
        y = y['random_read_samples'] if 'random_read_samples' in y else y
        return y


def vcf_to_ref(outfile, vcf_file, rec_file, 
               pop2sample, random_read_samples=[],
               pos_id='Physical_Pos', map_id='AA_Map',
               rec_rate = 1e-8
               ):

    data_cols = [f"{p}_{e}"  for p in pop2sample.keys() for e in EXT]

    #  get chromosomes
    with  VariantFile(vcf_file) as vcf:
        chroms = [i for i in vcf.header.contigs]

        sample2pop = defaultdict(list)
        for pop, v in pop2sample.items():
            for sample in v:
                if sample in vcf.header.samples:
                    sample2pop[sample].append(pop)

    samples = sample2pop.keys()
    pprint(sample2pop)

    with lzma.open(outfile, 'wt') as ref:
        ref.write("chrom,pos,ref,alt,map,")
        ref.write(",".join(data_cols))
        ref.write("\n")
        for chrom in chroms:

            # set up rec file
            if rec_file is not None:
                rec = pd.read_csv(rec_file.format(CHROM=chrom), sep=' ')
                if 'chrom' in rec:
                    rec = rec[rec.chrom==chrom]
                rec = rec[[pos_id, map_id]]
                rec_iter = rec.iterrows()
                R0 = next(rec_iter)[1]
                R1 = next(rec_iter)[1]
            
            with  VariantFile(vcf_file) as vcf:
                vcf.subset_samples(samples)
                for row in vcf.fetch(chrom):
                    if len(row.alleles) != 2: 
                        continue
                    D = defaultdict(int)
                    # rec stuff
                    if rec_file is None:
                        map_ = row.pos * rec_rate
                    else:
                        if R1 is None:
                            map_ = R0[map_id]
                        elif  row.pos <= R0[pos_id]:
                            map_ = R0[map_id]
                        elif  R0[pos_id] < row.pos <= R1[pos_id]:
                            slope = (R1[map_id] - R0[map_id])/(R1[pos_id] - R0[pos_id])
                            map_ = R0[map_id] + slope * (row.pos - R0[pos_id])/ (R1[pos_id] - R0[pos_id])
                        elif row.pos > R1[pos_id]:
                            try:
                                while row.pos > R1[pos_id]:
                                    R0, R1 = R1, next(rec_iter)[1]
                            except StopIteration:
                                R0, R1 = R1, None
                            if R1 is None:
                                map_ = R0[map_id]
                            else:
                                slope = (R1[map_id] - R0[map_id])/(R1[pos_id] - R0[pos_id])
                                map_ = R0[map_id] + slope * (row.pos - R0[pos_id])/ (R1[pos_id] - R0[pos_id])

                    ref.write(f'{row.chrom},{row.pos-1},{row.ref},{row.alts[0]},{map_},')

                    sample_data = row.samples
                    for s in sample_data:
                        if s in random_read_samples:
                            allele = sample_data[s]['GT'][0]
                            if allele is not None:
                                for pop in sample2pop[s]:
                                    D[f'{pop}_{EXT[allele]}'] += 1
                        else:
                            for allele in sample_data[s]['GT']:
                                if allele is not None:
                                    for pop in sample2pop[s]:
                                        D[f'{pop}_{EXT[allele]}'] += 1
        
                    ref.write(",".join((str(D[c]) for c in data_cols)))
                    ref.write("\n")
