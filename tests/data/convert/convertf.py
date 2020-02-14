import os
CHROMS = ['9', '15', 'X']
for CHROM in CHROMS:
    s = f"""genotypename:    /mnt/sequencedb/gendivdata/2_genotypes/reich_combined/v37.2.1240K.geno
snpname:         /mnt/sequencedb/gendivdata/2_genotypes/reich_combined/v37.2.1240K.snp
indivname:       /mnt/sequencedb/gendivdata/2_genotypes/reich_combined/v37.2.1240K.ind
outputformat:    PACKEDANCESTRYMAP
genotypeoutname: test_oase{CHROM}.geno
snpoutname:      test_oase{CHROM}.snp
indivoutname:    test_oase{CHROM}.ind
poplistname:     oase.pops
minchrom:		 {CHROM}
maxchrom:		 {CHROM}
"""
    print(s, file=open('convert/convert1.par', 'w'))
    os.system("convertf -p convert/convert1.par")


C1 = ['9', 'AUTO']
C2 = ['15', 'X']
C3 = ['AUTO', 'ALL']
for CHROM, CHROM2, CHROM12 in zip(C1, C2, C3):
    s= f'''
geno1: test_oase{CHROM}.geno
snp1:  test_oase{CHROM}.snp
ind1:  test_oase{CHROM}.ind
geno2: test_oase{CHROM2}.geno
snp2:  test_oase{CHROM2}.snp
ind2:  test_oase{CHROM2}.ind
genotypeoutname: test_oase{CHROM12}.geno
snpoutname:      test_oase{CHROM12}.snp
indivoutname:    test_oase{CHROM12}.ind
'''
    print(s, file=open('convert/merge1.par', 'w'))
    os.system("mergeit -p convert/convert1.par")

