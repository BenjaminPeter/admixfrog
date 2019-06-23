some simple tests. For now, they mostly cover correctly reading/writing files

we have the following formats we need to consider / test

1. reference
    - native (csv)
    - vcf
    - geno
    - some genotype likelihood (?) (NYI, but important)

2. target
    - native (csv)
    - bam
    - geno
        - random read
        - genotype
    - vcf
        - random read
        - genotype

3. assignment of individuals to pops
    - str
    - yaml

possible test cases


tests/test_basicio.py
1. 
    - ref: native 
    - target:  native
    - assignment: str
`admixfrog --ref data/ref_A1240k.csv.xz --infile data/oase_chr9.in.xz --states AFR NEA`
2. 
    - ref: native 
    - target:  bam
    - assignment: formatstr
`admixfrog --ref data/ref_A1240k.csv.xz --infile data/oase_chr9.bam --states AFK=AFR ANC=NEA+DEN --cont AFK`
3. 
    - ref: native 
    - target:  bam
    - assignment: yaml
`admixfrog --ref data/ref_A1240k.csv.xz --infile data/oase_chr9.in.xz --states AFK NEA --pop-file data/pops2.yaml --cont AFK`
=====
4. 
    - ref: geno
    - target: geno
    - assignment: formatstr
`admixfrog --geno data/oase --target Oase1_d --states AFR=Yoruba_d VIN=Vindija_DG`
5. 
    - ref: geno
    - target: native
    - assignment: str
`admixfrog --geno data/oase --infile data/oase_chr9.in.gz --states Yoruba_DG VIN=Vindija_DG`
6. 
    - ref: geno
    - target: bam
    - assignment: yaml
`admixfrog --geno data/oase --bam data/oase_chr9.bam --pop-file data/pops.yaml --states Yoruba_DG VIN=Vindija_DG`
7. 
    - ref: vcf
    - target:  bam
    - assignment: yaml
`admixfrog --vcf-ref data/oase.vcf.gz --bam data/oase_chr9.bam --pop-file data/pops2.yaml --states AFR NEA`
8. 
    - ref: vcf
    - target: vcf
    - assignment: yaml
`admixfrog --vcf-ref data/oase.vcf.gz --vcf-gt data/oase.vcf.gz --pop-file data/pops.yaml --states AFR NEA`
