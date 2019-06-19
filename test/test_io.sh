#test bam-stuff
samtools index data/oase_chr9.bam
admixfrog --bam data/oase_chr9.bam \
    --ref data/ref_A1240k.csv.xz \
    --out res/test_bam \
    --force-infile \
    -b 100000
Rscript scripts/test_bam_plot.R res/test_bam 0.1


#create reference
admixfrog-ref --vcf-file data/oase.vcf.gz \
    --rec-rate 1e-8                       \
    --out res/ref.xz                      \
    --chroms 9                            \
    --pop-file data/pops.yaml

admixfrog --bam data/oase_chr9.bam \
    --ref res/ref.xz               \
    --out res/test_bam2            \
    --force-infile                 \
    --states AFR NEA               \
    -b 100000 -P

Rscript scripts/test_bam_plot.R res/test_bam2 100000

#vcf-gt
bcftools index data/oase.vcf.gz
admixfrog-bam --vcf-gt data/oase.vcf.gz \
    --force-infile                      \
    --chroms 9                          \
    --random-read-sample                \
    --sample-id Oase1_d                 \
    --out res/oase_from_vcf.in.xz       \
    --ref data/ref_A1240k.csv.xz

admixfrog --infile res/oase_from_vcf.in.xz \
    -b 10000                               \
    --out res/test_vcf                     \
    -P                                     \
    --states AFR NEA                       \
    --ref res/ref.xz

Rscript scripts/test_bam_plot.R res/test_vcf 10000
