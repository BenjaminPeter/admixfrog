"""Test running using genotype from vcf
1. create input file from vcf and a reference
2. run program
"""


class TestVCFGT(object):
    out='res/oase_from_vcf.in.xz'
    out_yri='res/yri_from_vcf.in.xz'
    out_yri2='res/yri2_from_vcf.in.xz'
    ref='data/ref_A1240k.csv.xz'
    final='res/test_vcf'
    rle = 'res/test_vcf_rle.xz'

    def test_infile_from_vcfgt(self, script_runner):
        vcfgt='data/oase.vcf.gz'
        sid='Oase1_d'

        cmd = f'admixfrog-bam --vcf-gt {vcfgt} --sample-id {sid} --ref {self.ref} '
        cmd += f'--chroms 9 --random-read-sample --force-infile --out {self.out}'
        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success

    def test_infile_from_vcfgt_allele2(self, script_runner):
        vcfgt='data/oase.vcf.gz'
        sid='B_Yoruba-3.DG'

        cmd = f'admixfrog-bam --vcf-gt {vcfgt} --sample-id {sid} --ref {self.ref} '
        cmd += f'--chroms 9 --random-read-sample --force-infile --out {self.out_yri}'
        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
    def test_infile_from_vcfgt_allele3(self, script_runner):
        vcfgt='data/oase.vcf.gz'
        sid='B_Yoruba-3.DG'

        cmd = f'admixfrog-bam --vcf-gt {vcfgt} --sample-id {sid} --ref {self.ref} '
        cmd += f'--chroms 9 --force-infile --out {self.out_yri2}'
        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
    
    def test_created_infile(self, script_runner):
        cmd = f'admixfrog --infile {self.out} --states AFR NEA --gt-mode '
        cmd += f'-P --out {self.final} -b 100000 --ref {self.ref}'

        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success

    def test_rle(self, script_runner):
        cmd = f'admixfrog-rle --in {self.final}.bin.xz --run-penalty 0.23 '
        cmd += f'--out {self.rle} '

        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
