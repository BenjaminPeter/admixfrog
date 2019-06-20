"""Test running using genotype from vcf
1. create input file from vcf and a reference
2. run program
"""


class TestVCFGT(object):
    out='res/oase_from_vcf.in.xz'
    ref='data/ref_A1240k.csv.xz'
    final='res/test_vcf'

    def test_infile_from_vcfgt(self, script_runner):
        vcfgt='data/oase.vcf.gz'
        sid='Oase1_d'

        cmd = f'admixfrog-bam --vcf-gt {vcfgt} --sample-id {sid} --ref {self.ref} '
        cmd += f'--chroms 9 --random-read-sample --force-infile --out {self.out}'
        args = cmd.split()
        ret = script_runner.run(*args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
    
    def test_created_infile(self, script_runner):
        cmd = f'admixfrog --infile {self.out} --states AFR NEA '
        cmd += f'-P --out {self.final} -b 10000 --ref {self.ref}'

        args = cmd.split()
        ret = script_runner.run(*args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success

