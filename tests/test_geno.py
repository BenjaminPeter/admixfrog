"""Test running using genotype from vcf
1. create ref file from vcf
2. run program with bam file and created ref
"""


class TestGeno(object):
    data='data/oase'
    final='res/test_geno'

    def test_reffile_from_vcf(self, script_runner):
        vcfgt='data/oase.vcf.gz'
        sid='Oase1_d'

        cmd = f'admixfrog-ref --vcf-file {self.ref_raw} --rec-rate 1e-8 '
        cmd += f'--chroms 9 --pop-file {self.popfile} --out {self.ref} '
        args = cmd.split()
        ret = script_runner.run(*args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
    
    def test_runwithref(self, script_runner):
        cmd = f'admixfrog --bam {self.data} --states AFR NEA --force-infile '
        cmd += f'-P --out {self.final} -b 100000 --ref {self.ref}'

        args = cmd.split()
        ret = script_runner.run(*args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success

