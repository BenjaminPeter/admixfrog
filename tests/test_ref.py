"""Test running using genotype from vcf
1. create ref file from vcf
2. run program with bam file and created ref
"""


class TestRefVCF(object):
    data='data/oase_chr9.bam'
    ref='res/ref.xz'
    ref_raw='data/oase.vcf.gz'
    final='res/test_refvcf'
    popfile='data/pops.yaml'

    def test_reffile_from_vcf(self, script_runner):

        cmd = f'admixfrog-ref --vcf-ref {self.ref_raw} --rec-rate 1e-8 '
        cmd += f'--chroms 9 --pop-file {self.popfile} --out {self.ref} --states AFR NEA DEN '
        print(cmd)
        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
    
    def test_runwithref(self, script_runner):
        cmd = f'admixfrog --bam {self.data} --states AFR NEA --force-infile '
        cmd += f'-P --out {self.final} -b 100000 --ref {self.ref}'

        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success

class TestRefVCF2(TestRefVCF):
    ref='res/ref2.xz'
    rec='data/maps_chr.9'
    ref_raw='data/oase_small.vcf.gz'
    final='res/test2_refvcf'
    def test_reffile_from_vcf(self, script_runner):

        cmd = f'admixfrog-ref --out {self.ref} --vcf-ref {self.ref_raw} '
        cmd += f'--state-file {self.popfile} '
        cmd += f'--rec-file {self.rec} '
        cmd += f'--state-file {self.popfile} '
        cmd += f'--states AFR NEA=Altai_snpAD.DG ALT '
        cmd += f'--map-id AA_Map deCODE COMBINED_LD '
        cmd += f'--chroms 9 '
        print(cmd)
        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success

## RANDOM READ SAMPLE TEST
class TestRefVCFPseudoHAP(TestRefVCF):
    ref='res/ref2ph.xz'
    rec='data/maps_chr.9'
    ref_raw='data/oase_small.vcf.gz'
    final='res/test2_refvcf_pseudohap'
    def test_reffile_from_vcf(self, script_runner):

        cmd = f'admixfrog-ref --out {self.ref} --vcf-ref {self.ref_raw} '
        cmd += f'--state-file {self.popfile} '
        cmd += f'--rec-file {self.rec} '
        cmd += f'--state-file {self.popfile} '
        cmd += f'--states AFR NEA=Altai_snpAD.DG ALT '
        cmd += f'--map-id AA_Map deCODE COMBINED_LD '
        cmd += f'--chroms 9 '
        cmd += f'--pseudo-haploid Altai_snpAD.DG'
        print(cmd)
        args = cmd.split()
        ret = script_runner.run(args, cwd='tests/')
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
