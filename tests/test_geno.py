"""Test running using genotype from vcf
1. create ref file from vcf
2. run program with bam file and created ref
"""


class TestGeno(object):
    data = "data/oase"
    final = "res/test_geno"
    target = "Oase1_d"
    states = "AFR=Yoruba.DG NEA=Vindija.DG+Altai.DG"

    def test_geno(self, script_runner):
        sid = self.target

        cmd = f"admixfrog --gfile {self.data} --guess-ploidy -P --bin-size 100000"
        cmd += f" --states {self.states} --out {self.final} --target {self.target}"
        cmd += f" --dont-est-contamination  --ancestral Denisova.DG --c0 0 "
        args = cmd.split()
        ret = script_runner.run(*args, cwd="tests/")
        print(ret.stdout)
        print(ret.stderr)
        assert ret.success
