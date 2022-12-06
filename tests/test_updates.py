
def test_update_F(script_runner):
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_filter --seed 13 --force-infile '
    cmd += ' --states AFR NEA -b100000 -P'
    cmd += ' --est-F --est-error '
    args = cmd.split()
    print(args)
    ret = script_runner.run(*args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_update_Ftau(script_runner):
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_filter --seed 13 --force-infile '
    cmd += ' --states AFR NEA -b100000 -P'
    cmd += ' --est-F --est-tau --F0 0.3 0.4 --tau0 .5 .1'
    args = cmd.split()
    print(args)
    ret = script_runner.run(*args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_update_all(script_runner):
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_filter --seed 13 --force-infile '
    cmd += ' --states AFR NEA -b100000 -P'
    cmd += ' --est-F --est-tau --est-error'
    args = cmd.split()
    print(args)
    ret = script_runner.run(*args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
