def test_inbreeding(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_basic --seed 13 --force-infile --states AFR NEA -b 100000 -P'
    cmd += ' --est-inbreeding ' 
    args = cmd.split()
    print(args)
    ret = script_runner.run(*args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
