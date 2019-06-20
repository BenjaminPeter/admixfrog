def test_version(script_runner):
    ret = script_runner.run('admixfrog', '--version')
    assert ret.success

def test_run_from_bam_ref(script_runner):
    cmd = 'admixfrog --bam tests/data/oase_chr9.bam --ref tests/data/ref_A1240k.csv.xz '
    cmd += ' --out tests/res/test_bam2 --force-infile --states AFR NEA -b 100000 -P'
    args = cmd.split()
    print(args)
    ret = script_runner.run(*args)
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success



