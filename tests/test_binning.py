def test_frog_binning(script_runner):
    """test case with binned reads"""
    cmd = "admixfrog --infile data/oase_chr9_sfs.in.xz --ref data/ref_A1240k_9X.csv.xz "
    cmd += " --out res/test_binning_frog --seed 13   "
    cmd += "--states AFR NEA -b 100000 -P --bin-reads"
    cmd += " --deam-bin-size 100 --len-bin-size 150"
    args = cmd.split()
    print(args)
    ret = script_runner.run(*args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
