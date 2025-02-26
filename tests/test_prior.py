def test_prior_fix(script_runner):
    """test case 1"""
    cmd = "admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_basic --seed 13 --force-infile --states AFR NEA -b100000 -P"
    cmd += " --prior 1"
    args = cmd.split()
    print(" ".join(args))
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_prior_fix_anc(script_runner):
    """test case 1"""
    cmd = "admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_basic --seed 13 --force-infile --states AFR NEA -b100000 -P"
    cmd += " --prior 1 --ancestral-prior 1. --ancestral PAN"
    args = cmd.split()
    print(" ".join(args))
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_prior_eb(script_runner):
    """test case 1"""
    cmd = "admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_basic --seed 13 --force-infile --states AFR NEA -b100000 -P"
    cmd += " --ancestral PAN "
    args = cmd.split()
    print(" ".join(args))
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_prior_eb_anc(script_runner):
    """test case 1"""
    cmd = "admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_basic --seed 13 --force-infile --states AFR NEA -b100000 -P"
    cmd += " --ancestral PAN --ancestral-prior 1."
    args = cmd.split()
    print(" ".join(args))
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
