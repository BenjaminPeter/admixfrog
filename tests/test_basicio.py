def test_version(script_runner):
    ret = script_runner.run('admixfrog', '--version')
    assert ret.success

def test_basic(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k_9X.csv.xz '
    cmd += ' --out res/test_basic --seed 13 --force-infile  '
    cmd += '--states AFR NEA -b 100000 -P'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_autosomes(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k_9X.csv.xz '
    cmd += ' --out res/test_basic --seed 13 --force-infile  '
    cmd += '--states AFR NEA -b 100000 -P --autosomes'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_twostates(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k_9X.csv.xz '
    cmd += ' --out res/test_basic --seed 13 --force-infile  '
    cmd += '--states AFR NEA -b 100000 -P --homo-states AFR --het-states AFRNEA'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_snpmode(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k_9X.csv.xz '
    cmd += ' --out res/test_basic --seed 13 '
    cmd += '--states AFR NEA -b 100000 --snp-mode'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_run_from_bam_rename(script_runner):
    cmd = 'admixfrog --bam tests/data/oase_chr9.bam --ref tests/data/ref_A1240k.csv.xz '
    cmd += ' --out tests/res/test_bam2 --force-infile -b 100000 -P'
    cmd += ' --states AFK ARC=NEA+DEN --cont-id AFK'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args)
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_run_from_bam2_rename(script_runner):
    cmd = 'admixfrog --bam tests/data/oase_chr9.bam --ref tests/data/ref_A1240k.csv.xz '
    cmd += ' --out tests/res/test_bam22 --force-infile -b 100000 -P'
    cmd += ' --states AFK ARC=NEA+DEN --cont-id AFK --bin-reads --len-bin-size 100'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args)
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_run_from_bam_yaml(script_runner):
    cmd = 'admixfrog --bam data/oase_chr9.bam --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_bam3 --force-infile -b 100000 -P'
    cmd += ' --states AFK=AFR ARC=NEA+DEN --cont-id AFK'
    cmd += ' --pop-file data/pops2.yaml'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_run_from_bam_tmat(script_runner):
    cmd = 'admixfrog --bam data/oase_chr9.bam --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_bam3 --force-infile -b 100000 '
    cmd += ' --states AFK=AFR ARC=NEA+DEN --cont-id AFK'
    cmd += ' --pop-file data/pops2.yaml'
    cmd += ' --tmat data/tmat.txt'
    cmd += ' --dont-est-trans --dont-est-cont'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_fake_cont(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_basic --seed 13 --force-infile --states AFR NEA -b100000 -P'
    cmd += ' --fake-contamination 0.4 --downsample 0.1'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_filters(script_runner):
    """test case 1"""
    cmd = 'admixfrog --infile data/oase_chr9.in.xz --ref data/ref_A1240k.csv.xz '
    cmd += ' --out res/test_filter --seed 13 --force-infile '
    cmd += ' --states AFR NEA -b100000 -P'
    cmd += ' --filter-map 0 --filter-pos 50 --filter-delta .1'
    cmd += ' --filter-ancestral --ancestral PAN'
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd='tests')
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

