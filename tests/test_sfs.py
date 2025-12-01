def test_basic(script_runner):
    """test case 1"""
    cmd = "admixslug --infile data/oase_chr9_sfs.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_sfs --seed 13 --force-infile --states AFR NEA "
    cmd += " --output-vcf --autosomes"
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_jk(script_runner):
    """test case 1"""
    cmd = "admixslug --infile data/oase_chr9_sfs.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_sfs --seed 13 --force-infile --jk-resamples 3 --states AFR NEA  "
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_slug_input(script_runner):
    """case where we generate input file"""
    cmd = "admixfrog-bam2 --bam data/oase_chr{CHROM}.bam "
    cmd += " --ref data/ref_A1240k.csv.xz --out data/Oase.sfs.in.xz "
    cmd += " --chroms 9,X --force-bam"
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_slug_input_plus(script_runner):
    """case where we run admixslug from bam"""
    cmd = "admixslug --bam data/oase_chr{CHROM}.bam "
    cmd += " --ref data/ref_A1240k.csv.xz --out res/test_sfs2 "
    cmd += " --chroms 9,X --force-bam --states AFR NEA"
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success


def test_binning(script_runner):
    """test case 2 with binning"""
    cmd = "admixslug --infile data/oase_chr9_sfs.in.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_binning --seed 13 --force-infile --states AFR NEA  "
    cmd += " --deam-bin-size 100 --len-bin-size 20"
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success

def test_gtmode(script_runner):
    cmd = "admixslug --infile data/chag_gt.csv.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_sfsgt --seed 13 --force-infile --states CHA  --gt-mode "
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)

def test_gtmode_jk(script_runner):
    cmd = "admixslug --infile data/chag_gt.csv.xz --ref data/ref_A1240k.csv.xz "
    cmd += " --out res/test_sfsgtjk --seed 13 --force-infile --states CHA --gt-mode --jk 10 --output-f "
    args = cmd.split()
    print(args)
    ret = script_runner.run(args, cwd="tests")
    print(ret.stdout)
