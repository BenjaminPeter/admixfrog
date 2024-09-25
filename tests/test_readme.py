def test_version(script_runner):
    ret = script_runner.run(['admixfrog', '--version'])
    assert ret.success

def test_quickstart(script_runner):
    cmd = 'admixfrog --gfile data/oase --target Oase1_d --states NEA=Vindija.DG+Altai.DG YRI=Yoruba.DG  Denisova.DG  --cont YRI'
    args = cmd.split()
    ret = script_runner.run(args, cwd='tests')
    assert ret.success


