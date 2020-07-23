from .admixfrog import run_admixfrog
from .interface import run, bam, do_rle, profile, do_ref, run_sfs
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "v0"
