from .introgression_bb import run_hmm_bb
from .interface import run, bam, do_rle, profile
from pkg_resources import get_distribution, DistributionNotFound

try:
        __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "v0"

