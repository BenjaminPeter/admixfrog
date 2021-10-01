from .interface_slug import run_sfs
from .interface_frog import run_frog, do_rle
from .interface_io import bam, bam2, do_ref
from pkg_resources import get_distribution, DistributionNotFound
from .interface_frog import profile as profile_frog
from .interface_slug import profile as profile_slug

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "v0"
