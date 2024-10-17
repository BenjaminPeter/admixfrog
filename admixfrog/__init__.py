from .interface_slug import run_sfs
from .interface_frog import run_frog, do_rle
from .interface_io import bam, bam2, do_ref
from .interface_frog import profile as profile_frog
from .interface_slug import profile as profile_slug
import importlib


__version__ = importlib.metadata.version(__name__)


