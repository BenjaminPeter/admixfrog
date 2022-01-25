import sys
import logging

log_ = logging


def setup_log(filename=None):
    log_.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s]:  %(message)s",
        filename=filename,
        filemode="w",
        datefmt="%H:%M:%S",
    )
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
