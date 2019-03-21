from logging import Formatter
import logging
log_ = logging.getLogger(__name__)

def setup_log():
    #logger = logging.getLogger('admixfrog')
    log_.setLevel(logging.DEBUG)

    formatter = logging.Formatter('[%(asctime)s]: %(message)s')

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)

    log_.addHandler(ch)
    return log_
