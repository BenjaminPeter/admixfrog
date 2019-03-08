import logging
log_ = logging.getLogger(__name__)

def setup_log():
    logger = logging.getLogger('admixfrog')
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    formatter = logging.Formatter('[%(asctime)s]: %(message)s')
    ch.setFormatter(formatter)

    logger.addHandler(ch)
    return logger
