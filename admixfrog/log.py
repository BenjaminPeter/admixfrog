import logging
log_ = logging.getLogger(__name__)

def setup_log():
    logger = logging.getLogger('admixfrog')
    logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('[%(asctime)s]: %(message)s')
    ch.setFormatter(formatter)

    logger.addHandler(ch)
    return logger
