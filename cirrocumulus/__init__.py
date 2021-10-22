import logging
import os

from cirrocumulus.envir import CIRRO_LOG_LEVEL

logger = logging.getLogger('cirro')
log_level = os.environ.get(CIRRO_LOG_LEVEL, 'ERROR')
logger.propagate = False
logger.setLevel(log_level)
logger.addHandler(logging.StreamHandler())  # Logs go to stderr
logger.handlers[-1].setFormatter(logging.Formatter("%(message)s"))
logger.handlers[-1].setLevel(log_level)
