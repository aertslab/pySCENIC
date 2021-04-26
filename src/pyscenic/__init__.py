from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

import logging
from pyscenic.log import create_logging_handler

LOGGER = logging.getLogger(__name__)
# Set logging level.
logging_debug_opt = False
LOGGER.addHandler(create_logging_handler(logging_debug_opt))
LOGGER.setLevel(logging.DEBUG)
