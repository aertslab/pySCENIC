from . import _version

__version__ = _version.get_versions()['version']

import logging

from pyscenic.log import create_logging_handler

LOGGER = logging.getLogger(__name__)
# Set logging level.
logging_debug_opt = False
LOGGER.addHandler(create_logging_handler(logging_debug_opt))
LOGGER.setLevel(logging.DEBUG)
