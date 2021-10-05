import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%m-%d %H:%M:%S',
    stream=sys.stderr)

# Load all modules
from .userConfig import *
from .preparation import *
from .data import *
from .features import *
from .prediction import *
from .main import *
