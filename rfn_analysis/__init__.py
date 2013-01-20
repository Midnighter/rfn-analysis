# -*- coding: UTF-8 -*-


import subprocess

def get_version():
    version = "0.1.dev"
    version += subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).strip()
    return version

__version__ = get_version()

from .classes import *
from .extraction import *

