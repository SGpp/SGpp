import os
import sys

# add current directory to PYTHONPATH such that pysgpp_swig can be imported
sys.path.append(os.path.dirname(__file__))


is_windows = sys.platform.startswith('win')
if is_windows and sys.version_info >= (3, 8):
    # add current directory and such that pysgpp_swig can be imported
    os.add_dll_directory(os.path.dirname(__file__))
    # add sgpp lib path such that other sgpp dlls can be imported
    os.add_dll_directory(os.path.join(os.environ['WINSGPP_PATH'], "lib", "Release"))

# import pysgpp_swig and extensions
from .pysgpp_swig import *
from . import extensions
