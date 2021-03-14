import os
import sys

# add current directory to PYTHONPATH such that pysgpp_swig can be imported
sys.path.append(os.path.dirname(__file__))

# add current directory and sgpp lib path such that pysgpp_swig and other sgpp dlls can be imported
os.add_dll_directory(os.path.dirname(__file__))
os.add_dll_directory(os.path.join(os.environ['WINSGPP_PATH'], "lib", "Debug"))

# import pysgpp_swig and extensions
from .pysgpp_swig import *
from . import extensions
