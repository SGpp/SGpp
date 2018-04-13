"""
Quadrature package
==========================================

"""

__version__ = "0.1"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from sparse_grid import doQuadrature, getIntegral, getIntegralOfBasisFunction
from HashQuadrature import HashQuadrature

import strategies
import linearform
import bilinearform
import trilinearform
import marginalization
