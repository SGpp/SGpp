"""
Quadrature package
==========================================

"""
from __future__ import absolute_import

__version__ = "0.1"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from .sparse_grid import doQuadrature, getIntegral, getIntegralOfBasisFunction
from .HashQuadrature import HashQuadrature

from . import strategies
from . import linearform
from . import bilinearform
from . import trilinearform
from . import marginalization
