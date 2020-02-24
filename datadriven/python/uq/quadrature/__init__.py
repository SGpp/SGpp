# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
Quadrature package
==========================================

"""

__version__ = "0.1"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from pysgpp.extensions.datadriven.uq.quadrature.sparse_grid import doQuadrature, getIntegral, getIntegralOfBasisFunction
from pysgpp.extensions.datadriven.uq.quadrature.HashQuadrature import HashQuadrature

from pysgpp.extensions.datadriven.uq.quadrature import strategies
from pysgpp.extensions.datadriven.uq.quadrature import linearform
from pysgpp.extensions.datadriven.uq.quadrature import bilinearform
from pysgpp.extensions.datadriven.uq.quadrature import trilinearform
from pysgpp.extensions.datadriven.uq.quadrature import marginalization
