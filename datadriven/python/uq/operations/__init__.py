# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
Sparse Grid operations package
==========================================

"""

__version__ = "0.1"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from pysgpp.extensions.datadriven.uq.operations.discretization import discretize, discretizeFunction
from pysgpp.extensions.datadriven.uq.operations.discretizeProduct import discretizeProduct
from pysgpp.extensions.datadriven.uq.operations.epsilonComplexity import getL2EpsilonComplexity
from pysgpp.extensions.datadriven.uq.operations.general import (isNumerical,
                     isList,
                     isMatrix,
                     extend_grid_1d,
                     extend_grid,
                     join,
                     project)
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import (bsplineBoundaryGridTypes,
                         bsplineNoBoundaryGridTypes,
                         bsplineGridTypes,
                         polyBoundaryGridTypes,
                         polyNoBoundaryGridTypes,
                         polyGridTypes,
                         linearBoundaryGridTypes,
                         linearNoBoundaryGridTypes,
                         linearGridTypes,
                         getBoundsOfSupport,
                         sub,
                         add,
                         addConst,
                         getHierarchicalAncestors,
                         insertPoint,
                         insertHierarchicalAncestors,
                         insertTruncatedBorder,
                         hasBorder,
                         parent, parents,
                         getBasis, getDegree,
                         loadOperationMultiEval,
                         evalSGFunction,
                         evalSGFunctionMulti,
#                          evalSGFunctionMultiVectorized,
                         isValid,
                         isRefineable,
                         copyGrid, createGrid,
                         hierarchize, dehierarchize, dehierarchizeList,
                         hierarchizeBruteForce,
                         hierarchizeEvalHierToTop,
                         balance,
                         estimateSurplus, estimateConvergence,
                         checkInterpolation,
                         checkPositivity,
                         hasChildren,
                         hasAllChildren)

# from natafTransformation import NatafTransformation
