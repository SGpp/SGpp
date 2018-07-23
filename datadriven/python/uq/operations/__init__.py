"""
Sparse Grid operations package
==========================================

"""

__version__ = "0.1"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from discretization import discretize, discretizeFunction
from discretizeProduct import discretizeProduct
from epsilonComplexity import getL2EpsilonComplexity
from general import (isNumerical,
                     isList,
                     isMatrix,
                     extend_grid_1d,
                     extend_grid,
                     join,
                     project)
from sparse_grid import (bsplineBoundaryGridTypes,
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
