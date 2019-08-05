#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_pce_py PCE with Combigrids (Python)
## This simple example shows how to create a Polynomial Chaos Expansion from an
## adaptively refined combigrid.

import numpy as np
import pysgpp
import time

## First we have to define a model to approximate.
def expModel(x):
    return np.exp(-x[0] * x[1])

## Then we can create a refined combigrid 
def ct_to_pce():
    start_time = time.time()
    # initialize model function
    func = pysgpp.multiFunc(expModel)
    numDims = 2
    # regular sparse grid level q
    q = 6 
    # create polynomial basis
    config = pysgpp.OrthogonalPolynomialBasis1DConfiguration()
    config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_LEGENDRE
    basisFunction = pysgpp.OrthogonalPolynomialBasis1D(config) 
    # create sprarse grid interpolation operation
    op = pysgpp.CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(numDims, func)
    # start with regular level q and add some levels adaptively
    op.getLevelManager().addRegularLevels(q)
    op.getLevelManager().addLevelsAdaptiveByNumLevels(5)
    
##  and construct a PCE representation to easily calculate statistical features of our model.
    # create polynomial chaos surrogate from sparse grid
    surrogateConfig = pysgpp.CombigridSurrogateModelConfiguration()
    surrogateConfig.type = pysgpp.CombigridSurrogateModelsType_POLYNOMIAL_CHAOS_EXPANSION
    surrogateConfig.loadFromCombigridOperation(op)
    surrogateConfig.basisFunction = basisFunction
    pce = pysgpp.createCombigridSurrogateModel(surrogateConfig)
    # compute sobol indices
    sobol_indices = pysgpp.DataVector(1)
    total_indices = pysgpp.DataVector(1)
    pce.getComponentSobolIndices(sobol_indices)
    pce.getTotalSobolIndices(total_indices)
    # print results
    print("Mean: {} Variance: {}".format(pce.mean(), pce.variance()))
    print("Sobol indices {}".format(sobol_indices.toString()))
    print("Total Sobol indices {}".format(total_indices.toString()))
    print("Sum {}\n".format(sobol_indices.sum()))
    
    print("Elapsed time: {} s".format(time.time() - start_time))

## Output: 
## @verbatim
## Mean: 0.796599599298 Variance: 0.0250607565267
## Sobol indices [4.46123454530588048339e-01, 4.46123454530588825495e-01, 1.07753090938823264944e-01]
## Total Sobol indices [5.53876545469411341038e-01, 5.53876545469412118194e-01]
## Sum 1.0
## @endverbatim

 
if __name__ == "__main__":
    try:
        ct_to_pce()
    except RuntimeError as e:
        if "Eigen" in e.args[0]:
            print("SGpp was built without Eigen support.\nSkipping example...")
            exit(0)
        elif "DAKOTA" in e.args[0]:
            print("SGpp was built without DAKOTA support.\nSkipping example...")
            exit(0)
        else:
            raise(e)
