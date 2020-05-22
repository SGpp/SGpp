#!/usr/bin/env python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_generalisedGridsest_py Generalised Sparse Grids
##
## This example tests generalised sparse grids.       
## It generates a Friedman1 dataset and then compares 
## the performance of estimators with various grid    
## granularities.                                     
## The grid granularities are controlled by the parameter \f$T\f$; the number of grid
## points is then given by \f$G_n^T\f$ and  the approximation space is given by \f$V_n^T\f$.
## \f{align}
## G_n^T &= \bigcup_{\substack{\vert {\mathbf{l}} \vert_1 - T \vert \mathbf{i} \vert_\infty \\ \leq n + d - 1 - T n}} G_{\mathbf{l}},\\
## V_n^T &= \bigoplus_{\substack{\vert {\mathbf{l}} \vert_1 - T \vert \mathbf{i} \vert_\infty \\ \leq n + d - 1 - T n}} W_{\mathbf{l}}\nonumber \f}

try:
    ## We first import all pysgpp and other utility libraries.
    import numpy as np
    import pysgpp as sg; sg.omp_set_num_threads(4)
    #import pandas as pd
    import sklearn.datasets as data

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)

## This function generates the Friedman1 dataset on demand.
def generate_friedman1(seed):
    (X,y) = data.make_friedman1(n_samples=10000, random_state=seed, noise=1.0)
    y = sg.DataVector(y)   
    X = sg.DataMatrix(X)
    return X, y

## This function evaluates the performance of a learner with standard settings
## and different values of T.
def evaluate(X_tr, y_tr, X_te, y_te, T):
    gridConfig = sg.RegularGridConfiguration()
    gridConfig.dim_ = 10
    gridConfig.level_ = 4
    gridConfig.t_ = T
    gridConfig.type_ = sg.GridType_ModLinear

    adaptivityConfig = sg.AdaptivityConfiguration()
    adaptivityConfig.numRefinements_ = 5
    adaptivityConfig.numRefinementPoints_ = 3
    adaptivityConfig.numCoarseningPoints_ = 3

    solv = sg.SLESolverConfiguration()
    solv.maxIterations_ = 50
    solv.eps_ = 1e-5
    solv.threshold_ = 1e-5
    solv.type_ = sg.SLESolverType_CG

    final_solv = solv
    final_solv.maxIterations = 200

    regularizationConfig = sg.RegularizationConfiguration()
    regularizationConfig.type_ = sg.RegularizationType_Identity
    regularizationConfig.exponentBase_ = 1.0
    regularizationConfig.lambda_ = 1e-3

    ## Create the estimator, train it with the training data and then return the error
    ## for the testing set.
    estimator = sg.RegressionLearner(gridConfig, adaptivityConfig, solv, final_solv,regularizationConfig)
    estimator.train(X_tr,y_tr)
    print(estimator.getGridSize())
    return estimator.getMSE(X_te,y_te)

def main():
    ## First generate the training and test data.
    X_tr, y_tr = generate_friedman1(123456)    
    X_te, y_te = generate_friedman1(345678)

    ## Then we evaluate the testing error for \f$ T \in \{-0.5, 0, 0.5, 1.0\} \f$.
    Ts = [-0.5, 0, 0.5, 1.0]
    for T in Ts:
        mse = evaluate(X_tr, y_tr, X_te, y_te, T)
        print("The sparse grid with T={:2.1f} achieved a testing RMSE of {:2.4f}.".format(T, np.sqrt(mse)))

if __name__ == '__main__':
    main()
