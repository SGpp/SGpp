#!/usr/bin/env python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_diagonalRegularizationExample_py Gaussian Weight Priors 
##
## This example compares two different Gaussian priors for sparse grid regression.
## The first one is the standard prior that assumes a constant, identical variance
## for each weight.
## The second, improved prior uses the additional assumption that the true distribution
## function is sufficiently smooth and then imposes a multivariate Gaussian with the
## covariance matrix
## \f$\mathbf{\Gamma}_{i,i} = 0.25^{\vert \mathbf{l} \vert_1 - d}\f$, where
## \f$ d \f$ corresponds to the dimension of the grid and
## \f$ \mathbf{\vert \mathbf{l} \vert_1} \f$ to the level sum of the ith grid point.

from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
import requests as r
import numpy as np
import pandas as pd
import sklearn.preprocessing as pre
from sklearn.cross_validation import KFold
from scipy import stats
from zipfile import ZipFile
import io 
import pysgpp as sg; sg.omp_set_num_threads(4)

## This function scales all predictors so that they are suitable for sparse grids.
def scale(df, scaler=None):
    Y = df.ix[:,-1] # save Y (don't need to transform it/useless for cat. data!)
    X = df.values
    if scaler:
        X = scaler.transform(X)
    else:
        scaler = pre.MinMaxScaler()
        X = scaler.fit_transform(X)
    index = df.index
    columns = df.columns
    df = pd.DataFrame(data=X, index=index, columns=columns)
    df.ix[:,-1] = Y
    return scaler, df

## This function performs a Box-Cox transformation, which (hopefully) results in data
## that is better distributed.
def transform_cox(df, lambdas):
    scaler = pre.MinMaxScaler()
    for variable in lambdas:
        lamb = lambdas[variable]
        if lamb == 1:
            continue # identity transform
        data_trans = stats.boxcox(df[variable] + 10e-1)
        df[variable] = scaler.fit_transform(np.array(data_trans[0]).reshape(-1, 1))
    return df

## This function downloads the PowerPlant dataset and performs the necessary preprocessing steps.
def get_dataset():
    ## Download and unzip the dataset.
    data_url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00294/CCPP.zip"
    print("Loading power plant dataset from the UCI repository.")
    resp = r.get(data_url, stream=True)
    data = ZipFile(io.BytesIO(resp.content))
    with data.open('CCPP/Folds5x2_pp.xlsx') as xls:
        df =  pd.read_excel(xls)
    print("Preprocessing dataset.")
    ## Then scale and transform it.
    _, df = scale(df)
    lambdas = {'AP': 0,
            'AT': 1.3353837296219406,
            u'PE': 1,
            'RH': 2.4604945158589104,
            'V': -0.43989911518156471}
    df = transform_cox(df, lambdas)   
    return df

## This function returns a sparse grid regressing learner.
def make_estimator(lambda_reg, prior):
    grid = sg.RegularGridConfiguration()
    grid.dim_ = 4
    grid.level_ = 5
    grid.type_ = sg.GridType_ModLinear

    adapt = sg.AdpativityConfiguration()
    adapt.numRefinements_ = 5
    adapt.noPoints_ = 3

    solv = sg.SLESolverConfiguration()
    solv.maxIterations_ = 50
    solv.eps_ = 10e-6
    solv.threshold_ = 10e-6
    solv.type_ = sg.SLESolverType_CG

    final_solv = solv
    final_solv.maxIterations = 200

    regular = sg.RegularizationConfiguration()
    regular.type_ = sg.RegularizationType_Diagonal
    regular.exponentBase_ = prior
    regular.lambda_ = lambda_reg    

    estimator = sg.RegressionLearner(grid, adapt, solv, final_solv,regular)
    return estimator

## This function returns the testing error for one estimator.
def evaluate_one(estimator, X, y, train, test):
    train_X = sg.DataMatrix(X[train])
    train_y = sg.DataVector(y[train])
    test_X = sg.DataMatrix(X[test])
    test_y = sg.DataVector(y[test])
    estimator.train(train_X,train_y)
    error = estimator.getMSE(test_X, test_y)
    return error

## This function returns the cv-error for one estimator.
def evaluate(estimator, cv, X, y):
    errors = []
    for (train, test) in cv:
        error = evaluate_one(estimator, X, y, train, test)
        errors.append(error)
    return np.mean(errors)

## This function performs a grid search for the best regularization parameters.
def grid_search(X, y, prior):
    cv = KFold(X.shape[0], 10)
    lambda_grid = np.logspace(-9, -4, num=4)
    best_result = np.finfo(np.double).max
    for l in lambda_grid:
        estimator = make_estimator(l, 1.0)
        mse = evaluate(estimator, cv, X, y)
        print("Prior={} with lambda={:2.4e} achieved a RMSE of {:2.4e}.".format(prior, l, np.sqrt(mse)))
        best_result = min(best_result, mse)
    return np.sqrt(best_result)
    
def main():
    df = get_dataset()
    X = np.array(df.ix[:,0:-1])
    y = (df.ix[:,-1]).values
    ## A parameter of 1 corresponds to the standard ridge prior, one of 0.25 to the improved prior.
    best_identity = grid_search(X, y, 1.0)
    best_improved = grid_search(X, y, 0.25)
    print("\nThe identity matrix achieved an RMSE of {:3.5f}, the diagonal of {:3.5f}.".format(best_identity, best_improved))
    
if __name__ == '__main__':
    main()
