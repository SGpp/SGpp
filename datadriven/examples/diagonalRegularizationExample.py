#!/usr/bin/env python

######################################################################
# This examples tests the performance of the diagonal regularization #
# and the identity matrix for the power plant dataset.               #
# To do this, it performs a grid search for both Tikhonov matrices.  #
######################################################################

import requests as r
import numpy as np
import pandas as pd
import sklearn.preprocessing as pre
from sklearn.cross_validation import KFold
from scipy import stats
from zipfile import ZipFile
from StringIO import StringIO
import pysgpp as sg; sg.omp_set_num_threads(4)

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

def transform_cox(df, lambdas):
    scaler = pre.MinMaxScaler()
    for variable in lambdas:
        lamb = lambdas[variable]
        if lamb == 1:
            continue # identity transform
        data_trans = stats.boxcox(df[variable] + 10e-1)
        df[variable] = scaler.fit_transform(np.array(data_trans[0]).reshape(-1, 1))
    return df

def get_dataset():
    data_url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00294/CCPP.zip"
    print "Loading power plant dataset from the UCI repository."
    resp = r.get(data_url, stream=True)
    data = ZipFile(StringIO(resp.content))
    with data.open('CCPP/Folds5x2_pp.xlsx') as xls:
        df =  pd.read_excel(xls)

    print "Preprocessing dataset."
    _, df = scale(df)
    lambdas = {'AP': 0,
            'AT': 1.3353837296219406,
            u'PE': 1,
            'RH': 2.4604945158589104,
            'V': -0.43989911518156471}
    df = transform_cox(df, lambdas)   
    return df

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
    regular.regType_ = sg.RegularizationType_Diagonal
    regular.exponentBase_ = prior
    regular.lambda_ = lambda_reg    

    estimator = sg.RegressionLearner(grid, adapt, solv, final_solv,regular)
    return estimator

def evaluate_one(estimator, X, y, train, test):
    train_X = sg.DataMatrix(X[train])
    train_y = sg.DataVector(y[train])
    test_X = sg.DataMatrix(X[test])
    test_y = sg.DataVector(y[test])
    estimator.train(train_X,train_y)
    error = estimator.getMSE(test_X, test_y)
    return error

def evaluate(estimator, cv, X, y):
    errors = []
    for (train, test) in cv:
        error = evaluate_one(estimator, X, y, train, test)
        errors.append(error)
    return np.mean(errors)

def grid_search(X, y, prior):
    cv = KFold(X.shape[0], 10)
    lambda_grid = np.logspace(-9, -4, num=4)
    best_result = np.finfo(np.double).max
    for l in lambda_grid:
        estimator = make_estimator(l, 1.0)
        mse = evaluate(estimator, cv, X, y)
        print "Prior={} with lambda={:2.4e} achieved a RMSE of {:2.4e}.".format(prior, l, np.sqrt(mse))
        best_result = min(best_result, mse)
    return np.sqrt(best_result)
    
def main():
    df = get_dataset()
    X = np.array(df.ix[:,0:-1])
    y = (df.ix[:,-1]).values
    best_identity = grid_search(X, y, 1.0)
    best_improved = grid_search(X, y, 0.25)
    print "\nThe identity matrix achieved an RMSE of {:3.5f}, the diagonal of {:3.5f}.".format(best_identity, best_improved)

if __name__ == '__main__':
    main()
