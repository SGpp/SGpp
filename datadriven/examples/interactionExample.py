#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_interactionExample_py Interaction Terms Aware Sparse Grids
##
## This example compares standard sparse grids with sparse grids that only contain
## a subset of all possible interaction terms.
## It uses the optical digits dataset as an example.

import numpy as np
import pysgpp as sg; sg.omp_set_num_threads(4)
import pandas as pd
import sklearn.preprocessing as pre

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

## This function downloads the optical digits dataset and performs the necessary preprocessing steps.
def get_dataset():
    train_url = "https://archive.ics.uci.edu/ml/machine-learning-databases/optdigits/optdigits.tra"
    test_url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/optdigits/optdigits.tes'   
    print("Loading dataset from UCI repository.")
    columns = ["x{}".format(i) for i in range(0, 64)] + ['digit']
    df_train = pd.read_csv(train_url, header=None, index_col=None)
    df_test = pd.read_csv(train_url, header=None, index_col=None)
    df_train.columns=columns
    df_test.columns=columns
    print("Preprocessing dataset.")
    df_complete = df_train.append(df_test, ignore_index=True)
    scaler , _ = scale(df_complete)
    _, df_train = scale(df_train, scaler)
    _, df_test = scale(df_test, scaler)
    return df_train, df_test

## This function evaluates a sparse grid learner for a different set of
## interaction terms.
## To do this, it first trains a classification learner with the training set
## and then evaluates it using the testing part of the dataset.
def evaluate(X_tr, y_tr, X_te, y_te, interactions=None):
    grid = sg.RegularGridConfiguration()
    grid.dim_ = 64
    grid.level_ = 2
    grid.type_ = sg.GridType_ModLinear

    adapt = sg.AdaptivityConfiguration()
    adapt.numRefinements_ = 0
    adapt.noPoints_ = 0

    solv = sg.SLESolverConfiguration()
    solv.maxIterations_ = 50
    solv.eps_ = 10e-6
    solv.threshold_ = 10e-6
    solv.type_ = sg.SLESolverType_CG

    final_solv = solv
    final_solv.maxIterations = 200

    regular = sg.RegularizationConfiguration()
    regular.type_ = sg.RegularizationType_Identity
    regular.exponentBase_ = 1.0
    regular.lambda_ = 0.1    

    X_tr = sg.DataMatrix(X_tr)
    y_tr = sg.DataVector(y_tr)
    X_te = sg.DataMatrix(X_te)
    y_te = sg.DataVector(y_te)
    
    if interactions is None:
        estimator = sg.ClassificationLearner(grid, adapt, solv, final_solv,regular)
    else:
        estimator = sg.ClassificationLearner(grid, adapt, solv, final_solv,regular, interactions)
    estimator.train(X_tr,y_tr)
    return estimator.getAccuracy(X_te,y_te)

def main():
    df_tr, df_te = get_dataset()
    X_tr = np.array(df_tr.ix[:,0:-1])
    y_tr = (df_tr.ix[:,-1]).values
    X_te = np.array(df_te.ix[:,0:-1])
    y_te = (df_te.ix[:,-1]).values

    ## We first create all possible interactions between pixels whose pairwise \f$L_2\f$ distance is
    ## smaller than \f$\sqrt{2}\f$.
    nn = sg.NearestNeighbors(8,8)
    interactions = nn.getAllInteractions(3, 2**0.5)

    ## We then compare a standard sparse grid with a sparse grid learner that only contains the
    ## aforementioned interaction terms.
    standard_accuracy = evaluate(X_tr, y_tr, X_te, y_te)
    print("The standard sparse grid achieved an accuracy of {:2.3f}".format(standard_accuracy))
    ia_accuracy = evaluate(X_tr, y_tr, X_te, y_te, interactions)
    print("The interaction aware grid achieved an accuracy of {:2.3f}".format(ia_accuracy))

if __name__ == '__main__':
    main()
