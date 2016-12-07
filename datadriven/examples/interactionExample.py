#!/usr/bin/env python

################################################################
# This example tests the performance of interaction-term aware #
# sparse grids for the optical digits datasets.                #
# It calculates the testing accuracy of both a normal and a    #
# ia-term aware grid.                                          #
################################################################

import numpy as np
import pysgpp as sg; sg.omp_set_num_threads(4)
import pandas as pd
import sklearn.preprocessing as pre

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

def get_dataset():
    train_url = "https://archive.ics.uci.edu/ml/machine-learning-databases/optdigits/optdigits.tra"
    test_url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/optdigits/optdigits.tes'   
    print "Loading dataset from UCI repository."
    columns = ["x{}".format(i) for i in range(0, 64)] + ['digit']
    df_train = pd.read_csv(train_url, header=None, index_col=None)
    df_test = pd.read_csv(train_url, header=None, index_col=None)
    df_train.columns=columns
    df_test.columns=columns
    print "Preprocessing dataset."
    df_complete = df_train.append(df_test, ignore_index=True)
    scaler , _ = scale(df_complete)
    _, df_train = scale(df_train, scaler)
    _, df_test = scale(df_test, scaler)
    return df_train, df_test

def evaluate(X_tr, y_tr, X_te, y_te, interactions=None):
    grid = sg.RegularGridConfiguration()
    grid.dim_ = 64
    grid.level_ = 2
    grid.type_ = sg.GridType_ModLinear

    adapt = sg.AdpativityConfiguration()
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
    regular.regType_ = sg.RegularizationType_Identity
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

    nn = sg.NearestNeighbors(8,8)
    interactions = nn.getAllInteractions(3, 2**0.5)

    standard_accuracy = evaluate(X_tr, y_tr, X_te, y_te)
    print "The standard sparse grid achieved an accuracy of {:2.3f}".format(standard_accuracy)
    ia_accuracy = evaluate(X_tr, y_tr, X_te, y_te, interactions)
    print "The interaction aware grid achieved an accuracy of {:2.3f}".format(ia_accuracy)

if __name__ == '__main__':
    main()
