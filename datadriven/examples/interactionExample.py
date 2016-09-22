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
    df_tr = pd.read_csv('../../datasets/optdigits/optdigits_train.csv')
    X_tr = np.array(df_tr.ix[:,0:-1])
    y_tr = (df_tr.ix[:,-1]).values
    df_te = pd.read_csv('../../datasets/optdigits/optdigits_test.csv')
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
