#!/usr/bin/env python

## \page example_sparsePathExample_py Calculating the regularization path
##
## This example generates a regularization path for sparsity-inducing penalties.
## The output format is a comma seperated file.
## It accepts one argument that determines the desired
## regularization penalty.       
## The weight path calculation follows the discussion in
## \verbatim
## "Regularization paths for generalized linear models via coordinate descent."
## Friedman, Jerome, Trevor Hastie, and Rob Tibshirani, Journal of statistical software, 2010.
## \endverbatim

import numpy as np
import pysgpp as sg; sg.omp_set_num_threads(4)
import pandas as pd
import sklearn.datasets as data
from scipy.sparse.linalg import LinearOperator, svds
import sys

## This function generates the Friedman1 dataset.
def generate_friedman1(seed):
    (X,y) = data.make_friedman1(n_samples=10000, random_state=seed, noise=1.0)
    y = sg.DataVector(y)   
    X = sg.DataMatrix(X)
    return X, y

## This function calculates the design matrix for the linear model.
def get_Phi(X_train):
    def eval_op(x, op, size):
        result_vec = sg.DataVector(size)
        x = sg.DataVector(np.array(x).flatten())
        op.mult(x, result_vec)
        return result_vec.array().copy()

    def eval_op_transpose(x, op, size):
        result_vec = sg.DataVector(size)
        x = sg.DataVector(np.array(x).flatten())
        op.multTranspose(x, result_vec)
        return result_vec.array().copy()

    num_elem = X_train.array().shape[0]
    
    grid = sg.Grid.createModLinearGrid(10)
    gen = grid.getGenerator()
    gen.regular(2)
    
    op = sg.createOperationMultipleEval(grid, X_train)
    matvec = lambda x: eval_op(x, op, num_elem)
    rmatvec = lambda x: eval_op_transpose(x, op, grid.getSize())

    shape = (num_elem, grid.getSize())
    linop = LinearOperator(shape, matvec, rmatvec, dtype='float64')

    Phi = linop.matmat(np.matrix(np.identity(grid.getSize())))
    return Phi

## This function calculates the value of \f$ \lambda \f$ for which all weights are zero.
def get_max_lambda(Phi, y, num_rows, l1_ratio=1.0):
    max_prod = 0
    for i in range(0, Phi.shape[1]):
        a = np.asarray(Phi[:,i]).flatten()
        prod = np.inner(a, y)
        max_prod = max(max_prod, prod)
    max_lambda = max_prod / (l1_ratio * num_rows)
    return max_lambda

## This function calculates the weights for different lambdas.
## The result is printed to the standard output; it then resembles
## a comma seperated value file.
def calculate_weight_path(X, y, max_lambda,penalty, l1_ratio):
    epsilon=0.001
    num_lambdas=25
    estimator = make_estimator(penalty, max_lambda, l1_ratio) 
    estimator.train(X, y) 

    ## We create a grid such that
    ## \f$ \lambda \in \{ \varepsilon \cdot \lambda_\text{max} \text{ and } \lambda_\text{max} \}\f$.
    min_lambda = epsilon * max_lambda
    lambda_grid = np.logspace(np.log10(max_lambda), np.log10(min_lambda), num=num_lambdas)
    print("no_learner, lambda, weights")
    for i, lamb in enumerate(lambda_grid):
        last_weights = estimator.getWeights()
        estimator = make_estimator(penalty, l1_ratio, lamb,)
        estimator.setWeights(last_weights) # reuse old weights
        estimator.train(X, y)
        print("{}, {:2.6f}, {}".format(i, lamb, np.array2string(estimator.getWeights().array(), separator=',', max_line_width=float('inf'))))

## This function returns an estimator that uses the given penalty, l1_ratio
## and \f$\lambda\f$
def make_estimator(penalty, l1_ratio, lambda_reg):
    grid = sg.RegularGridConfiguration()
    grid.dim_ = 10
    grid.level_ = 2
    grid.type_ = sg.GridType_ModLinear

    adapt = sg.AdaptivityConfiguration()
    adapt.numRefinements_ = 0
    adapt.noPoints_ = 0

    solv = sg.SLESolverConfiguration()
    solv.maxIterations_ = 500
    solv.threshold_ = 10e-10
    solv.type_ = sg.SLESolverType_FISTA

    final_solv = solv

    regular = sg.RegularizationConfiguration()
    regular.type_ = penalty
    regular.exponentBase_ = 1.0
    regular.lambda_ = lambda_reg
    regular.l1_ratio_ = l1_ratio

    estimator = sg.RegressionLearner(grid, adapt, solv, final_solv,regular)
    return estimator

def main():
    penalties = {'lasso': sg.RegularizationType_Lasso,
                 'elasticNet': sg.RegularizationType_ElasticNet,
                 'groupLasso': sg.RegularizationType_GroupLasso}
    if (len(sys.argv) >= 2) and (sys.argv[1] not in penalties):
        print("Call this script by ./sparsePathExample.py <regularization method> <l1_ratio>")
        print("Acceptable regularization methods are:")
        print(list(penalties.keys()))
        return
    reg_method = penalties[sys.argv[1] if len(sys.argv) >= 2 else 'lasso']
    if len(sys.argv) > 2 and sys.argv[1]=='elasticNet':
        l1_ratio = float(sys.argv[2])
    else:
        l1_ratio = 1.0

    X, y = generate_friedman1(123456)    

    ## Create the design matrix, calculate \f$ \lambda_\text{max} \f$ and calculate the weight path.
    Phi = get_Phi(X)
    max_lambda = get_max_lambda(Phi, y.array(), X.array().shape[0], l1_ratio)
    calculate_weight_path(X, y, max_lambda, reg_method, l1_ratio)
    
if __name__ == '__main__':
    main()
