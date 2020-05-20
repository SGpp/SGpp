# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

try:
    import sklearn.datasets as data
    import numpy as np
    import matplotlib.pyplot as plt
    from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
    from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d
    from pysgpp.extensions.datadriven.learner import Types
    import pysgpp as sg

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)


## This function generates the Friedman1 dataset on demand
def generate_friedman1(seed):
    (X,y) = data.make_friedman1(n_samples=2000, random_state=seed, noise=1.0)
    
    print(X.shape, y.shape)
    
    # transform values to DataMatrix/DataVector types
    X = sg.DataMatrix(X)
    y = sg.DataVector(y)
    
    return (X,y)


## This function defines a vector of uniformly random numbers
def randu_vec(size):
    return np.random.rand(size)


## This function defines a matrix of uniformly random numbers
def randu_mat(nrows, ncols):
    return np.random.rand(nrows, ncols)


## ~~~ Main ~~~ ##
def main():
    # Generate data
    print("generate dataset... ", end=' ')
    data_tr,_ = generate_friedman1(123456)
    print("Done")
    print("generated a friedman1 dataset (10D) with 2000 samples")
    
    # Config grid
    print("create grid config... ", end=' ')
    gridConfig = sg.RegularGridConfiguration()
    gridConfig.dim_ = 10
    gridConfig.level_ = 3
    gridConfig.type_ = sg.GridType_Linear
    print("Done")

    # Config adaptivity
    print("create adaptive refinement config... ", end=' ')
    adaptivityConfig = sg.AdaptivityConfiguration()
    adaptivityConfig.numRefinements_ = 0
    adaptivityConfig.numRefinementPoints_ = 10
    print("Done")
    
    # Config solver
    print("create solver config... ", end=' ')
    solv = sg.SLESolverConfiguration()
    solv.maxIterations_ = 1000
    solv.eps_ = 1e-14
    solv.threshold_ = 1e-14
    solv.type_ = sg.SLESolverType_CG
    print("Done")

    # Config regularization
    print("create regularization config... ", end=' ')
    regularizationConfig = sg.RegularizationConfiguration()
    regularizationConfig.regType_ = sg.RegularizationType_Laplace  
    print("Done")

    # Config cross validation for learner
    print("create learner config... ", end=' ')
    #crossValidationConfig = sg.CrossvalidationConfiguration()
    crossValidationConfig = sg.CrossvalidationConfiguration()
    crossValidationConfig.enable_ = False
    crossValidationConfig.kfold_ = 3
    crossValidationConfig.lambda_ = 3.16228e-06
    crossValidationConfig.lambdaStart_ = 1e-1
    crossValidationConfig.lambdaEnd_ = 1e-10
    crossValidationConfig.lambdaSteps_ = 3
    crossValidationConfig.logScale_ = True
    crossValidationConfig.shuffle_ = True
    crossValidationConfig.seed_ = 1234567
    crossValidationConfig.silent_ = False
    print("Done")

    #
    # Create the learner with the given configuration
    #
    print("create the learner... ")
    learner = sg.LearnerSGDE(gridConfig, adaptivityConfig, solv, regularizationConfig, crossValidationConfig)
    learner.initialize(data_tr)
    
    # Train the learner
    print("start training... ")
    learner.train()
    print("done training")
    
    #
    # Estimate the probability density function (pdf) via a Gaussian kernel 
    # density estimation (KDE) and print the corresponding values
    #
    kde = sg.KernelDensityEstimator(data_tr)
    x = sg.DataVector(learner.getDim())
    x.setAll(0.5)
    
    print("-----------------------------------------------")
    print(learner.getSurpluses().getSize(), " -> ", learner.getSurpluses().sum())
    print("pdf_SGDE(x) = ", learner.pdf(x), " ~ ", kde.pdf(x), " = pdf_KDE(x)")
    print("mean_SGDE = ", learner.mean(), " ~ ", kde.mean(), " = mean_KDE")
    print("var_SGDE = ", learner.variance(), " ~ ", kde.variance(), " = var_KDE")
    
    # Print the covariances
    C = sg.DataMatrix(gridConfig.dim_, gridConfig.dim_)
    print("----------------------- Cov_SGDE -----------------------")
    learner.cov(C)
    print(C)
    print("----------------------- Cov_KDE -----------------------")
    kde.cov(C)
    print(C)
    
    #
    # Apply the inverse Rosenblatt transformatio to a matrix of random points. To 
    # do this, first generate the random points uniformly, then initialize an 
    # inverse Rosenblatt transformation operation and apply it to the points.
    # Finally print the calculated values
    #
    print("-----------------------------------------------")
    opInvRos = sg.createOperationInverseRosenblattTransformation(learner.getGrid())
    points = sg.DataMatrix(randu_mat(12, gridConfig.dim_))
    print(points)
    
    pointsCdf = sg.DataMatrix(points.getNrows(), points.getNcols())
    opInvRos.doTransformation(learner.getSurpluses(), points, pointsCdf)
    
    #
    # To check whether the results are correct perform a Rosenform transformation on
    # the data that has been created by the inverse Rosenblatt transformation above
    # and print the calculated values
    #
    points.setAll(0.0)
    opRos = sg.createOperationRosenblattTransformation(learner.getGrid())
    opRos.doTransformation(learner.getSurpluses(), pointsCdf, points)
    
    print("-----------------------------------------------")
    print(pointsCdf)
    print("-----------------------------------------------")
    print(points)
    
    
if __name__ == '__main__':
    main()

