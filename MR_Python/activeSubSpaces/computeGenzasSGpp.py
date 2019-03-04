import os
import time
import as1DIntegral
import asMain
import numpy as np


# approximates integrals of genz Oscillary and Corner Peak functions for alpha and u loaded from files
# using the active subspace 1D integration
def calculateIntegralsasSGpp(model, numThreads=4):
    numSamples = 3 
    for i in range(numSamples):
        minPoints = 10
        maxPoints = 1000
        numSteps = 3
        saveFlag = 1
        numShadow1DPoints = 100
        numRefine = 10
        initialLevel = 2
        doResponse = 1
        doIntegral = 1
        gridType = 'nakbsplinemodified'
        degree = 3
        responseType = 'adaptive'
        asmType = 'adaptive'
        integralType = 'Spline'
        appSplineLevel = 0;appSplineDegree = 0;minDataPoints = 0;maxDataPoints = 0;numDataSteps = 0
        
        start = time.time()
        asMain.executeMain(model, 'asSGpp', numThreads, minPoints, maxPoints, numSteps, saveFlag,
                           numShadow1DPoints, numRefine, initialLevel, doResponse, doIntegral,
                           gridType, degree, responseType, asmType, integralType, appSplineLevel,
                           appSplineDegree, minDataPoints, maxDataPoints, numDataSteps, genzIndex=i)
        print("\n index {} done in {}s \n".format(i, time.time() - start))


##### Main ####
dim = 2
numThreads = 4

# model = 'genzOscillatory{}D'.format(dim)
model = 'genzCornerPeak{}D'.format(dim)
calculateIntegralsasSGpp(model, numThreads)

