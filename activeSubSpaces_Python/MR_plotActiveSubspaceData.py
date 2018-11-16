from argparse import ArgumentParser
import os
import time

import activeSubspaceFunctions
import active_subspaces as ac
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pysgpp

#------------------------------------ main ---------------------------------------
# arguments
model = 'exp2D'
# gridType = 'nakbsplineextended'
# degree = 3
# maxPoints = 300

# load data 
objFunc = activeSubspaceFunctions.getFunction(model)
resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces_Python/results"
resultsPath = os.path.join(resultsPath, objFunc.getName())
# if gridType in ['OLS', 'QPHD', 'AS']:
#     folder = gridType + '_' + str(maxPoints)
# else: 
#     folder = gridType + '_' + str(degree) + '_' + str(maxPoints)
folders = ['AS_300', 'nakbsplineextended_3_300']

for folder in folders:        
    path = os.path.join(resultsPath, folder)
    with open(os.path.join(path, 'data.pkl'), 'rb') as fp:
        data = pickle.load(fp)
    
    # calculate error of first eigenvector
    eivecReference = objFunc.getEigenvec()
    eival = data["eigenvalues"]
    eivec = data["eigenvectors"]
    sampleRange = data["sampleRange"]
    gridType = data["gridType"]
    err = np.zeros(np.shape(eivec)[2])
    for i in range(len(err)):
        err[i] = np.linalg.norm(abs(eivec[:, :, i][0]) - abs(eivecReference[0]))
        
    if gridType in ['OLS', 'QPHD', 'AS']:
        label = gridType
    else:
        degree = data["degree"]
        label = gridType + '_' + str(degree)
    plt.loglog(sampleRange, err, label=label)
    
plt.legend()
plt.xlabel('# points')
plt.ylabel('error in 1st eigenvec')
plt.title(model)
plt.show()
