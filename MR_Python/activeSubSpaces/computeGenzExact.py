import os
import time
import as1DIntegral
import numpy as np


# calculates the exact integrals of genz Oscillary and Corner Peak functions for alpha and u loaded from files
def calculateIntegrals(alphas, us, model):
    [dim, numSamples] = np.shape(alphas)
    integrals = [0] * numSamples
    for i in range(numSamples):
        start = time.time()
        
        alpha = alphas[:, i]
        u = us[:, i]
        W1 = alpha / np.linalg.norm(alpha)
        if model == 'genzOscillatory{}D'.format(dim):
            g = lambda x: np.cos(2 * np.pi * u[0] + np.linalg.norm(alpha) * x)
        elif model == 'genzCornerPeak{}D'.format(dim):
            g = lambda x: (1 + np.linalg.norm(alpha) * x) ** (-dim - 1)
        integrals[i] = as1DIntegral.integrateASg(g, W1, dim)
        print(i)
        print(integrals[i])
        print('{}s'.format(time.time() - start))
        print(' ')
    return integrals


dim = 2  
print("starting execution")
pathO = '/home/rehmemk/git/SGpp/activeSubSpaces/results/genzOscillatory{}D'.format(dim)
pathCP = '/home/rehmemk/git/SGpp/activeSubSpaces/results/genzCornerPeak{}D'.format(dim)
alphaO = np.loadtxt(os.path.join(pathO, 'alpha.txt'))
uO = np.loadtxt(os.path.join(pathO, 'u.txt'))
alphaCP = np.loadtxt(os.path.join(pathCP, 'alpha.txt'))
uCP = np.loadtxt(os.path.join(pathCP, 'u.txt'))

integralsO = calculateIntegrals(alphaO, uO, model='genzOscillatory{}D'.format(dim))
integralsCP = calculateIntegrals(alphaCP, uCP, model='genzCornerPeak{}D'.format(dim))

# np.savetxt(os.path.join(pathO, 'integrals.txt'), integralsO)
# np.savetxt(os.path.join(pathCP, 'integrals.txt'), integralsCP)
