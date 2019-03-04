import os

import numpy as np

# calculates random alpha and u for the genz test functions and saves them for later use
# The difficulty of the integration is determined via ||alpha||_1 = sum alpha_i
# all u_i must be in [0,1]
# in this particular case we create alpha,u for the genz Oscillatory and Corner Peak function
# The difficulties are chosen as in "CUBA- a library for multidimensional numerical integration" 
# by T. Hahn as 6.0 for oscillatory and 2.2 for corner peak.
dim = 8
numSamples = 20

difficultyO = 6.0
difficultyCP = 2.2

alphaO = np.random.rand(dim, numSamples)
uO = np.random.rand(dim, numSamples)
alphaCP = np.random.rand(dim, numSamples)
uCP = np.random.rand(dim, numSamples)

for i in range(numSamples):
    alphaO[:, i] = alphaO[:, i] / np.sum(alphaO[:, i]) * difficultyO
    alphaCP[:, i] = alphaCP[:, i] / np.sum(alphaCP[:, i]) * difficultyCP

# pathO = '/home/rehmemk/git/SGpp/activeSubSpaces/results/genzOscillatory{}D'.format(dim)
# pathCP = '/home/rehmemk/git/SGpp/activeSubSpaces/results/genzCornerPeak{}D'.format(dim)
# np.savetxt(os.path.join(pathO, 'alpha.txt'), alphaO)
# np.savetxt(os.path.join(pathO, 'u.txt'), uO)
# np.savetxt(os.path.join(pathCP, 'alpha.txt'), alphaCP)
# np.savetxt(os.path.join(pathCP, 'u.txt'), uCP)
