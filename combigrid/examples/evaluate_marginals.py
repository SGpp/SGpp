from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler import MCSampler

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# define parameters
parameterBuilder = ParameterBuilder()
up = parameterBuilder.defineUncertainParameters()

up.new().isCalled('phi').withLognormalDistribution(0.15, 0.2, alpha=.01)
up.new().isCalled('e').withBetaDistribution(3, 3, 0, 2)
up.new().isCalled('K_L').withLognormalDistribution(1e-12, np.exp(-1), alpha=.01)

# design parameter
up.new().isCalled('Q_CO2').withLognormalDistribution(8.87, 0.2, alpha=.01).hasValue(8.87)
up.new().isCalled('W_z').withBetaDistribution(2, 2, 0, 30).hasValue(30)

params = parameterBuilder.andGetResult()

# load marginals and transformations as list
marginals = params.getDistributions()
transformations = params.getTransformations()

x = np.linspace(0, 1, 100)
for i, (dist, trans) in enumerate(zip(marginals, transformations)):
    # evaluate them on the unit interval
    y = [trans.vol() * dist.pdf(trans.unitToProbabilistic(xi)) for xi in x]

    # plot the result
    plt.figure()
    plt.plot(x, y)
    plt.title("%s" % dist)

plt.show()
