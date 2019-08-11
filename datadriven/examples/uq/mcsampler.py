# coding: utf-8
try:
    from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
    from pysgpp.extensions.datadriven.uq.sampler import MCSampler

    from pysgpp import *

    import numpy as np
    import matplotlib.pyplot as plt
    from pysgpp.extensions.datadriven.uq.dists.Lognormal import Lognormal

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)


parameterBuilder = ParameterBuilder()
up = parameterBuilder.defineUncertainParameters()

up.new().isCalled('phi').withLognormalDistribution(0.15, 0.2, alpha=.01)
up.new().isCalled('e').withBetaDistribution(3, 3, 0, 2)
up.new().isCalled('K_L').withLognormalDistribution(1e-12, np.exp(-1), alpha=.01)

# design parameter
up.new().isCalled('Q_CO2').withLognormalDistribution(8.87, 0.2, alpha=.01)  # .hasValue(8.87)
up.new().isCalled('W_z').withBetaDistribution(2, 2, 0, 30)  # .hasValue(30)

params = parameterBuilder.andGetResult()

n = 9

a = MCSampler(params).nextSamples(n).ndarray()
b = MCSampler(params).nextSamples(n + 1).ndarray()


print(a)
print("-" * 80)
print(b)

print(a[:n, :] - b[:n, :])
