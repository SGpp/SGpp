from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
import matplotlib.pyplot as plt


builder = ParameterBuilder()
up = builder.defineUncertainParameters()
up.new().isCalled('x').withTNormalDistribution(0, .1, -2, 2).hasValue(0.)
up.new().isCalled('y').withUniformDistribution(0, 1)
up.new().isCalled('z').withUniformDistribution(0, 1)  # .hasValue(0)
params = builder.andGetResult()

n = 10

samplers = {'naive': MCSampler.withNaiveSampleGenerator(params),
            'sobol': MCSampler.withSobolSampleGenerator(params),
            'scrambled sobol': MCSampler.withScrambledSobolSampleGenerator(params),
            'latin hypercube': MCSampler.withLatinHypercubeSampleGenerator(params, n)
#             'stratified': MCSampler.withStratifiedSampleGenerator(params, 10)
            }

for name, sampler in samplers.items():
    samples = sampler.nextSamples(n)
    samples.selectProbabilisticSpace().selectActiveElements()
    ndarray = samples.ndarray()
    fig = plt.figure()
    plt.scatter(ndarray[:, 0], ndarray[:, 1])
    plt.title(name)
    fig.show()
plt.show()
