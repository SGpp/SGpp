from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.learner import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder

import numpy as np
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder

# define random variables
builder = ParameterBuilder()
up = builder.defineUncertainParameters()

up.new().isCalled('y1').withUniformDistribution(-1, 1)
up.new().isCalled('y2').withUniformDistribution(-1, 1)
up.new().isCalled('y3').withUniformDistribution(-1, 1)

params = builder.andGetResult()

# define model function
def g(x, **kws):
    return np.prod([4 * xi * (1 - xi) for xi in x])

uqSetting = UQBuilder().withSimulation(g)\
                       .andGetResult()

# define learner
builder = SimulationLearnerBuilder()
builder.buildInterpolant()
builder.withGrid().withLevel(4)
builder.withSpecification().withParameters(params)

learner = builder.andGetResult()

# define sampling strategy
builder = ASGCSamplerBuilder()
builder.withLearner(learner)\
       .withSpecification().withParameters(params)
# stop policy
builder.withStopPolicy().withAdaptiveIterationLimit(0)
sampler = builder.andGetResult()

# ----------------------------------------------
# first run
while sampler.hasMoreSamples():
    samples = sampler.nextSamples()
    uqSetting.runSamples(samples)
    sampler.learnData(uqSetting)

# ----------------------------------------------
# build analysis
builder = ASGCAnalysisBuilder()
builder.withLearner(learner)\
       .withAnalyticEstimationStrategy()
