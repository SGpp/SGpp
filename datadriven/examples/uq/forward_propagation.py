from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.learner import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder

import numpy as np
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d

# define random variables
builder = ParameterBuilder()
up = builder.defineUncertainParameters()

up.new().isCalled('y1').withUniformDistribution(0, 1)
up.new().isCalled('y2').withUniformDistribution(0, 1)

params = builder.andGetResult()

# define UQ setting
builder = ASGCUQManagerBuilder()
builder.withParameters(params)\
       .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                              KnowledgeTypes.SQUARED])\
       .useInterpolation()

# define model function
def g(x, **kws):
    return np.sum(x)
#     return np.prod([4 * xi * (1 - xi) for xi in x])

builder.defineUQSetting().withSimulation(g)

samplerSpec = builder.defineSampler()
samplerSpec.withGrid().withLevel(3)  # .withPolynomialBase(2)
samplerSpec.withRefinement()\
           .withAdaptThreshold(1e-10)\
           .withAdaptPoints(5)\
           .withBalancing()\
           .refineMostPromisingNodes().withSquaredSurplusRanking()\
                                      .createAllChildrenOnRefinement()
#         ref.withBalancing()\
#            .addMostPromisingChildren().withLinearSurplusEstimationRanking()

samplerSpec.withStopPolicy().withAdaptiveIterationLimit(0)

uqManager = builder.andGetResult()

# ----------------------------------------------
# first run
while uqManager.hasMoreSamples():
    uqManager.runNextSamples()

# ----------------------------------------------
# build analysis
analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                .withAnalyticEstimationStrategy()\
                                .andGetResult()

print analysis.computeMoments()['data']

distKDE = analysis.estimateDensity(dtype="gaussianKDE")
config = {"grid_dim": 2,
          "grid_level": 8,
          "grid_type": "Linear",
          "refinement_numSteps": 0,
          "refinement_numPoints": 10,
          "regularization_type": "Laplace",
          "crossValidation_lambda": 5e-4,
          "crossValidation_enable": False,
          "crossValidation_kfold": 5,
          "crossValidation_silent": False}
distSGDE = analysis.estimateDensity(dtype="sgde", config=config)

import matplotlib.pyplot as plt
plt.figure()
plotDensity1d(distKDE, label="KDE")
plotDensity1d(distSGDE, label="SGDE")
samples = distSGDE.getSamples()
plt.scatter(samples, np.zeros(samples.shape[0]))
plt.legend()

x = np.linspace(0, 1, 100)
plt.figure()
plt.plot(x, [distKDE.cdf(xi) for xi in x], label="KDE")
plt.plot(x, [distSGDE.cdf(xi) for xi in x], label="SGDE")
plt.legend()

plt.show()

