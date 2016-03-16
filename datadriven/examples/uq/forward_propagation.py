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
    return np.prod([4 * xi * (1 - xi) for xi in x])

builder.defineUQSetting().withSimulation(g)

samplerSpec = builder.defineSampler()
samplerSpec.withGrid().withLevel(3).withPolynomialBase(2)
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

import matplotlib.pyplot as plt
for dtype, config in [("gaussianKDE", {}),
                      ("sgde", {"grid_type": "Linear",
                                "grid_level": 5,
                                "crossValidation_enable": True})]:
    U = analysis.estimateDensity(dtype=dtype, config=config)
    plt.figure()
    plotDensity1d(U)
    plt.title(dtype)

#     x = np.linspace(0, 1, 100)
#     y = [U.cdf(xi) for xi in x]
#     plt.figure()
#     plt.plot(x, y)
#     plt.title(dtype)

    plt.show()

