from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder

from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d, \
    plotSobolIndices
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotFunction3d, plotSG3d
from pysgpp.extensions.datadriven.learner.Types import BorderTypes

import numpy as np
import matplotlib.pyplot as plt

# define random variables
builder = ParameterBuilder()
up = builder.defineUncertainParameters()

up.new().isCalled('x_0').withUniformDistribution(-2, 1)
up.new().isCalled('x_1').withUniformDistribution(0, 1)

params = builder.andGetResult()

# define UQ setting
builder = ASGCUQManagerBuilder()
builder.withParameters(params)\
       .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                              KnowledgeTypes.SQUARED])\
       .useInterpolation()

# define model function
def g(x, **kws):
    return np.arctan(50 * (x[0] - .35)) + np.pi / 2 + 4 * x[1] ** 3 + 10 * np.exp(x[0] * x[1] - 1)

builder.defineUQSetting().withSimulation(g)

samplerSpec = builder.defineSampler()
samplerSpec.withGrid().withLevel(3)\
                      .withPolynomialBase(2)\
                      .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
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

# ----------------------------------------------
fig, _, _ = plotSG3d(analysis.getGrid(),
                     analysis.getSurpluses())
fig.show()
# ----------------------------------------------
# show sobol indices

# anova decomposition
anova = analysis.getAnovaDecomposition(nk=len(params))

# main effects
me = anova.getSobolIndices()
te = anova.getTotalEffects()

names = anova.getSortedPermutations(me.keys())
values = [me[name] for name in names]
fig = plotSobolIndices(values, legend=True, names=names)
fig.show()

# ---------------------------------------------------------------------------
# distKDE = analysis.estimateDensity(dtype="gaussianKDE")
# config = {"grid_dim": 2,
#           "grid_level": 6,
#           "grid_type": "Linear",
#           "refinement_numSteps": 0,
#           "refinement_numPoints": 10,
#           "regularization_type": "Identity",
#           "crossValidation_lambda": 0.000562341,
#           "crossValidation_enable": False,
#           "crossValidation_kfold": 5,
#           "crossValidation_silent": False}
# distSGDE = analysis.estimateDensity(dtype="sgde", config=config)
#
# # ---------------------------------------------------------------------------
# print distSGDE.getBounds(), distSGDE.getSamples().shape, distSGDE.pdf(1.79753700978)
# y = analysis.eval(analysis.generateUnitSamples())
#
# plt.figure()
# plotDensity1d(distSGDE, label="SGDE")
# plotDensity1d(distKDE, label="KDE")
# samples = distSGDE.getSamples()
# plt.scatter(samples, np.zeros(samples.shape[0]))
# plt.hist(y, normed=True, cumulative=False, label="hist")
# plt.legend()
#
# x = np.linspace(distSGDE.getBounds()[0, 0],
#                 distSGDE.getBounds()[0, 1],
#                 100)
#
# plt.figure()
# plt.plot(x, [distSGDE.cdf(xi) for xi in x], label="SGDE")
# plt.plot(x, [distKDE.cdf(xi) for xi in x], label="KDE")
# plt.hist(y, normed=True, cumulative=True, label="hist")
# plt.legend()

plt.show()

