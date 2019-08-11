try:
    from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
    from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
    from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder

    from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
    from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
    from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
    from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d, \
        plotSobolIndices
    from pysgpp.extensions.datadriven.uq.plot.plot3d import plotFunction3d, plotSG3d, \
        plotDensity3d
    from pysgpp.extensions.datadriven.learner.Types import BorderTypes

    import numpy as np
    import matplotlib.pyplot as plt
    from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
    import matplotlib
    from matplotlib import rc

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
        'size'   : 22}
matplotlib.rc('font', **font)
rc('text', usetex=True)

# define random variables
builder = ParameterBuilder()
up = builder.defineUncertainParameters()

up.new().isCalled('x_0').withTNormalDistribution(0.0, 0.5, -2, 1)
up.new().isCalled('x_1').withTNormalDistribution(0.5, 0.2, 0, 1).withInverseCDFTransformation()
up.new().isCalled('x_2,x_3,x_4,...').withSGDEDistribution(dist).withInverseCDFTransformation()

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
samplerSpec.withGrid().withLevel(4)\
                      .withPolynomialBase(2)\
                      .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
samplerSpec.withRefinement()\
           .withAdaptThreshold(1e-10)\
           .withAdaptPoints(5)\
           .withBalancing()\
           .refineMostPromisingNodes().withSquaredSurplusRanking()\
                                      .createAllChildrenOnRefinement()\
                                      .refineInnerNodes()
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

analysis.computeMoments()['data']

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

names = anova.getSortedPermutations(list(me.keys()))
values = [me[name] for name in names]
fig = plotSobolIndices(values, legend=True, names=names)
fig.show()

# ---------------------------------------------------------------------------
distKDE = analysis.estimateDensity(dtype="gaussianKDE")
config = {"grid_level": 6,
          "grid_type": "Linear",
          "refinement_numSteps": 0,
          "refinement_numPoints": 10,
          "regularization_type": "Laplace",
          "crossValidation_lambda": 3.16228e-06,
          "crossValidation_enable": True,
          "crossValidation_kfold": 5,
          "crossValidation_silent": False}
distSGDE = analysis.estimateDensity(dtype="sgde", config=config)

print("mean(u) = %g ~ %g (KDE) ~ %g (SGDE)" % (analysis.mean()[0], distKDE.mean(), distSGDE.mean()))
print("var(u) = %g ~ %g (KDE) ~ %g (SGDE)" % (analysis.var()[0], distKDE.var(), distSGDE.var()))
# ---------------------------------------------------------------------------
y = analysis.eval(analysis.generateUnitSamples())

fig = plt.figure()
plotDensity1d(distSGDE, label="SGDE", linewidth=3)
plt.vlines(distSGDE.mean(), 0, distSGDE.pdf(distSGDE.mean()),
           linewidth=3, color="red")
plt.vlines(distSGDE.mean() - distSGDE.std(), 0, distSGDE.pdf(distSGDE.mean() - distSGDE.std()),
           linewidth=3, color="red")
plt.vlines(distSGDE.mean() + distSGDE.std(), 0, distSGDE.pdf(distSGDE.mean() + distSGDE.std()),
           linewidth=3, color="red")
# plotDensity1d(distKDE, label="KDE")
# plt.vlines(distKDE.mean(), 0, distKDE.pdf(distKDE.mean()))
# samples = distKDE.getSamples()
# plt.scatter(samples, np.zeros(samples.shape[0]))
plt.hist(y, normed=True, cumulative=False, label="histogram")
plt.legend()
fig.show()

x = np.linspace(distSGDE.getBounds()[0, 0],
                distSGDE.getBounds()[0, 1],
                200)

fig = plt.figure()
plt.plot(x, [distSGDE.cdf(xi) for xi in x], label="SGDE",
         linewidth=3)
# plt.plot(x, [distKDE.cdf(xi) for xi in x], label="KDE")
plt.vlines(distSGDE.mean(), 0, distSGDE.cdf(distSGDE.mean()),
           linewidth=3, color="red")
plt.vlines(distSGDE.mean() - distSGDE.std(), 0, distSGDE.cdf(distSGDE.mean() - distSGDE.std()),
           linewidth=3, color="red")
plt.vlines(distSGDE.mean() + distSGDE.std(), 0, distSGDE.cdf(distSGDE.mean() + distSGDE.std()),
           linewidth=3, color="red")
plt.hist(y, normed=True, cumulative=True, label="histogram")
plt.legend()
fig.show()

plt.show()
