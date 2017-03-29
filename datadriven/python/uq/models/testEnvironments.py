import numpy as np

from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.helper import findSetBits, sortPermutations

from work.probabilistic_transformations_for_inference.sampling import TensorQuadratureSampleGenerationStrategy, \
    ApproximateFeketeSampleGeneratorStrategy, LejaSampleGeneratorStrategy
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder

class TestEnvironmentSG(object):

    def buildSetting(self,
                     params,
                     f=None,
                     level=None,
                     gridType=None,
                     grid=None,
                     deg=1,
                     maxGridSize=1000,
                     isFull=False,
                     boundaryLevel=None,
                     balancing=True,
                     epsilon=1e-15,
                     adaptive=None,
                     refinementTechnique="refinement",
                     adaptPoints=3,
                     saveAfterN=-1,
                     adaptRate=None,
                     uqSetting=None,
                     uqSettingRef=None,
                     knowledgeFilename=None,
                     knowledgeTypes=[KnowledgeTypes.SIMPLE,
                                     KnowledgeTypes.SQUARED,
                                     KnowledgeTypes.EXPECTATIONVALUE],
                     qoi=None,
                     toi=None):
        builder = ASGCUQManagerBuilder()

        builder.withParameters(params)\
               .withTypesOfKnowledge(knowledgeTypes)\
               .useInterpolation()

        if qoi is not None:
            builder.withQoI(qoi)
        if toi is not None:
            builder.withTimeStepsOfInterest(toi)
        if uqSetting is not None:
            builder.useUQSetting(uqSetting)
        else:
            builder.defineUQSetting().withSimulation(f).saveAfterEachRun(saveAfterN)

        if uqSettingRef is not None and len(uqSettingRef) > 0:
            builder.withTestSet(uqSettingRef)\
                   .learnWithTest()

        if knowledgeFilename is not None:
            builder.withKnowledge(knowledgeFilename)

        samplerSpec = builder.defineSampler()
        gridSpec = samplerSpec.withGrid()
        if grid is not None:
            gridSpec.fromGrid(grid)
        else:
            assert level is not None
            assert gridType is not None
            gridSpec.withLevel(level).hasType(gridType)
            if deg > 1:
                gridSpec.withDegree(deg)
            if isFull:
                gridSpec.isFull()
            if boundaryLevel is not None:
                gridSpec.withBoundaryLevel(boundaryLevel)

        if adaptive is not None:
            # specify the refinement
            refinement = samplerSpec.withRefinement()
            
            refinement.withAdaptThreshold(epsilon)\
                      .withAdaptPoints(adaptPoints)

            if balancing:
                refinement.withBalancing()
                
            if adaptRate is not None:
                refinement.withAdaptRate(adaptRate)

            if refinementTechnique == "refinement":
                refineNodes = refinement.refineMostPromisingNodes()
                refineNodes.createAllChildrenOnRefinement()
                if adaptive == "simple":
                    refineNodes.withSurplusRanking()
                elif adaptive == "weighted":
                    refineNodes.withWeightedSurplusRanking()
                elif adaptive == "l2":
                    refineNodes.withWeightedL2OptimizationRanking()
                elif adaptive == "exp":
                    refineNodes.withExpectationValueOptimizationRanking()
                elif adaptive == "var":
                    refineNodes.withVarianceOptimizationRanking()
                elif adaptive == "mean_squared":
                    refineNodes.withMeanSquaredOptRanking()
                elif adaptive == "squared":
                    refineNodes.withSquaredSurplusRanking()
                elif adaptive == "anchored_exp":
                    refineNodes.withAnchoredExpectationValueOptimizationRanking()
                elif adaptive == "anchored_var":
                    refineNodes.withAnchoredVarianceOptimizationRanking()
                elif adaptive == "anchored_l2":
                    refineNodes.withAnchoredWeightedL2OptimizationRanking()
                elif adaptive == "anchored_mean_squared":
                    refineNodes.withAnchoredMeanSquaredOptRanking()
                else:
                    raise AttributeError("unknown ranking method: refinement, %s" % adaptive)
            else:
                addNodes = refinement.addMostPromisingChildren()
                if adaptive == "weighted":
                    addNodes.withWeightedSurplusOptimizationRanking()
                elif adaptive == "l2":
                    addNodes.withWeightedL2OptimizationRanking()
                else:
                    raise AttributeError("unknown ranking method: predictive, %s" % adaptive)
                    
            if toi is not None:
                refinement.withAverageWeightening()

            samplerSpec.withStopPolicy().withGridSizeLimit(maxGridSize)

        uqManager = builder.andGetResult()

        # update the stats, which are not stored with the knowledge
        if knowledgeFilename is not None:
            uqManager.recomputeStats()

        return uqManager

    def runSampler(self, uqManager, label, alabel, blabel):
        uqSetting = self.uqSettings[label]
        # ----------------------------------------------
        # prepare folders
        pathResults = os.path.join(self.pathResults, alabel, blabel)

        if not os.path.exists(pathResults):
            os.makedirs(pathResults)

        for newdir in [os.path.join(pathResults, 'checkpoints'),
                       os.path.join(pathResults, 'grids'),
                       os.path.join(pathResults, 'samples')]:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        # write the setting to file
        uqManager.uqSetting.writeToFile()


class TestEnvironmentMC(object):

    def buildSetting(self,
                     params,
                     simulation,
                     N):
        mcSampler = MCSampler.withLatinHypercubeSampleGenerator(params, N)
        uqSetting = UQBuilder().withSimulation(simulation).andGetResult()
        return mcSampler, uqSetting


class ProbabilisticSpaceSGpp(object):
    
    def __init__(self, numDims):
        self.numDims = numDims
    
    def uniform(self, a=0, b=1):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        for idim in xrange(self.numDims):
            up.new().isCalled("x%i" % idim).withUniformDistribution(a, b)
        return up.andGetResult()

    def normal(self, mu=0, sigma=1, alpha=0.001):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        for idim in xrange(self.numDims):
            up.new().isCalled("x%i" % idim).withNormalDistribution(mu, sigma, 0.001)
        return builder.andGetResult()

    def multivariate_normal(self, mu=0, cov=None, a=0, b=1):
            # set distributions of the input parameters
        mu = np.array([mu] * numDims)
        if cov is None:
            # use standard values
            diag = np.eye(numDims) * 0.005
            offdiag = np.abs(np.eye(numDims) - 1) * 0.001
            cov = diag + offdiag
        # estimate the density
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        names = ", ".join(["x%i" for i in xrange(numDims)])
        up.new().isCalled(names).withMultivariateNormalDistribution(mu, cov, 0, 1)
        return builder.andGetResult()


class PCEBuilderHeat(object):

    def __init__(self, numDims):
        self.numDims = numDims


    def define_expansion(self, pce, expansion, degree_1d):
        if expansion == "full_tensor":
            pce.define_full_tensor_expansion(degree_1d)
        elif expansion == "total_degree":
            pce.define_isotropic_expansion(degree_1d, 1.0)
        else:
            raise AttributeError("expansion '%s' is not supported" % expansion)

    def define_full_tensor_samples(self, sample_type, rv_trans, expansion):
        return TensorQuadratureSampleGenerationStrategy("uniform", rv_trans, expansion)
    
    def define_approximate_fekete_samples(self, samples, pce, rv_trans):
        return ApproximateFeketeSampleGeneratorStrategy(samples, pce, rv_trans)
    
    def define_approximate_leja_samples(self, samples, pce, rv_trans):
        return LejaSampleGeneratorStrategy(samples, pce, rv_trans)

    def eval_samples(self, samples, rv_trans, f):
        trans_samples = rv_trans.map_from_canonical_distributions(samples)
        values = np.ndarray(trans_samples.shape[1])
        for i, sample in enumerate(trans_samples.T):
            values[i] = f(sample)
        return trans_samples, values

    def getSortedSobolIndices(self, pce):
        sobol_indices = pce.sobol_indices()
        indices = [findSetBits(i + 1) for i in xrange(len(sobol_indices))]
        indices, ixs = sortPermutations(indices, index_return=True)
        sobol_indices_dict = {}
        for index, i in zip(indices, ixs):
            sobol_indices_dict[index] = sobol_indices[i][0]
        return sobol_indices_dict
