# -------------------------------------------------------------------------------
# DataDist tests
# -------------------------------------------------------------------------------
import unittest
import os
import matplotlib.pyplot as plt
import numpy as np
import json

from matplotlib import rc
from scipy.integrate.quadpack import dblquad

from pysgpp import createOperationDensityMarginalize, \
    createOperationLTwoDotExplicit, createOperationQuadrature, \
    createOperationMakePositive, DataVector, Grid, \
    BandwidthOptimizationType_SILVERMANSRULE, \
    KernelType_GAUSSIAN
import pysgpp.extensions.datadriven.uq.dists as dists
from pysgpp.extensions.datadriven.uq.dists import J, Normal, Uniform, SGDEdist
from pysgpp.extensions.datadriven.uq.plot import plotDensity2d
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d, plotSG1d
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSGDE2d, plotSG2d
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotDensity3d, plotSG3d
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation
from pysgpp.extensions.datadriven.uq.dists.MultivariateNormal import MultivariateNormal
from pysgpp.extensions.datadriven.uq.dists.KDEDist import KDEDist
from pysgpp.extensions.datadriven.uq.quadrature.sparse_grid import doQuadrature
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist


class SGDEdistTest(unittest.TestCase):

    def testExpPoly2d(self):
        trainSamples = np.loadtxt("exp_2d.csv").T
        # build parameter set
        dist_sgde = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                                 config={"grid_level": 4,
                                                         "grid_type": "modpoly",
                                                         "grid_maxDegree": 6,
                                                         "refinement_numSteps": 0,
                                                         "refinement_numPoints": 10,
                                                         "solver_threshold": 1e-10,
                                                         "solver_verbose": True,
                                                         "regularization_type": "Laplace",
                                                         "crossValidation_lambda": 0.000562341,
                                                         "crossValidation_enable": False,
                                                         "crossValidation_kfold": 5,
                                                         "crossValidation_silent": True,
                                                         "sgde_makePositive": False,
                                                         "sgde_makePositive_candidateSearchAlgorithm": "joined",
                                                         "sgde_makePositive_interpolationAlgorithm": "setToZero",
                                                         "sgde_makePositive_verbose": True,
                                                         "sgde_unitIntegrand": True})

        # build parameter set
        dist_kde = dists.KDEDist(trainSamples,
                                 kernelType=KernelType_GAUSSIAN,
                                 bandwidthOptimizationType=BandwidthOptimizationType_SILVERMANSRULE)

        # fig = plt.figure()
        # plotSG2d(dist.grid, dist.alpha, show_grid_points=True)
        # plt.scatter(trainSamples[:, 0], trainSamples[:, 1], np.zeros(trainSamples.shape[0]))
        # plt.title("%.12f" % dist.vol)

        fig, _, _ = plotDensity3d(dist_sgde)
        plt.title("SGDE: vol=%g" % dist_sgde.vol)

        fig, _, _ = plotDensity3d(dist_kde)
        plt.title("KDE: vol=1.0")
        plt.show()

    def testExp2d(self):
        trainSamples = np.loadtxt("exp_2d.csv").T
        # build parameter set
        dist = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                            config={"grid_level": 7,
                                                    "grid_type": "linear",
                                                    "grid_maxDegree": 1,
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "solver_threshold": 1e-10,
                                                    "solver_verbose": False,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 0.000562341,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": False,
                                                    "sgde_makePositive": True,
                                                    "sgde_makePositive_candidateSearchAlgorithm": "joined",
                                                    "sgde_makePositive_interpolationAlgorithm": "interpolateBoundaries1d",
                                                    "sgde_makePositive_verbose": True,
                                                    "sgde_unitIntegrand": True})

        fig, ax, _ = plotDensity3d(dist)
        ax.scatter(trainSamples[:, 0], trainSamples[:, 1], np.zeros(trainSamples.shape[0]))
        ax.set_title("vol=%.12f" % dist.vol)
        fig.show()
        plt.show()

    def test1DNormalDist(self):
        # prepare data
        U = dists.TNormal(0.5, .2, -1, 2)
        np.random.seed(1234567)
        trainSamples = np.array([U.rvs(1000)]).T
        testSamples = np.array([U.rvs(1000)]).T

        # build parameter set
        dist = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                            config={"grid_level": 6,
                                                    "grid_type": "modlinear",
                                                    "grid_maxDegree": 3,
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "solver_threshold": 1e-10,
                                                    "solver_verbose": True,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_enable": True,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": False,
                                                    "sgde_makePositive": False,
                                                    "sgde_makePositive_candidateSearchAlgorithm": "fullGrid",
                                                    "sgde_makePositive_interpolationAlgorithm": "setToZero",
                                                    "sgde_makePositive_verbose": True,
                                                    "sgde_unitIntegrand": False},
                                            bounds=np.array([U.getBounds()]))

        fig = plt.figure()
        plotDensity1d(U, label="analytic")
        plotDensity1d(dist, label="sgde")
        plt.legend()
#         plt.title("mean = %g ~ %g (err=%g), var = %g ~ %g (err=%g)" % (np.mean(trainSamples),
#                                                                        dist.mean(),
#                                                                        np.abs(np.mean(trainSamples) - dist.mean()) / np.mean(trainSamples),
#                                                                        np.var(trainSamples),
#                                                                        dist.var(),
#                                                                        np.abs(np.var(trainSamples) - dist.var()) / np.var(trainSamples)
#                                                                        ))

        print("1d: mean = %g ~ %g (err=%g)" % (np.mean(trainSamples), dist.mean(), (np.abs(np.mean(trainSamples) - dist.mean()) / np.mean(trainSamples))))
        print("1d: var = %g ~ %g (err=%g)" % (np.var(trainSamples), dist.var(), (np.abs(np.var(trainSamples) - dist.var()) / np.var(trainSamples))))
        print("KL = %g" % U.klDivergence(dist, testSamples, testSamples))
        print("CE = %g" % dist.crossEntropy(testSamples))
        print("MSE = %g" % dist.l2error(U, testSamples, testSamples))
        plt.show()

    def test2DNormalDist(self):
        # prepare data
        U = dists.J([dists.Normal(2.0, .5, -1, 4),
                     dists.Normal(1.0, .5, -1, 3)])

        U = dists.J([dists.Normal(0.5, .5, -1, 2),
                     dists.Normal(0.5, .4, -1, 2)])

        np.random.seed(1234567)
        trainSamples = U.rvs(300)
        testSamples = U.rvs(1000)

        # build parameter set
        dist = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                            config={"grid_level": 5,
                                                    "grid_type": "modlinear",
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 0.000562341,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": False,
                                                    "sgde_makePositive": False,
                                                    "sgde_makePositive_candidateSearchAlgorithm": "joined",
                                                    "sgde_makePositive_interpolationAlgorithm": "setToZero",
                                                    "sgde_makePositive_generateConsistentGrid": False,
                                                    "sgde_makePositive_verbose": True,
                                                    "sgde_unitIntegrand": True},
                                            bounds=U.getBounds())
        fig = plt.figure()
        plotDensity2d(U)
        fig.show()

        fig = plt.figure()
        plotSG2d(dist.grid, dist.alpha, addContour=True,
                 show_negative=True, show_grid_points=True)
        fig.show()

        print("2d: mean = %g ~ %g" % (U.mean(), dist.mean()))
        print("2d: var = %g ~ %g" % (U.var(), dist.var()))
        plt.show()
        print("KL = %g" % U.klDivergence(dist, testSamples, testSamples))
        print("CE = %g" % dist.crossEntropy(testSamples))
        print("MSE = %g" % dist.l2error(U, testSamples, testSamples))

    def test2DNormalMoments(self):
        mean = 0
        var = 0.5

        U = dists.J([dists.Normal(mean, var, -2, 2),
                     dists.Normal(mean, var, -2, 2)])

        np.random.seed(1234567)
        trainSamples = U.rvs(1000)
        dist = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                            config={"grid_level": 5,
                                                    "grid_type": "linear",
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 0.000562341,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": True,
                                                    "sgde_makePositive": True},
                                            bounds=U.getBounds())
        samples_dist = dist.rvs(1000, shuffle=True)
        kde = KDEDist(trainSamples)
        samples_kde = kde.rvs(1000, shuffle=True)
        # -----------------------------------------------
        self.assertTrue(np.abs(U.mean() - dist.mean()) < 1e-2, "SGDE mean wrong")
        self.assertTrue(np.abs(U.var() - dist.var()) < 4e-2, "SGDE variance wrong")
        # -----------------------------------------------

        # print the results
        print("E(x) ~ %g ~ %g" % (kde.mean(), dist.mean()))
        print("V(x) ~ %g ~ %g" % (kde.var(), dist.var()))
        print("log  ~ %g ~ %g" % (kde.crossEntropy(trainSamples),
                                  dist.crossEntropy(trainSamples)))
        print("-" * 60)

        print(dist.cov())
        print(kde.cov())

        sgde_x1 = dist.marginalizeToDimX(0)
        kde_x1 = kde.marginalizeToDimX(0)

        plt.figure()
        plotDensity1d(U.getDistributions()[0], label="analytic")
        plotDensity1d(sgde_x1, label="sgde")
        plotDensity1d(kde_x1, label="kde")
        plt.title("mean: sgde=%g, kde=%g; var: sgde=%g, kde=%g" % (sgde_x1.mean(),
                                                                   kde_x1.mean(),
                                                                   sgde_x1.var(),
                                                                   kde_x1.var()))
        plt.legend()

        fig = plt.figure()
        plotDensity2d(U, addContour=True)
        plt.title("analytic")

        fig = plt.figure()
        plotDensity2d(kde, addContour=True)
        plt.scatter(samples_kde[:, 0], samples_kde[:, 1])
        plt.title("kde")

        fig = plt.figure()
        plotDensity2d(dist, addContour=True)
        plt.scatter(samples_dist[:, 0], samples_dist[:, 1])
        plt.title("sgde (I(f) = %g)" % (np.prod(U.getBounds()) * doQuadrature(dist.grid, dist.alpha),))

        plt.show()
#         print dist.toJson()
#         print SGDEdist.fromJson(json.loads(dist.toJson()))

    def test1DCDFandPPF(self):
        # prepare data
        U = Normal(0.5, 0.1, 0, 1)
        train_samples = U.rvs(1000).reshape(1000, 1)

        dist = SGDEdist.byLearnerSGDEConfig(train_samples,
                                            config={"grid_level": 5,
                                                    "grid_type": "poly",
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 0.000562341,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": True},
                                            bounds=U.getBounds())

        fig = plt.figure()
        plt.hist(train_samples, bins=10, normed=True)
        plotDensity1d(U)
        plotDensity1d(dist)
        plt.title("original space")
        fig.show()

        transformed_samples = dist.cdf(train_samples)

        fig = plt.figure()
        plt.hist(transformed_samples, bins=10, normed=True)
        plt.title("uniform space")
        fig.show()

        transformed_samples = dist.ppf(transformed_samples)

        fig = plt.figure()
        plt.hist(transformed_samples, bins=10, normed=True)
        plotDensity1d(U)
        plotDensity1d(dist)
        plt.title("original space")
        fig.show()
        plt.show()

    def test2DPPF(self):
        # prepare data
        C = np.array([[0.1, 0.08], [0.08, 0.1]]) / 10.
        U = dists.MultivariateNormal([0.5, 0.5], C, 0, 1)

        train_samples = U.rvs(1000)

        fig = plt.figure()
        plotDensity2d(U)
        plt.title('true density')
        fig.show()

        dist = SGDEdist.byLearnerSGDEConfig(train_samples,
                                            config={"grid_level": 5,
                                                    "grid_type": "linear",
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 0.000562341,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": True},
                                            bounds=U.getBounds())
        fig = plt.figure()
        plotDensity2d(dist)
        plt.title('estimated SGDE density')
        fig.show()

        samples = dists.J([dists.Uniform(0, 1),
                           dists.Uniform(0, 1)]).rvs(1000)

        fig = plt.figure()
        plt.plot(samples[:, 0], samples[:, 1], "o ")
        plt.title('uniformly drawn samples')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()

        transformed_samples = dist.ppf(samples)

        fig = plt.figure()
        plt.plot(transformed_samples[:, 0], transformed_samples[:, 1], "o ")
        plt.title('transformed samples')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()
        plt.show()

    def test2DCDFandPPF(self, plot=True):
        # prepare data
        C = np.array([[0.1, 0.08], [0.08, 0.1]]) / 10.
        U = dists.MultivariateNormal([0.5, 0.5], C, 0, 1)
        train_samples = U.rvs(1000)

        if plot:
            fig = plt.figure()
            plotDensity2d(U)
            plt.title('true density')
            fig.show()

        dist = SGDEdist.byLearnerSGDEConfig(train_samples,
                                            config={"grid_level": 5,
                                                    "grid_type": "polyClenshawCurtis",
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 0.000562341,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": True,
                                                    "sgde_makePositive": False},
                                            bounds=U.getBounds())

        if plot:
            fig = plt.figure()
            plotDensity2d(dist)
            plt.title('estimated SGDE density')
            fig.show()

        samples = dists.J([dists.Uniform(0, 1),
                           dists.Uniform(0, 1)]).rvs(500)

        if plot:
            fig = plt.figure()
            plt.plot(samples[:, 0], samples[:, 1], "o ")
            plt.title('u space')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            fig.show()
        else:
            print("-" * 80)
            print(samples)

        transformed_samples = dist.ppf(samples, shuffle=False)

        if plot:
            fig = plt.figure()
            plt.plot(transformed_samples[:, 0], transformed_samples[:, 1], "o ")
            plt.title('x space (transformed)')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            fig.show()
        else:
            print("-" * 80)
            print(transformed_samples)

        samples = dist.cdf(transformed_samples, shuffle=False)

        if plot:
            fig = plt.figure()
            plt.plot(samples[:, 0], samples[:, 1], "o ")
            plt.title('u space (transformed)')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            fig.show()

            plt.show()
        else:
            print("-" * 80)
            print(samples)

    def test2DCovarianceMatrix(self):
        # prepare data
        np.random.seed(1234567)
        C = np.array([[0.3, 0.09], [0.09, 0.3]]) / 10.

        U = dists.MultivariateNormal([0.5, 0.5], C, 0, 1)
        samples = U.rvs(2000)
        kde = KDEDist(samples)

        sgde = SGDEdist.byLearnerSGDEConfig(samples,
                                            bounds=U.getBounds(),
                                            config={"grid_level": 5,
                                                    "grid_type": "linear",
                                                    "grid_maxDegree": 1,
                                                    "refinement_numSteps": 0,
                                                    "refinement_numPoints": 10,
                                                    "solver_threshold": 1e-10,
                                                    "solver_verbose": False,
                                                    "regularization_type": "Laplace",
                                                    "crossValidation_lambda": 3.16228e-06,
                                                    "crossValidation_enable": False,
                                                    "crossValidation_kfold": 5,
                                                    "crossValidation_silent": False,
                                                    "sgde_makePositive": True,
                                                    "sgde_makePositive_candidateSearchAlgorithm": "joined",
                                                    "sgde_makePositive_interpolationAlgorithm": "setToZero",
                                                    "sgde_makePositive_verbose": True,
                                                    "sgde_generateConsistentGrid": True,
                                                    "sgde_unitIntegrand": True})

        sgde_x1 = sgde.marginalizeToDimX(0)
        sgde_x2 = sgde.marginalizeToDimX(1)

        plt.figure()
        plotDensity1d(sgde_x1, label="x1")
        plotDensity1d(sgde_x2, label="x2")
        plt.title("mean: x1=%g, x2=%g; var: x1=%g, x2=%g" % (sgde_x1.mean(),
                                                             sgde_x2.mean(),
                                                             sgde_x1.var(),
                                                             sgde_x2.var()))
        plt.legend()

        jsonStr = sgde.toJson()
        jsonObject = json.loads(jsonStr)
        sgde = Dist.fromJson(jsonObject)

        fig = plt.figure()
        plotDensity2d(U, addContour=True)
        plt.title("analytic")

        fig = plt.figure()
        plotDensity2d(kde, addContour=True)
        plt.title("kde")

        fig = plt.figure()
        plotDensity2d(sgde, addContour=True)
        plt.title("sgde (I(f) = %g)" % (doQuadrature(sgde.grid, sgde.alpha),))

        # print the results
        print("E(x) ~ %g ~ %g" % (kde.mean(), sgde.mean()))
        print("V(x) ~ %g ~ %g" % (kde.var(), sgde.var()))
        print("-" * 60)

        print(kde.cov())
        print(sgde.cov())

        self.assertTrue(np.linalg.norm(C - kde.cov()) < 1e-2, "KDE cov wrong")
        self.assertTrue(np.linalg.norm(np.corrcoef(samples.T) - kde.corrcoeff()) < 1e-1, "KDE corrcoef wrong")
        plt.show()

    def test_1DNormalDist_variance(self):
        # prepare data
        U = dists.Normal(1, 2, -8, 8)
#         U = dists.Normal(0.5, .2, 0, 1)

        # define linear transformation
        trans = JointTransformation()
        a, b = U.getBounds()
        trans.add(LinearTransformation(a, b))

        # get a sparse grid approximation
        grid = Grid.createPolyGrid(U.getDim(), 10)
        grid.getGenerator().regular(5)
        gs = grid.getStorage()

        # now refine adaptively 5 times
        p = DataVector(gs.getDimension())
        nodalValues = np.ndarray(gs.getSize())

        # set function values in alpha
        for i in range(gs.getSize()):
            gs.getPoint(i).getStandardCoordinates(p)
            nodalValues[i] = U.pdf(trans.unitToProbabilistic(p.array()))

        # hierarchize
        alpha = hierarchize(grid, nodalValues)
        dist = SGDEdist(grid, alpha, bounds=U.getBounds())

        fig = plt.figure()
        plotDensity1d(U,
                      alpha_value=0.1,
                      mean_label="$\mathbb{E}",
                      interval_label="$\alpha=0.1$")
        fig.show()

        fig = plt.figure()
        plotDensity1d(dist,
                      alpha_value=0.1,
                      mean_label="$\mathbb{E}",
                      interval_label="$\alpha=0.1$")
        fig.show()

        print("1d: mean = %g ~ %g" % (U.mean(), dist.mean()))
        print("1d: var = %g ~ %g" % (U.var(), dist.var()))
        plt.show()

    def test_2DNormalDist_variance(self):
        # prepare data
        U = dists.J([dists.Normal(2.0, .5, -1, 4),
                     dists.Normal(1.0, .5, -1, 3)])
#         U = dists.J([dists.Normal(0.5, .5, -1, 2),
#                      dists.Normal(0.5, .4, -1, 2)])

        # define linear transformation
        trans = JointTransformation()
        for a, b in U.getBounds():
            trans.add(LinearTransformation(a, b))

        # get a sparse grid approximation
        grid = Grid.createPolyGrid(U.getDim(), 10)
        grid.getGenerator().regular(5)
        gs = grid.getStorage()

        # now refine adaptively 5 times
        p = DataVector(gs.getDimension())
        nodalValues = np.ndarray(gs.getSize())

        # set function values in alpha
        for i in range(gs.getSize()):
            gs.getPoint(i).getStandardCoordinates(p)
            nodalValues[i] = U.pdf(trans.unitToProbabilistic(p.array()))

        # hierarchize
        alpha = hierarchize(grid, nodalValues)

#         # make positive
#         alpha_vec = DataVector(alpha)
#         createOperationMakePositive().makePositive(grid, alpha_vec)
#         alpha = alpha_vec.array()

        dist = SGDEdist(grid, alpha, bounds=U.getBounds())

        fig = plt.figure()
        plotDensity2d(U)
        fig.show()

        fig = plt.figure()
        plotSG2d(dist.grid, dist.alpha, addContour=True,
                 show_negative=True, show_grid_points=True)
        fig.show()

        print("2d: mean = %g ~ %g" % (U.mean(), dist.mean()))
        print("2d: var = %g ~ %g" % (U.var(), dist.var()))
        plt.show()


# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
    plt.show()
