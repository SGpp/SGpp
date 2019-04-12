# -------------------------------------------------------------------------------
# DataDist tests
# -------------------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
    import numpy as np
    import json
    from scipy.integrate.quadpack import dblquad
    from pysgpp import createOperationDensityMarginalize, \
        createOperationLTwoDotExplicit, createOperationQuadrature, \
        createOperationMakePositive, DataVector, Grid, \
        BandwidthOptimizationType_SILVERMANSRULE, \
        KernelType_GAUSSIAN
    import pysgpp.extensions.datadriven.uq.dists as dists
    from pysgpp.extensions.datadriven.uq.dists import J, Normal, Uniform, SGDEdist, Lognormal
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

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)

def test_sgdeLaplace():
    l2_samples = 10000
    # sample_range = np.arange(10, 500, 50)
    sample_range = [10, 20, 50, 100, 200, 500]
    points = {}
    grids = ["linear",
             "modlinear", # keine OperationQuadrature
             "poly",
             "modpoly",
             "polyBoundary",
             "polyClenshawCurtis",
             "modPolyClenshawCurtis",
             "polyClenshawCurtisBoundary",
             "bsplineClenshawCurtis",
             "modBsplineClenshawCurtis" # keine OperationMultipleEval
    ]

    U = dists.J([dists.Lognormal.by_alpha(0.5, 0.1, 0.001),
                 dists.Lognormal.by_alpha(0.5, 0.1, 0.001)])
    l2_errors = {}
    for grid in grids:
        l2_errors[grid] = []
        points[grid] = []

    l2_errors["kde"] = []
    samples = 1000
    for samples in sample_range:
    # for lvl in range(5, 6):
        trainSamples = U.rvs(samples)
        # testSamples = U.rvs(l2_samples)
        for grid_name in grids:
            # build parameter set
            print("--------------------Samples: {} Grid: {}--------------------".format(samples, grid_name))
            dist_sgde = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                                     bounds=U.getBounds(),
                                                     unitIntegrand=True,
                                                     config={"grid_level": 1,
                                                             "grid_type": grid_name,
                                                             "grid_maxDegree": 6,
                                                             "refinement_numSteps": 0,
                                                             "refinement_numPoints": 10,
                                                             "solver_threshold": 1e-10,
                                                             "solver_verbose": False,
                                                             "regularization_type": "Laplace",
                                                             "crossValidation_lambda": 1e-6,
                                                             "crossValidation_enable": True,
                                                             "crossValidation_kfold": 4,
                                                             "crossValidation_lambdaSteps": 10,
                                                             "crossValidation_silent": False})
            points[grid_name].append(dist_sgde.grid.getSize())
            # l2_errors[grid_name].append(dist_sgde.l2error(U, testSamplesUnit=testSamples))
            l2_errors[grid_name].append(dist_sgde.l2error(U, n=l2_samples))
            # plt.figure()
            # plotDensity2d(U, levels=(10, 20, 40, 50, 60))
            # plt.figure()
            # plotDensity2d(dist_sgde, levels=(10, 20, 40, 50, 60))
            # plt.show()

        dist_kde = dists.KDEDist(trainSamples,
                                 kernelType=KernelType_GAUSSIAN,
                                 bandwidthOptimizationType=BandwidthOptimizationType_SILVERMANSRULE)
        l2_errors["kde"].append(dist_kde.l2error(U, testSamplesUnit=testSamples))

    for grid_name in grids:
        plt.plot(sample_range, l2_errors[grid_name], label=grid_name)
        # plt.plot(points[grid], l2_errors[grid_name],".-", label=grid_name)

    plt.plot(sample_range, l2_errors["kde"], label="KDE")

    # plt.plot([x for x in range(1,300, 100)], [l2_errors["kde"][0] for i in range(1,4)], label="KDE")

    plt.xlabel("# Gitterpunkte")
    plt.ylabel("L2-Fehler")
    plt.yscale("log")
    plt.legend()
    plt.show()

test_sgdeLaplace()
