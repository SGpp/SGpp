from pysgpp import (RegularGridConfiguration, Linear, LinearBoundary,
                    AdpativityConfiguration,
                    SLESolverConfiguration,
                    RegularizationConfiguration, Laplace,
                    LearnerSGDEConfiguration,
                    DataVector, DataMatrix, Grid, LearnerSGDE)
from pysgpp.extensions.datadriven.uq.dists import SGDEdist, TNormal, J
from pysgpp import DataMatrix
import numpy as np

# -------------------- prepare data
U = J([TNormal(0.5, 0.06, 0, 1)])
np.random.seed(12345)
samples = DataMatrix(U.rvs(1000))
# ---------- using SGDE from SG++ ------------------------
gridConfig = RegularGridConfiguration()
gridConfig.dim_ = U.getDim()
gridConfig.level_ = 6
gridConfig.type_ = LinearBoundary

adaptConfig = AdpativityConfiguration()
adaptConfig.noPoints_ = 10
adaptConfig.numRefinements_ = 0

solverConfig = SLESolverConfiguration()
solverConfig.maxIterations_ = 100
solverConfig.eps_ = 1e-10
solverConfig.threshold_ = 1e-10

regularizationConfig = RegularizationConfiguration()
regularizationConfig.regType_ = Laplace

learnerConfig = LearnerSGDEConfiguration()
learnerConfig.doCrossValidation_ = False
learnerConfig.kfold_ = 0
learnerConfig.lambdaStart_ = 1e-1
learnerConfig.lambdaEnd_ = 1e-10
learnerConfig.lambdaSteps_ = 0
learnerConfig.logScale_ = True
learnerConfig.shuffle_ = False
learnerConfig.seed_ = 1234567
learnerConfig.silent_ = False

learner = LearnerSGDE(gridConfig, adaptConfig, solverConfig,
                      regularizationConfig, learnerConfig)
learner.initialize(samples)
grid = learner.getGrid()
print learner.getDim(), learner.getNsamples()
print learner.getGrid()
print learner.getSurpluses()
# dist = SGDEdist.byLearnerSGDE(samples, gridConfig, adaptConfig, solverConfig,
#                               regularizationConfig, learnerConfig)
