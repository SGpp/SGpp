from pysgpp.extensions.datadriven.uq.dists import SGDEdist

dist = SGDEdist.byLearnerSGDE(samples, gridConfig, adaptConfig, solverConfig,
                              regularizationConfig, learnerConfig)
