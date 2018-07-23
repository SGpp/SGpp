// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATADRIVEN_HPP
#define DATADRIVEN_HPP

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#ifdef USE_GSL
#include <sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#endif /* USE_GSL */

#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>

#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>

#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/datadriven/application/KernelDensityEstimator.hpp>
#include <sgpp/datadriven/application/SparseGridDensityEstimator.hpp>

#include <sgpp/datadriven/application/ClassificationLearner.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/application/LearnerSVM.hpp>
#include <sgpp/datadriven/application/PrimalDualSVM.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>

#ifdef USE_MPI
#include <sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIRequestPool.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/PendingMPIRequest.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RefinementHandler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp>
#endif /* USE_MPI */

#include <sgpp/datadriven/tools/NearestNeighbors.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonal.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationLimitFunctionValueRange.hpp>

#include <sgpp/datadriven/tools/TypesDatadriven.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/datadriven/configuration/BatchConfiguration.hpp>
#include <sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/configuration/DatabaseConfiguration.hpp>

#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

/* ************************
 * datamining
 * ************************/
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>

#include <sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/MinerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp>

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/configuration/DensityEstimationTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GridTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/RegularizationTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/SLESolverTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GeneralGridTypeParser.hpp>

#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformationConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataSourceShufflingTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp>

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceFileTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>

#include <sgpp/datadriven/datamining/modules/scoring/CrossValidation.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Accuracy.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/NegativeLogLikelihood.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/RandomShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerShufflingTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SequentialShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SplittingScorer.hpp>

#endif /* DATADRIVEN_HPP */
