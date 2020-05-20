// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// Shared pointers.
%shared_ptr(sgpp::datadriven::DMSystemMatrixBase)
%shared_ptr(sgpp::datadriven::DMSystemMatrix)
%shared_ptr(sgpp::datadriven::DMSystemMatrixDRE)
%shared_ptr(sgpp::datadriven::DensitySystemMatrix)
%shared_ptr(sgpp::datadriven::DensityDifferenceSystemMatrix)
%shared_ptr(sgpp::datadriven::OperationRegularizationDiagonal)

%{
#include <sgpp/solver/TypesSolver.hpp>
%}

// base class is not exported from the configuration
%warnfilter(401) sgpp::datadriven::LearnerSGDEConfiguration;
%warnfilter(401) sgpp::datadriven::SparseGridDensityEstimatorConfiguration;

// overloaded methods in DataMatrixDistributed and DataVectorDistributed might be shadowed
%warnfilter(509, 516) sgpp::datadriven::DataMatrixDistributed;
%warnfilter(509, 516) sgpp::datadriven::DataVectorDistributed;
%ignore sgpp::datadriven::DataMatrixDistributed::operator();
%ignore sgpp::datadriven::DataVectorDistributed::operator();

// The Good, i.e. without any modifications
#ifdef SG_DATADRIVEN
%include "datadriven/src/sgpp/datadriven/algorithm/test_dataset.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DMSystemMatrixDRE.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DensitySystemMatrix.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DensityDifferenceSystemMatrix.hpp"

%include "datadriven/src/sgpp/datadriven/tools/Dataset.hpp"

%include "datadriven/src/sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/DatabaseConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/DensityEstimationTypeParser.hpp"
//%include "datadriven/src/sgpp/datadriven/configuration/GeometryConfiguration.hpp"
//%include "datadriven/src/sgpp/datadriven/configuration/GeometryConfigurationParser.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/LearnerConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/MatrixDecompositionTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/ParallelConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/RegularizationConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/configuration/RegularizationTypeParser.hpp"

%include "datadriven/src/sgpp/datadriven/algorithm/RefinementMonitor.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp"

%include "datadriven/src/sgpp/datadriven/scalapack/BlacsProcessGrid.hpp"
%include "datadriven/src/sgpp/datadriven/scalapack/DataMatrixDistributed.hpp"
%include "datadriven/src/sgpp/datadriven/scalapack/DataVectorDistributed.hpp"

%ignore sgpp::datadriven::VectorHash::operator();
%include "datadriven/src/sgpp/datadriven/algorithm/CombiScheme.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp"
%ignore *::operator=;
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOffline.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineGE.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnline.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDE.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp"


%feature("flatnested");
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatObjectStore.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDatabase.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/GridFactory.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp"

#ifdef USE_GSL
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSEigen.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp"
%ignore sgpp::datadriven::DBMatOfflineLU::DBMatOfflineLU(DBMatOfflineLU &&);
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineLU.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEEigen.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMS_SMW.hpp"
#endif /* USE_GSL */

#ifdef __AVX__
%include "datadriven/src/sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp"
#endif

%include "OpFactory.i"

%include "datadriven/src/sgpp/datadriven/tools/TypesDatadriven.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerBase.hpp"
%include "datadriven/src/sgpp/datadriven/application/DensityEstimator.hpp"

%include "datadriven/src/sgpp/datadriven/application/KernelDensityEstimator.hpp"
%newobject sgpp::datadriven::KernelDensityEstimator::margToDimX(size_t idim);
%newobject sgpp::datadriven::KernelDensityEstimator::marginalize(size_t idim);
%include "datadriven/src/sgpp/datadriven/application/SparseGridDensityEstimator.hpp"
%newobject sgpp::datadriven::SparseGridDensityEstimator::margToDimX(size_t idim);
%newobject sgpp::datadriven::SparseGridDensityEstimator::marginalize(size_t idim);

%include "datadriven/src/sgpp/datadriven/application/LearnerSGDE.hpp"
%include "datadriven/src/sgpp/datadriven/application/RegressionLearner.hpp"
%include "datadriven/src/sgpp/datadriven/application/ClassificationLearner.hpp"
%include "datadriven/src/sgpp/datadriven/tools/NearestNeighbors.hpp"

%include "datadriven/src/sgpp/datadriven/application/LearnerSGD.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerSVM.hpp"
%include "datadriven/src/sgpp/datadriven/application/PrimalDualSVM.hpp"

%include "datadriven/src/sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/ClassificationRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp"

#ifdef USE_MPI
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/RefinementHandler.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp"
//%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/MPIRequestPool.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/PendingMPIRequest.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp"
#endif

%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp"
%ignore sgpp::datadriven::DataShufflingFunctor::operator();
%ignore sgpp::datadriven::DataShufflingFunctorRandom::operator();
%ignore sgpp::datadriven::DataShufflingFunctorSequential::operator();
%ignore sgpp::datadriven::DataShufflingFunctorCrossValidation::operator();
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/shuffling/DataSourceShufflingTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorCrossValidation.hpp"

%ignore *::operator=;

%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp"
#ifdef ZLIB
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp"
#endif /* ZLIB */
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformationConfig.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp"
%ignore sgpp::datadriven::DataTransformation::DataTransformation(DataTransformation &&);
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp"

%ignore sgpp::datadriven::DataSource::begin;
%ignore sgpp::datadriven::DataSource::end;
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp"

%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClassification.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityRatioEstimation.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationCG.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationOnOff.hpp"

%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/VisualizationParameters.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/VisualizerDummy.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/visualization/VisualizerClassification.hpp"

%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/Metric.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/MSE.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/Accuracy.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/NegativeLogLikelihood.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/Scorer.hpp"

%ignore sgpp::datadriven::SparseGridMiner::print;
%include "datadriven/src/sgpp/datadriven/datamining/base/SparseGridMiner.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/base/SparseGridMinerSplittingTwoDatasets.hpp"

%include "datadriven/src/sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/ScorerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/MinerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/DensityDifferenceEstimationMinerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/DensityRatioEstimationMinerFactory.hpp"

//TODO(lettrich) : parser not wrapable because of unwrapped JSON, but json parsing inside python is easy.
//%include "datadriven/src/sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp"
#endif

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point };

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationTest.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonal.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityConditionalKDE.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp"

%ignore sgpp::datadriven::HashGridPointCompare::operator();
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationLimitFunctionValueRange.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"

// --------------------------------------
// renaming ambiguous function declarations for python
%ignore base::datadriven::OperationRosenblattTransformation::doTransformation(base::DataVector* alpha, base::DataMatrix* points, base::DataMatrix* pointscdf, size_t dim_start);

//-     namespace datadriven ------------------------------------------
namespace datadriven {
/*%nodefaultdtor sgpp::datadriven::OperationTest;
class OperationTest
{
public:
  virtual double test(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes) = 0;
  virtual double testMSE(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& refValues) = 0;
  virtual double testWithCharacteristicNumber(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes, sgpp::base::DataVector& charaNumbers) = 0;
  virtual void calculateROCcurve(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes, sgpp::base::DataVector& thresholds, sgpp::base::DataMatrix& ROC_curve) = 0;
};

class OperationRegularizationDiagonal
{
public:
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) = 0;
  static const int HKMIX = 1;
  static const int H0HKLAPLACE = 2;
  static const int ISOTROPIC_PENALTY = 3;
  static const int ANISOTROPIC_PENALTY = 4;
};

class OperationRosenblattTransformation
{
public:
  virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* points, base::DataMatrix* pointscdf) = 0;
  virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* points, base::DataMatrix* pointscdf, size_t dim_start) = 0;
};

class OperationInverseRosenblattTransformation
{
public:
  virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf, base::DataMatrix* points) = 0;
  virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf, base::DataMatrix* points, size_t dim_start) = 0;
};
*/
}

//- end namespace datadriven ------------------------------------------
