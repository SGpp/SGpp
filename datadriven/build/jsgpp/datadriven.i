// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%{
#include <sgpp/solver/TypesSolver.hpp>
%}

// base class is not exported from the configuration
%warnfilter(401) sgpp::datadriven::LearnerSGDEConfiguration;

// The Good, i.e. without any modifications
#ifdef SG_DATADRIVEN
%include "datadriven/src/sgpp/datadriven/algorithm/test_dataset.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DensitySystemMatrix.hpp"
%include "datadriven/src/sgpp/datadriven/application/RegularizationConfiguration.hpp"
%rename (getConstTargets) sgpp::datadriven::Dataset::getTargets() const;
%rename (getConstData) sgpp::datadriven::Dataset::getData() const;
%include "datadriven/src/sgpp/datadriven/tools/Dataset.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/ConvergenceMonitor.hpp"

%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp"
%ignore *::operator=;
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOffline.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineGE.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp"
%rename (getConstOfflineObject) sgpp::datadriven::DBMatOnline::getOfflineObject() const;
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnline.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDE.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp"

#ifdef USE_GSL
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSEigen.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp"
%ignore sgpp::datadriven::DBMatOfflineLU::DBMatOfflineLU(DBMatOfflineLU &&);
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOfflineLU.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDEEigen.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp"
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
%include "datadriven/src/sgpp/datadriven/application/GaussianKDE.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerSGDE.hpp"
%include "datadriven/src/sgpp/datadriven/application/RegressionLearner.hpp"
%include "datadriven/src/sgpp/datadriven/application/ClassificationLearner.hpp"
%include "datadriven/src/sgpp/datadriven/tools/NearestNeighbors.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerSGDEOnOff.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerSGD.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerSVM.hpp"
%include "datadriven/src/sgpp/datadriven/application/PrimalDualSVM.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp"
%ignore sgpp::datadriven::MPIMethods::waitForIncomingMessageType(*);
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/MPIRequestPool.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/PendingMPIRequest.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/RefinementHandler.hpp"
%include "datadriven/src/sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp"


%include "datadriven/src/sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp"

%ignore  sgpp::datadriven::SampleProvider::operator=(SampleProvider&&);
%rename(__assign__) sgpp::datadriven::SampleProvider::operator =;
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp"
%ignore  sgpp::datadriven::FileSampleDecorator::operator=(FileSampleDecorator&&);
%rename(__assign__) sgpp::datadriven::FileSampleDecorator::operator =;
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp"
#ifdef ZLIB
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp"
#endif /* ZLIB */
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp"


%ignore sgpp::datadriven::DataSource::begin;
%ignore sgpp::datadriven::DataSource::end;
%include "datadriven/src/sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp"
%ignore  sgpp::datadriven::FitterConfiguration::operator=(FitterConfiguration&&);
%rename(__assign__) sgpp::datadriven::FitterConfiguration::operator =;
%rename (getConstGridConfig) sgpp::datadriven::FitterConfiguration::getGridConfig() const;
%rename (getConstRefinementConfig) sgpp::datadriven::FitterConfiguration::getRefinementConfig() const;
%rename (getConstSolverRefineConfig) sgpp::datadriven::FitterConfiguration::getSolverRefineConfig() const;
%rename (getConstSolverFinalConfig) sgpp::datadriven::FitterConfiguration::getSolverFinalConfig() const;
%rename (getConstRegularizationConfig) sgpp::datadriven::FitterConfiguration::getRegularizationConfig() const;
%rename (getConstMultipleEvalConfig) sgpp::datadriven::FitterConfiguration::getMultipleEvalConfig() const;
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp"
%ignore  sgpp::datadriven::ModelFittingBase::operator=(ModelFittingBase&&);
%rename(__assign__) sgpp::datadriven::ModelFittingBase::operator =;
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp"

%ignore  sgpp::datadriven::Metric::operator=(Metric&&);
%rename(__assign__) sgpp::datadriven::Metric::operator =;
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/Metric.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/MSE.hpp"
%ignore  sgpp::datadriven::ShufflingFunctor::operator=(ShufflingFunctor&&);
%rename(__assign__) sgpp::datadriven::ShufflingFunctor::operator =;
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/SequentialShufflingFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/RandomShufflingFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/ScorerShufflingTypeParser.hpp"
%ignore  sgpp::datadriven::Scorer::operator=(Scorer&&);
%rename(__assign__) sgpp::datadriven::Scorer::operator =;
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/Scorer.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/SplittingScorer.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/modules/scoring/CrossValidation.hpp"

%ignore  sgpp::datadriven::SparseGridMiner::operator=(SparseGridMiner&&);
%include "datadriven/src/sgpp/datadriven/datamining/base/SparseGridMiner.hpp"

%include "datadriven/src/sgpp/datadriven/datamining/configuration/GridTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/configuration/RegularizationTypeParser.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/configuration/SLESolverTypeParser.hpp"

%include "datadriven/src/sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/ScorerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/MinerFactory.hpp"
%include "datadriven/src/sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp"


//TODO(lettrich): parser not wrapable because of unwrapped JSON
//%include "datadriven/src/sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp"


#endif

//%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

//%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
//%apply std::vector<double> *INPUT { std::vector<double>& point };

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationTest.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonal.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityConditionalKDE.hpp"



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
