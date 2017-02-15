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
%include "datadriven/src/sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DensitySystemMatrix.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/ConvergenceMonitor.hpp"
#ifdef USE_GSL
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOffline.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnline.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatOnlineDE.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSEigen.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDMSChol.hpp"
%include "datadriven/src/sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp"
#endif

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
%include "datadriven/src/sgpp/datadriven/application/LearnerSGDE.hpp"
%newobject sgpp::datadriven::LearnerSGDE::margToDimX(size_t idim);

#ifdef USE_GSL
%include "datadriven/src/sgpp/datadriven/application/LearnerSGDEOnOff.hpp"
#endif

%include "datadriven/src/sgpp/datadriven/application/LearnerSGD.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerSVM.hpp"
%include "datadriven/src/sgpp/datadriven/application/PrimalDualSVM.hpp"

%include "datadriven/src/sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp"
%include "datadriven/src/sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp"
#endif

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point };

%include "datadriven/src/sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationTest.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityConditionalKDE.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp"

%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationLimitFunctionValueRange.hpp"

%include "datadriven/src/sgpp/datadriven/application/RegularizationConfiguration.hpp"

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
