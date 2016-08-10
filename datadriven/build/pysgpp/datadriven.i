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
%include "datadriven/src/sgpp/datadriven/application/LearnerSGDE.hpp"
#endif

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point };

%include "datadriven/src/sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"

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

%rename("makePositive_cpp") makePositive(base::Grid*&, base::DataVector&, bool);
%ignore makePositive(base::Grid*&, base::DataVector&, bool);
%include "datadriven/src/sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp"
%rename("doLowerLimitation_cpp") doLowerLimitation(base::Grid*&, base::DataVector&, double, bool);
%rename("doUpperLimitation_cpp") doUpperLimitation(base::Grid*&, base::DataVector&, double, bool);
%rename("doLimitation_cpp") doLimitation(base::Grid*&, base::DataVector&, double, double);
%ignore doLowerLimitation(base::Grid*&, base::DataVector&, double, bool);
%ignore doUpperLimitation(base::Grid*&, base::DataVector&, double, bool);
%ignore doLimitation(base::Grid*&, base::DataVector&, double, double);
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

%extend sgpp::datadriven::OperationMakePositive {
  sgpp::base::Grid* makePositive(base::DataVector& newAlpha) {
    sgpp::base::Grid* grid = nullptr;
    $self->makePositive(grid, newAlpha);
    return grid;
  }
}

%extend sgpp::datadriven::OperationLimitFunctionValueRange {
  sgpp::base::Grid* doLowerLimitation(base::DataVector& newAlpha, double ylower, bool resetGrid) {
    sgpp::base::Grid* grid = nullptr;
    $self->doLowerLimitation(grid, newAlpha, ylower, resetGrid);
    return grid;
  }

  sgpp::base::Grid* doUpperLimitation(sgpp::base::DataVector& newAlpha, double yupper, bool resetGrid)  {
    sgpp::base::Grid* grid = nullptr;
    $self->doUpperLimitation(grid, newAlpha, yupper, resetGrid);
    return grid;
  }

  sgpp::base::Grid* doLimitation(sgpp::base::DataVector& newAlpha, double ylower, double yupper)  {
    sgpp::base::Grid* grid = nullptr;
    $self->doLimitation(grid, newAlpha, ylower, yupper);
    return grid;
  }
}
