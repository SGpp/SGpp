// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

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

%include "datadriven/src/sgpp/datadriven/DatadrivenOpFactory.hpp"
%include "datadriven/src/sgpp/datadriven/tools/TypesDatadriven.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerBase.hpp"
%include "datadriven/src/sgpp/datadriven/application/LearnerDensityCluster.hpp"
#endif

//%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

//%apply std::vector<std::pair<size_t, float_t> > *OUTPUT { std::vector<std::pair<size_t, float_t> >& result };
//%apply std::vector<float_t> *INPUT { std::vector<float_t>& point }; 

%include "datadriven/src/sgpp/datadriven/operation/hash/OperationTest.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/OperationRegularizationDiagonal.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/OperationRosenblattTransformation.hpp"
%include "datadriven/src/sgpp/datadriven/operation/hash/OperationInverseRosenblattTransformation.hpp"

//-     namespace datadriven ------------------------------------------
namespace datadriven {
/*%nodefaultdtor SGPP::datadriven::OperationTest;
class OperationTest
{
public:
  virtual float_t test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes) = 0;
  virtual float_t testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues) = 0;
  virtual float_t testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers) = 0;
  virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) = 0;
};

class OperationRegularizationDiagonal
{
public:
  virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;
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
