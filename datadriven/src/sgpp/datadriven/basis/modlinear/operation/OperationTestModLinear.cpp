/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include <sgpp/datadriven/basis/modlinear/operation/OperationTestModLinear.hpp>
#include <sgpp/base/basis/modlinear/ModifiedLinearBasis.hpp>

#include <sgpp/datadriven/algorithm/test_dataset.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    double OperationTestModLinear::test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes) {
      base::ModifiedLinearBasis<unsigned int, unsigned int> base;
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    double OperationTestModLinear::testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues) {
      base::ModifiedLinearBasis<unsigned int, unsigned int> base;
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    double OperationTestModLinear::testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers) {
      base::ModifiedLinearBasis<unsigned int, unsigned int> base;
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestModLinear::calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) {
      base::ModifiedLinearBasis<unsigned int, unsigned int> base;
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}
