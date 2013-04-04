/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "datadriven/basis/modwavelet/operation/OperationTestModWavelet.hpp"

#include "datadriven/algorithm/test_dataset.hpp"


namespace sg {
  namespace datadriven {

    double OperationTestModWavelet::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes) {
      base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    double OperationTestModWavelet::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues) {
      base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    double OperationTestModWavelet::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers) {
      base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestModWavelet::calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) {
      base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}
