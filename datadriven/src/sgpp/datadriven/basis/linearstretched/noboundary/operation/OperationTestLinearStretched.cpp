/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de), Dirk Pflueger (pflueged@in.tum.de)


#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
#include "datadriven/algorithm/test_dataset.hpp"

#include "datadriven/basis/linearstretched/noboundary/operation/OperationTestLinearStretched.hpp"


namespace sg {
  namespace datadriven {

    double OperationTestLinearStretched::test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      return test_dataset(this->storage, base, alpha, data, classes);
    }

    double OperationTestLinearStretched::testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      return test_dataset_mse(this->storage, base, alpha, data, refValues);
    }

    double OperationTestLinearStretched::testWithCharacteristicNumber(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers, 0.0);
    }

    void OperationTestLinearStretched::calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) {
      base::LinearStretchedBasis<unsigned int, unsigned int> base;
      test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
    }

  }
}

