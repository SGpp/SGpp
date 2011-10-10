/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/basis.hpp"
#include "basis/modwavelet/operation/datadriven/OperationTestModWavelet.hpp"
#include "basis/modwavelet/ModifiedWaveletBasis.hpp"

#include "algorithm/datadriven/test_dataset.hpp"

#include "exception/operation_exception.hpp"

namespace sg
{
namespace datadriven
{

  double OperationTestModWavelet::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes)
  {
    base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
    return test_dataset(this->storage, base, alpha, data, classes);
  }

  double OperationTestModWavelet::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues)
  {
    base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
    return test_dataset_mse(this->storage, base, alpha, data, refValues);
  }

  double OperationTestModWavelet::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers)
  {
    base::ModifiedWaveletBasis<unsigned int, unsigned int> base;
    return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
  }

}
}
