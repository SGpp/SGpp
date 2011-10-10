/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/basis.hpp"
#include "basis/prewavelet/PrewaveletBasis.hpp"
#include "algorithm/datadriven/test_dataset.hpp"

#include "basis/prewavelet/operation/datadriven/OperationTestPrewavelet.hpp"

#include "data/DataVector.hpp"

#include "data/DataMatrix.hpp"

namespace sg
{
namespace datadriven
{

  double OperationTestPrewavelet::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes)
  {
    base::PrewaveletBasis<unsigned int, unsigned int> base;
    return test_dataset(this->storage, base, alpha, data, classes);
  }

  double OperationTestPrewavelet::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues)
  {
    base::PrewaveletBasis<unsigned int, unsigned int> base;
    return test_dataset_mse(this->storage, base, alpha, data, refValues);
  }

  double OperationTestPrewavelet::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers)
  {
    base::PrewaveletBasis<unsigned int, unsigned int> base;
    return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
  }

}
}

