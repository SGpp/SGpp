/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "basis/basis.hpp"
#include "basis/linear/boundary/linearboundaryBase.hpp"
#include "algorithm/datadriven/test_dataset.hpp"

#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"

#include "data/DataVector.hpp"

namespace sg
{
namespace datadriven
{

  double OperationTestLinearBoundary::test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes)
  {
    base::linearboundaryBase<unsigned int, unsigned int> base;
    return test_dataset(this->storage, base, alpha, data, classes);
  }

  double OperationTestLinearBoundary::testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues)
  {
    base::linearboundaryBase<unsigned int, unsigned int> base;
    return test_dataset_mse(this->storage, base, alpha, data, refValues);
  }

  double OperationTestLinearBoundary::testWithCharacteristicNumber(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers)
  {
    base::linearboundaryBase<unsigned int, unsigned int> base;
    return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
  }

}
}

