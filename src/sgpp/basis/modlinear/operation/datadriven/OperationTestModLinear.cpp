/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"
#include "basis/modlinear/ModifiedLinearBasis.hpp"

#include "algorithm/datadriven/test_dataset.hpp"


namespace sg
{
namespace datadriven
{

  double OperationTestModLinear::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes)
  {
    base::ModifiedLinearBasis<unsigned int, unsigned int> base;
    return test_dataset(this->storage, base, alpha, data, classes);
  }

  double OperationTestModLinear::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues)
  {
    base::ModifiedLinearBasis<unsigned int, unsigned int> base;
    return test_dataset_mse(this->storage, base, alpha, data, refValues);
  }

  double OperationTestModLinear::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers)
  {
    base::ModifiedLinearBasis<unsigned int, unsigned int> base;
    return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
  }

}
}
