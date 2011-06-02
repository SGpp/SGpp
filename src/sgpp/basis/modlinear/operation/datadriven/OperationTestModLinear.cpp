/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"

#include "exception/operation_exception.hpp"

#include "data/DataVector.hpp"

namespace sg
{
namespace datadriven
{

double OperationTestModLinear::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes)
{
	modified_linear_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestModLinear::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues)
{
	modified_linear_base<unsigned int, unsigned int> base;
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestModLinear::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers)
{
	modified_linear_base<unsigned int, unsigned int> base;
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
}
