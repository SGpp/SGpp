/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/basis/poly/operation/OperationTestPoly.hpp"

#include "datadriven/algorithm/test_dataset.hpp"

#include "base/exception/operation_exception.hpp"

namespace sg
{
namespace datadriven
{

double OperationTestPoly::test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes)
{
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestPoly::testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues)
{
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestPoly::testWithCharacteristicNumber(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers)
{
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
}
