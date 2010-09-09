/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/poly/operation/datadriven/OperationTestPoly.hpp"

#include "data/DataVector.hpp"

#include "exception/operation_exception.hpp"

namespace sg
{

double OperationTestPoly::test(DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestPoly::testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestPoly::testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
