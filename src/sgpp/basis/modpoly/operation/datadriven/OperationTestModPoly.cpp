/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modpoly/operation/datadriven/OperationTestModPoly.hpp"

#include "exception/operation_exception.hpp"

#include "data/DataVector.hpp"
using namespace sg::base;

namespace sg
{

double OperationTestModPoly::test(DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestModPoly::testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestModPoly::testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
