/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modbspline/operation/datadriven/OperationTestModBspline.hpp"

#include "exception/operation_exception.hpp"

#include "data/DataVector.hpp"
using namespace sg::base;

namespace sg
{

double OperationTestModBspline::test(DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestModBspline::testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestModBspline::testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
