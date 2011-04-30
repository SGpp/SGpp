/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "sgpp.hpp"

#include "basis/basis.hpp"

#include "basis/linearstretched/noboundary/operation/datadriven/OperationTestLinearStretched.hpp"

#include "data/DataVector.hpp"

namespace sg
{
namespace base
{

double OperationTestLinearStretched::test(DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	linearstretched_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestLinearStretched::testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	linearstretched_base<unsigned int, unsigned int> base;
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestLinearStretched::testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	linearstretched_base<unsigned int, unsigned int> base;
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
}

