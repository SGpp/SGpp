/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"

#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"

#include "data/DataVector.hpp"

namespace sg
{
namespace datadriven
{

double OperationTestLinearBoundary::test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes)
{
	linearboundaryBase<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestLinearBoundary::testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues)
{
	linearboundaryBase<unsigned int, unsigned int> base;
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestLinearBoundary::testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers)
{
	linearboundaryBase<unsigned int, unsigned int> base;
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
}

