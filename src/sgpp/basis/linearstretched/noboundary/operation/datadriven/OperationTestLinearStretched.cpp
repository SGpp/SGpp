/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de), Dirk Pflueger (pflueged@in.tum.de)


#include "basis/basis.hpp"
#include "basis/linearstretched/noboundary/linearstretched_base.hpp"
#include "algorithm/datadriven/test_dataset.hpp"

#include "basis/linearstretched/noboundary/operation/datadriven/OperationTestLinearStretched.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{
namespace datadriven
{

double OperationTestLinearStretched::test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes)
{
  base::linearstretched_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestLinearStretched::testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues)
{
  base::linearstretched_base<unsigned int, unsigned int> base;
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestLinearStretched::testWithCharacteristicNumber(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers)
{
  base::linearstretched_base<unsigned int, unsigned int> base;
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
}

