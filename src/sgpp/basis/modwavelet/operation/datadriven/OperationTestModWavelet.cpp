/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/basis.hpp"
#include "basis/modwavelet/operation/datadriven/OperationTestModWavelet.hpp"

#include "sgpp.hpp"

#include "data/DataVector.hpp"

#include "exception/operation_exception.hpp"
using namespace sg::base;

namespace sg
{

double OperationTestModWavelet::test(DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	modified_wavelet_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestModWavelet::testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	modified_wavelet_base<unsigned int, unsigned int> base;
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestModWavelet::testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	modified_wavelet_base<unsigned int, unsigned int> base;
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
