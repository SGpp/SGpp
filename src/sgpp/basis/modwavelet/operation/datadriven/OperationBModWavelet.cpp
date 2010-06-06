/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/basis.hpp"
#include "basis/modwavelet/operation/datadriven/OperationBModWavelet.hpp"

#include "sgpp.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

#include "exception/operation_exception.hpp"

namespace sg
{

void OperationBModWavelet::mult(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SModWaveletBase> op;
	modified_wavelet_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

void OperationBModWavelet::multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SModWaveletBase> op;
	modified_wavelet_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
}

}
