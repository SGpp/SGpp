/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modlinear/operation/datadriven/OperationBModLinear.hpp"

#include "exception/operation_exception.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{

void OperationBModLinear::mult(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SModLinearBase> op;
	modified_linear_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

void OperationBModLinear::multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SModLinearBase> op;
	modified_linear_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
}

}
