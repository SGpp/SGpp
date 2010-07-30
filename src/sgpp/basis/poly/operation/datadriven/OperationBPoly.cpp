/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/poly/operation/datadriven/OperationBPoly.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

#include "exception/operation_exception.hpp"

namespace sg
{

void OperationBPoly::mult(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SPolyBase> op;

	op.mult(storage, base, alpha, data, result);
}

void OperationBPoly::multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SPolyBase> op;

	op.mult_transpose(storage, base, alpha, data, result);
}

}
