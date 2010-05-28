/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modpoly/operation/datadriven/OperationBModPoly.hpp"

#include "exception/operation_exception.hpp"

#include "data/DataVector.hpp"

namespace sg
{

void OperationBModPoly::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmDGEMV<SModPolyBase> op;

	op.mult(storage, base, alpha, data, result);
}

void OperationBModPoly::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmDGEMV<SModPolyBase> op;

	op.mult_transpose(storage, base, alpha, data, result);
}

}
