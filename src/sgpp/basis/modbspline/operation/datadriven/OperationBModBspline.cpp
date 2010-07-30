/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modbspline/operation/datadriven/OperationBModBspline.hpp"

#include "exception/operation_exception.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{

void OperationBModBspline::mult(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SModBsplineBase> op;

	op.mult(storage, base, alpha, data, result);
}

void OperationBModBspline::multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	AlgorithmDGEMV<SModBsplineBase> op;

	op.mult_transpose(storage, base, alpha, data, result);
}

}
