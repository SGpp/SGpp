/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/AlgorithmDGEMV.hpp"

#include "basis/modlinear/ModifiedLinearBasis.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalModLinear.hpp"



namespace sg
{
namespace base
{

void OperationMultipleEvalModLinear::mult(DataVector& alpha, DataVector& result)
{
	AlgorithmDGEMV<SModLinearBase> op;
	ModifiedLinearBasis<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, *(this->dataset_), result);
}

void OperationMultipleEvalModLinear::multTranspose(DataVector& source, DataVector& result)
{
	AlgorithmDGEMV<SModLinearBase> op;
	ModifiedLinearBasis<unsigned int, unsigned int> base;

	op.mult_transposed(storage, base, source, *(this->dataset_), result);
}

}
}
