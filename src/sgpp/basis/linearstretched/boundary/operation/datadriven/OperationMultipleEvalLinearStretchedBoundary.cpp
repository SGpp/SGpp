/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/AlgorithmDGEMV.hpp"

#include "basis/basis.hpp"

#include "basis/linearstretched/boundary/operation/datadriven/OperationMultipleEvalLinearStretchedBoundary.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{
namespace base
{

void OperationMultipleEvalLinearStretchedBoundary::mult(DataVector& alpha, DataVector& result)
{
	AlgorithmDGEMV<SLinearStretchedBoundaryBase> op;
	linearstretchedboundaryBase<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, *(this->dataset_), result);
}

void OperationMultipleEvalLinearStretchedBoundary::multTranspose(DataVector& source, DataVector& result)
{
	AlgorithmDGEMV<SLinearStretchedBoundaryBase> op;
	linearstretchedboundaryBase<unsigned int, unsigned int> base;

	op.mult_transposed(storage, base, source, *(this->dataset_), result);
}

}
}
