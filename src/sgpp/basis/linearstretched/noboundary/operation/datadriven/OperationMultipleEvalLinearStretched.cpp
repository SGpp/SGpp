/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/AlgorithmDGEMV.hpp"
#include "basis/basis.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalLinearStretched.hpp"
#include "exception/operation_exception.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"


namespace sg
{
namespace base
{

void OperationMultipleEvalLinearStretched::mult(DataVector& alpha, DataVector& result)
{
	AlgorithmDGEMV<SLinearStretchedBase> op;
	linearstretched_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, *(this->dataset_), result);
}

void OperationMultipleEvalLinearStretched::multTranspose(DataVector& source, DataVector& result)
{
	AlgorithmDGEMV<SLinearStretchedBase> op;
	linearstretched_base<unsigned int, unsigned int> base;

	op.mult_transposed(storage, base, source, *(this->dataset_), result);
}

}
}
