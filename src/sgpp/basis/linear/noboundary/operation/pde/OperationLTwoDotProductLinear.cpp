/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

OperationLTwoDotProductLinear::OperationLTwoDotProductLinear(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLTwoDotProductLinear::~OperationLTwoDotProductLinear()
{
}

void OperationLTwoDotProductLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinear func(this->storage);
	sweep<PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLTwoDotProductLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinear func(this->storage);
	sweep<PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
