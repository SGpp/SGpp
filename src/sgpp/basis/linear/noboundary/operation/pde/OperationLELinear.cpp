/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/OperationLELinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/DPhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/DPhidPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLELinear::OperationLELinear(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLELinear::~OperationLELinear()
{
}

void OperationLELinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * dphi
	detail::DPhidPhiUpBBLinear func(this->storage);
	sweep<detail::DPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLELinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * dphi
	detail::DPhidPhiDownBBLinear func(this->storage);
	sweep<detail::DPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
