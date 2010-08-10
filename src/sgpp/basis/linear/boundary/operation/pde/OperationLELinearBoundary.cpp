/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/OperationLELinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/DPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/DPhidPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLELinearBoundary::OperationLELinearBoundary(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLELinearBoundary::~OperationLELinearBoundary()
{
}

void OperationLELinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * dphi
	detail::DPhidPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::DPhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLELinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * dphi
	detail::DPhidPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::DPhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
