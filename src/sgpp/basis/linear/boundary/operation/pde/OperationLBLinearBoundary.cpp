/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/OperationLBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLBLinearBoundary::OperationLBLinearBoundary(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLBLinearBoundary::~OperationLBLinearBoundary()
{
}

void OperationLBLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	detail::DPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::DPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLBLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	detail::DPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::DPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
