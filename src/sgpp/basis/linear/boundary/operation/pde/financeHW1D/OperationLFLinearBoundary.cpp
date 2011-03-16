/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLFLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::base;

namespace sg
{

OperationLFLinearBoundary::OperationLFLinearBoundary(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLFLinearBoundary::~OperationLFLinearBoundary()
{
}

void OperationLFLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// X * dphi * phi
	detail::XdPhiPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLFLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// X * dphi * phi
	detail::XdPhiPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
