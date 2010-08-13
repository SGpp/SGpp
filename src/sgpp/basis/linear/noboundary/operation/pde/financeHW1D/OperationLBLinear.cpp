/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLBLinear::OperationLBLinear(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLBLinear::~OperationLBLinear()
{
}

void OperationLBLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	detail::DPhiPhiUpBBLinear func(this->storage);
	sweep<detail::DPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLBLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	detail::DPhiPhiDownBBLinear func(this->storage);
	sweep<detail::DPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
