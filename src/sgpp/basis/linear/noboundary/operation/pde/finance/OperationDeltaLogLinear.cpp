/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLogLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationDeltaLogLinear::OperationDeltaLogLinear(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLogLinear::~OperationDeltaLogLinear()
{
}

void OperationDeltaLogLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinear::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	detail::PhidPhiUpBBLinear func(this->storage);
	sweep<detail::PhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinear::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	detail::PhidPhiDownBBLinear func(this->storage);
	sweep<detail::PhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
