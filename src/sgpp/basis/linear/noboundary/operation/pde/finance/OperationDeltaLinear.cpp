/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationDeltaLinear::OperationDeltaLinear(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLinear::~OperationDeltaLinear()
{
}

void OperationDeltaLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinear::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiUpBBLinear func(this->storage);
	sweep<detail::XPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLinear::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiDownBBLinear func(this->storage);
	sweep<detail::XPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
