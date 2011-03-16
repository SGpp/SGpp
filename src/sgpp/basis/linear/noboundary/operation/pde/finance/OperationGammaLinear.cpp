/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>
using namespace sg::base;

namespace sg
{

OperationGammaLinear::OperationGammaLinear(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLinear::~OperationGammaLinear()
{
}

void OperationGammaLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiUpBBLinear func(this->storage);
	sweep<detail::XPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiDownBBLinear func(this->storage);
	sweep<detail::XPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiUpBBLinear func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiDownBBLinear func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiUpBBLinear func(this->storage);
	sweep<detail::SqXdPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiDownBBLinear func(this->storage);
	sweep<detail::SqXdPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
