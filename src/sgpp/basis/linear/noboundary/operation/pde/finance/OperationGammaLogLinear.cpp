/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLogLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/common/DowndPhidPhiBBIterativeLinear.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>

namespace sg
{

OperationGammaLogLinear::OperationGammaLogLinear(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLogLinear::~OperationGammaLogLinear()
{
}

void OperationGammaLogLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	detail::PhidPhiUpBBLinear func(this->storage);
	sweep<detail::PhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	detail::PhidPhiDownBBLinear func(this->storage);
	sweep<detail::PhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	detail::DPhiPhiUpBBLinear func(this->storage);
	sweep<detail::DPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	detail::DPhiPhiDownBBLinear func(this->storage);
	sweep<detail::DPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
}

void OperationGammaLogLinear::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinear myDown(this->storage);
	myDown(alpha, result, dim);
}

}
