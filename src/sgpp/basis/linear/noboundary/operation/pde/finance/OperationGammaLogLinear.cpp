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
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
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
	PhiPhiUpBBLinear func(this->storage);
	sweep<PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinear func(this->storage);
	sweep<PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	PhidPhiUpBBLinear func(this->storage);
	sweep<PhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	PhidPhiDownBBLinear func(this->storage);
	sweep<PhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinear func(this->storage);
	sweep<DPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinear func(this->storage);
	sweep<DPhiPhiDownBBLinear> s(func, this->storage);

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
}
