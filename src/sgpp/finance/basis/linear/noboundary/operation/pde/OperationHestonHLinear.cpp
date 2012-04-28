/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonHLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace finance
{

OperationHestonHLinear::OperationHestonHLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationHestonHLinear::~OperationHestonHLinear()
{
}

void OperationHestonHLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	XPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<XPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	XPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<XPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonHLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
