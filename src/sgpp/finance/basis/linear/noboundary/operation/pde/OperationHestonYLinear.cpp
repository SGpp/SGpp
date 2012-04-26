/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonYLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace finance
{

OperationHestonYLinear::OperationHestonYLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationHestonYLinear::~OperationHestonYLinear()
{
}

void OperationHestonYLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	XPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<XPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	XPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<XPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	sg::finance::SqXdPhidPhiUpBBLinear func(this->storage);
	sg::base::sweep<sg::finance::SqXdPhidPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonYLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	sg::finance::SqXdPhidPhiDownBBLinear func(this->storage);
	sg::base::sweep<sg::finance::SqXdPhidPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
