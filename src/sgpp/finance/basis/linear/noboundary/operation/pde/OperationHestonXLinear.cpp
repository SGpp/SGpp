/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonXLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace finance
{

OperationHestonXLinear::OperationHestonXLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationHestonXLinear::~OperationHestonXLinear()
{
}

void OperationHestonXLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	sg::finance::XPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<sg::finance::XPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	sg::finance::XPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<sg::finance::XPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	sg::finance::XdPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<sg::finance::XdPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonXLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	sg::finance::XdPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<sg::finance::XdPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
