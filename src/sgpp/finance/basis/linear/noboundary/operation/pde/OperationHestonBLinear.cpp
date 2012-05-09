/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonBLinear.hpp"

//#include "pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp"
//#include "pde/basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/DPhidPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/DPhidPhiUpBBLinear.hpp"
#include "pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace finance
{

OperationHestonBLinear::OperationHestonBLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationHestonBLinear::~OperationHestonBLinear()
{
}

void OperationHestonBLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	XPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<XPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonBLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * phi * phi
	XPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<XPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonBLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
//	// dphi * dphi
//	DPhidPhiUpBBLinear func(this->storage);
//	sg::base::sweep<DPhidPhiUpBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);

}

void OperationHestonBLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// Dphi * dphi
		sg::pde::DowndPhidPhiBBIterativeLinear myDown(this->storage);
		myDown(alpha, result, dim);

//	// dphi * dphi
//	DPhidPhiDownBBLinear func(this->storage);
//	sg::base::sweep<DPhidPhiDownBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);

}

}
}
