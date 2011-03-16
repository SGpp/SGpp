/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "basis/linear/noboundary/common/DowndPhidPhiBBIterativeLinear.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::base;

namespace sg
{

OperationLaplaceLinear::OperationLaplaceLinear(GridStorage* storage) : UpDownOneOpDim(storage)
{
}

OperationLaplaceLinear::~OperationLaplaceLinear()
{
}

void OperationLaplaceLinear::specialOP(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
{
	// In direction gradient_dim we only calculate the norm of the gradient
	// The up-part is empty, thus omitted
	if(dim > 0)
	{
		DataVector temp(alpha.getSize());
		updown(alpha, temp, dim-1, gradient_dim);
		downOpDim(temp, result, gradient_dim);
	}
	else
	{
		// Terminates dimension recursion
		downOpDim(alpha, result, gradient_dim);
	}
}

void OperationLaplaceLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinear::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinear myDown(this->storage);
	myDown(alpha, result, dim);
}

void OperationLaplaceLinear::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
}

}
