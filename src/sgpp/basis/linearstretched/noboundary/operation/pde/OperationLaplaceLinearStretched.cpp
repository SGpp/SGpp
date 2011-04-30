/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/operation/pde/OperationLaplaceLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/common/DowndPhidPhiBBIterativeLinearStretched.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{
namespace pde
{

OperationLaplaceLinearStretched::OperationLaplaceLinearStretched(GridStorage* storage) : UpDownOneOpDim(storage)
{
}

OperationLaplaceLinearStretched::~OperationLaplaceLinearStretched()
{
}

void OperationLaplaceLinearStretched::specialOP(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
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

void OperationLaplaceLinearStretched::up(DataVector& alpha, DataVector& result, size_t dim)
{
	PhiPhiUpBBLinearStretched func(this->storage);
	sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinearStretched::down(DataVector& alpha, DataVector& result, size_t dim)
{
	PhiPhiDownBBLinearStretched func(this->storage);
	sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinearStretched::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinearStretched myDown(this->storage);
	myDown(alpha, result, dim);
}

void OperationLaplaceLinearStretched::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
}

}
}
