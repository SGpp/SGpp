/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLogLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/PhidPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/PhidPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/common/DowndPhidPhiBBIterativeLinearStretched.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationGammaLogLinearStretched::OperationGammaLogLinearStretched(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLogLinearStretched::~OperationGammaLogLinearStretched()
{
}

void OperationGammaLogLinearStretched::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretched func(this->storage);
	sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretched func(this->storage);
	sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	PhidPhiUpBBLinearStretched func(this->storage);
	sweep<PhidPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	PhidPhiDownBBLinearStretched func(this->storage);
	sweep<PhidPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinearStretched func(this->storage);
	sweep<DPhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinearStretched func(this->storage);
	sweep<DPhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
}

void OperationGammaLogLinearStretched::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinearStretched myDown(this->storage);
	myDown(alpha, result, dim);
}

}
}
