/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLogLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhidPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhidPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/common/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/common/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationGammaLogLinearStretchedBoundary::OperationGammaLogLinearStretchedBoundary(GridStorage* storage, DataMatrix& coef) : UpDownTwoOpDims(storage, coef)
{
}

OperationGammaLogLinearStretchedBoundary::~OperationGammaLogLinearStretchedBoundary()
{
}

void OperationGammaLogLinearStretchedBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::upOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	PhidPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<PhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::downOpDimOne(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * dphi
	PhidPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<PhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
	myUp(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
	myDown(alpha, result, dim);
}

}
}
