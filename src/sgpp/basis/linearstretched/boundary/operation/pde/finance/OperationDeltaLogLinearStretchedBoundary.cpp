/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLogLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationDeltaLogLinearStretchedBoundary::OperationDeltaLogLinearStretchedBoundary(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLogLinearStretchedBoundary::~OperationDeltaLogLinearStretchedBoundary()
{
}

void OperationDeltaLogLinearStretchedBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
	sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
	sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
