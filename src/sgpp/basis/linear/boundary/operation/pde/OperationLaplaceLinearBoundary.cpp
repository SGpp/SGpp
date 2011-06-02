/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/common/DowndPhidPhiBBIterativeLinearBoundary.hpp"
#include "basis/linear/boundary/common/UpdPhidPhiBBIterativeLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include "grid/common/BoundingBox.hpp"

namespace sg
{
namespace pde
{

OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(sg::base::GridStorage* storage) : UpDownOneOpDim(storage)
{
}

OperationLaplaceLinearBoundary::~OperationLaplaceLinearBoundary()
{
}

void OperationLaplaceLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	PhiPhiUpBBLinearBoundary func(this->storage);
	sg::base::sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	PhiPhiDownBBLinearBoundary func(this->storage);
	sg::base::sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
	myDown(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
	myUp(alpha, result, dim);
}

}
}
