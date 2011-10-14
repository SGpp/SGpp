/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLELinearBoundary.hpp"


#include "basis/linear/boundary/common/DowndPhidPhiBBIterativeLinearBoundary.hpp"
#include "basis/linear/boundary/common/UpdPhidPhiBBIterativeLinearBoundary.hpp"


namespace sg
{
namespace finance
{

OperationLELinearBoundary::OperationLELinearBoundary(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage)
{
}

OperationLELinearBoundary::~OperationLELinearBoundary()
{
}

void OperationLELinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// Dphi * dphi
	sg::pde::UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
	myUp(alpha, result, dim);
}

void OperationLELinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// Dphi * dphi
	sg::pde::DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
	myDown(alpha, result, dim);
}

}
}
