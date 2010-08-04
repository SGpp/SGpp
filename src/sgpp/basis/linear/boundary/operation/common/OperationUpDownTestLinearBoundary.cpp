/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/common/OperationUpDownTestLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/DPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/DPhidPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include "grid/common/BoundingBox.hpp"

namespace sg
{

OperationUpDownTestLinearBoundary::OperationUpDownTestLinearBoundary(GridStorage* storage) : storage(storage)
{
}

OperationUpDownTestLinearBoundary::~OperationUpDownTestLinearBoundary()
{
}

void OperationUpDownTestLinearBoundary::mult(DataVector& alpha, DataVector& result)
{
	this->updown(alpha, result);
}

void OperationUpDownTestLinearBoundary::updown(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());

	this->updown(alpha, beta, storage->dim() - 1);

	result.add(beta);
}

void OperationUpDownTestLinearBoundary::updown(DataVector& alpha, DataVector& result, size_t dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		up(alpha, temp, dim);
		updown(temp, result, dim-1);


		// Same from the other direction
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1);
		down(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		up(alpha, result, dim);

		DataVector temp(alpha.getSize());
		down(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationUpDownTestLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	//detail::PhiPhiUpBBLinearBoundary func(this->storage);
	//sweep<detail::PhiPhiUpBBLinearBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
	//detail::SqXdPhidPhiUpBBLinearBoundary func(this->storage);
	//sweep<detail::SqXdPhidPhiUpBBLinearBoundary> s(func, this->storage);

	// x * dphi * phi
	//detail::XdPhiPhiUpBBLinearBoundary func(this->storage);
	//sweep<detail::XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

	// x * phi * dphi
	//detail::XPhidPhiUpBBLinearBoundary func(this->storage);
	//sweep<detail::XPhidPhiUpBBLinearBoundary> s(func, this->storage);

	// dphi * dphi
	detail::DPhidPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::DPhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationUpDownTestLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	//detail::PhiPhiDownBBLinearBoundary func(this->storage);
	//sweep<detail::PhiPhiDownBBLinearBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
	//detail::SqXdPhidPhiDownBBLinearBoundary func(this->storage);
	//sweep<detail::SqXdPhidPhiDownBBLinearBoundary> s(func, this->storage);

	// x * dphi * phi
	//detail::XdPhiPhiDownBBLinearBoundary func(this->storage);
	//sweep<detail::XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

	// x * phi * dphi
	//detail::XPhidPhiDownBBLinearBoundary func(this->storage);
	//sweep<detail::XPhidPhiDownBBLinearBoundary> s(func, this->storage);

	// dphi * dphi
	//detail::DPhidPhiDownBBLinearBoundary func(this->storage);
	//sweep<detail::DPhidPhiDownBBLinearBoundary> s(func, this->storage);

	detail::DPhidPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::DPhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
