/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/linear/boundary/operation/OperationUpDownTestLinearBoundary.hpp"

#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/PhidPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/PhidPhiUpBBLinearBoundary.hpp"

#include "base/algorithm/sweep.hpp"

#include "base/grid/common/BoundingBox.hpp"

namespace sg
{
namespace pde
{

OperationUpDownTestLinearBoundary::OperationUpDownTestLinearBoundary(sg::base::GridStorage* storage) : storage(storage)
{
}

OperationUpDownTestLinearBoundary::~OperationUpDownTestLinearBoundary()
{
}

void OperationUpDownTestLinearBoundary::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	this->updown(alpha, result);
}

void OperationUpDownTestLinearBoundary::updown(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector beta(result.getSize());

	this->updown(alpha, beta, storage->dim() - 1);

	result.add(beta);
}

void OperationUpDownTestLinearBoundary::updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		sg::base::DataVector temp(alpha.getSize());
		up(alpha, temp, dim);
		updown(temp, result, dim-1);


		// Same from the other direction
		sg::base::DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1);
		down(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		up(alpha, result, dim);

		sg::base::DataVector temp(alpha.getSize());
		down(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationUpDownTestLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	//PhiPhiUpBBLinearBoundary func(this->storage);
	//sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
	//SqXdPhidPhiUpBBLinearBoundary func(this->storage);
	//sweep<sg::finance::SqXdPhidPhiUpBBLinearBoundary> s(func, this->storage);

	// x * dphi * phi
	//XdPhiPhiUpBBLinearBoundary func(this->storage);
	//sweep<sg::finance::XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

	// x * phi * dphi
	sg::finance::XPhidPhiUpBBLinearBoundary func(this->storage);
	sg::base::sweep<sg::finance::XPhidPhiUpBBLinearBoundary> s(func, this->storage);

	// dphi * phi
	//DPhiPhiUpBBLinearBoundary func(this->storage);
	//sweep<sg::finance::DPhiPhiUpBBLinearBoundary> s(func, this->storage);

	// phi * dphi
	//PhidPhiUpBBLinearBoundary func(this->storage);
	//sweep<sg::finance::PhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationUpDownTestLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	//PhiPhiDownBBLinearBoundary func(this->storage);
	//sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);

	// x^2 * dphi * dphi
	//SqXdPhidPhiDownBBLinearBoundary func(this->storage);
	//sweep<sg::finance::SqXdPhidPhiDownBBLinearBoundary> s(func, this->storage);

	// x * dphi * phi
	//XdPhiPhiDownBBLinearBoundary func(this->storage);
	//sweep<sg::finance::XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

	// x * phi * dphi
	sg::finance::XPhidPhiDownBBLinearBoundary func(this->storage);
	sg::base::sweep<sg::finance::XPhidPhiDownBBLinearBoundary> s(func, this->storage);

	//  dphi * phi
	//DPhiPhiDownBBLinearBoundary func(this->storage);
	//sweep<sg::finance::DPhiPhiDownBBLinearBoundary> s(func, this->storage);

	//  phi * dphi
	//PhidPhiDownBBLinearBoundary func(this->storage);
	//sweep<sg::finance::PhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
}
