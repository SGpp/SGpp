/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de

#include "basis/linear/noboundary/common/DowndPhidPhiBBIterativeLinear.hpp"
#include "grid/common/BoundingBox.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

DowndPhidPhiBBIterativeLinear::DowndPhidPhiBBIterativeLinear(GridStorage* storage) : storage(storage)
{
}

DowndPhidPhiBBIterativeLinear::~DowndPhidPhiBBIterativeLinear()
{
}

void DowndPhidPhiBBIterativeLinear::operator()(DataVector& alpha, DataVector& result, size_t dim)
{
	// Bounding Box handling
	BoundingBox* boundingBox = this->storage->getBoundingBox();
	double q = boundingBox->getIntervalWidth(dim);
	double Qqout = 1.0/q;

	// init the coefficients of the ansatz functions with boundary
	result.setAll(0.0);

	if (q != 1.0)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			//only affects the diagonal of the stiffness matrix
			result[i] = alpha[i]*(Qqout*pow(2.0, static_cast<int>(level+1)));
		}
	}
	else
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			//only affects the diagonal of the stiffness matrix
			result[i] = alpha[i]*pow(2.0, static_cast<int>(level+1));
		}
	}
}

}
}
