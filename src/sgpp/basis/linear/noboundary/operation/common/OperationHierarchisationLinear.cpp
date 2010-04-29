/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/linear/noboundary/operation/common/OperationHierarchisationLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/HierarchisationLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/DehierarchisationLinear.hpp"

#include "algorithm/common/sweep.hpp"

#include "basis/basis.hpp"
#include "data/DataVector.hpp"

namespace sg
{

void OperationHierarchisationLinear::doHierarchisation(DataVector& node_values)
{
	detail::HierarchisationLinear func(this->storage);
	sweep<detail::HierarchisationLinear> s(func, this->storage);

	// Execute hierarchisation in every dimension of the grid
	for (size_t i = 0; i < this->storage->dim(); i++)
	{
		s.sweep1D(node_values, node_values, i);
	}
}

void OperationHierarchisationLinear::doDehierarchisation(DataVector& alpha)
{
	detail::DehierarchisationLinear func(this->storage);
	sweep<detail::DehierarchisationLinear> s(func, this->storage);

	// Execute hierarchisation in every dimension of the grid
	for (size_t i = 0; i < this->storage->dim(); i++)
	{
		s.sweep1D(alpha, alpha, i);
	}
}

}
