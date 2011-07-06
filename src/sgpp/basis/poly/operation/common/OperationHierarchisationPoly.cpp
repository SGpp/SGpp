/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

// #include "basis/poly/poly_base.hpp"
#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/poly/operation/common/OperationHierarchisationPoly.hpp"
#include "basis/poly/algorithm_sweep/HierarchisationPoly.hpp"
#include "basis/poly/algorithm_sweep/DehierarchisationPoly.hpp"

#include "data/DataVector.hpp"
#include "algorithm/common/sweep.hpp"

namespace sg
{
namespace base
{

void OperationHierarchisationPoly::doHierarchisation(DataVector& node_values)
{
	
	HierarchisationPoly func(this->storage, &this->base);
	sweep<HierarchisationPoly> s(func, this->storage);

	// Execute hierarchisation in every dimension of the grid
	for (size_t i = 0; i < this->storage->dim(); i++)
	{
		s.sweep1D(node_values, node_values, i);
	}
}

void OperationHierarchisationPoly::doDehierarchisation(DataVector& alpha)
{
	DehierarchisationPoly func(this->storage, &this->base);
	sweep<DehierarchisationPoly> s(func, this->storage);

	// Execute hierarchisation in every dimension of the grid
	for (size_t i = 0; i < this->storage->dim(); i++)
	{
		s.sweep1D(alpha, alpha, i);
	}
}

}
}
