/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"
#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/HierarchisationModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/DehierarchisationModLinear.hpp"

#include "algorithm/common/sweep.hpp"

#include "data/DataVector.hpp"

namespace sg
{
/**
 * Implements the hierarchisation on a sprase grid with mod linear base functions
 *
 * @param node_values the functions values in the node base
 *
 * @todo (heinecke, nice) Implement the hierarchisation on the sparse grid with mod linear base functions
 */
void OperationHierarchisationModLinear::doHierarchisation(DataVector& node_values)
{
	detail::HierarchisationModLinear func(this->storage);
	sweep<detail::HierarchisationModLinear> s(func, this->storage);

	// Execute hierarchisation in every dimension of the grid
	for (size_t i = 0; i < this->storage->dim(); i++)
	{
		s.sweep1D(node_values, node_values, i);
	}
}

/**
 * Implements the dehierarchisation on a sprase grid with mod linear base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 *
 * @todo (heinecke, nice) Implement the dehierarchisation on the sparse grid with mod linear base functions
 */
void OperationHierarchisationModLinear::doDehierarchisation(DataVector& alpha)
{
	detail::DehierarchisationModLinear func(this->storage);
	sweep<detail::DehierarchisationModLinear> s(func, this->storage);

	// Execute hierarchisation in every dimension of the grid
	for (size_t i = 0; i < this->storage->dim(); i++)
	{
		s.sweep1D(alpha, alpha, i);
	}
}

}
