/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "basis/linear/operation/common/OperationHierarchisationLinear.hpp"
#include "basis/linear/algorithm_sweep/HierarchisationLinear.hpp"
#include "basis/linear/algorithm_sweep/DehierarchisationLinear.hpp"

#include "sgpp.hpp"

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
