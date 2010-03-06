 /*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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

#include "sgpp.hpp"

#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/HierarchisationLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/DehierarchisationLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include "basis/basis.hpp"
#include "data/DataVector.hpp"

namespace sg
{

void OperationHierarchisationLinearBoundary::doHierarchisation(DataVector& node_values)
{
	detail::HierarchisationLinearBoundary func(this->storage);
	sweep<detail::HierarchisationLinearBoundary> s(func, this->storage);

	// N D case
	if (this->storage->dim() > 1)
	{
		for (size_t i = 0; i < this->storage->dim(); i++)
		{
			s.sweep1D_Boundary(node_values, node_values, i);
		}
	}
	// 1 D case
	else
	{
		s.sweep1D(node_values, node_values, 0);
	}
}

void OperationHierarchisationLinearBoundary::doDehierarchisation(DataVector& alpha)
{
	detail::DehierarchisationLinearBoundary func(this->storage);
	sweep<detail::DehierarchisationLinearBoundary> s(func, this->storage);

	// N D case
	if (this->storage->dim() > 1)
	{
		for (size_t i = 0; i < this->storage->dim(); i++)
		{
			s.sweep1D_Boundary(alpha, alpha, i);
		}
	}
	// 1 D case
	else
	{
		s.sweep1D(alpha, alpha, 0);
	}
}

}
