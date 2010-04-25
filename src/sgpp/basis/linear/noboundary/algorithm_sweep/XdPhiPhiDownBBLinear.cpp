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

#include "basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"

namespace sg
{

namespace detail
{

XdPhiPhiDownBBLinear::XdPhiPhiDownBBLinear(GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox())
{
}

XdPhiPhiDownBBLinear::~XdPhiPhiDownBBLinear()
{
}

void XdPhiPhiDownBBLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	double q = boundingBox->getIntervalWidth(dim);
	double t = boundingBox->getIntervalOffset(dim);

	bool useBB = false;

	if (q != 1.0 || t != 0.0)
	{
		useBB = true;
	}

	if (useBB)
	{
		recBB(source, result, index, dim, 0.0, 0.0, q, t);
	}
	else
	{
		rec(source, result, index, dim, 0.0, 0.0);
	}
}

void XdPhiPhiDownBBLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);

	double hhalf = 1.0/pow(2.0,static_cast<int>(l+1));
	double i_dbl = static_cast<double>(i);

	// integration
	result[seq] = (  ( (fl * ((hhalf * i_dbl) - hhalf)) + (fr * (((-1.0)*(hhalf * i_dbl)) - hhalf)) ) - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l)))) * alpha_value)) );    // diagonal entry

	// dehierarchisation
	double fm = ((fl+fr)/2.0) + alpha_value;

	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fl, fm);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fm, fr);
		}

		index.up(dim);
	}
}

void XdPhiPhiDownBBLinear::recBB(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr, double q, double t)
{
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);

	double hhalf = 1.0/pow(2.0,static_cast<int>(l+1));
	double i_dbl = static_cast<double>(i);

	// integration
	result[seq] = (  ( (fl * ((q*((hhalf * i_dbl) - hhalf))+(0.5*t))) + (fr * ((q*(((-1.0)*(hhalf * i_dbl)) - hhalf))-(0.5*t))) )
						  - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l))) * q) * alpha_value)) );    // diagonal entry

	// dehierarchisation
	double fm = ((fl+fr)/2.0) + alpha_value;

	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			recBB(source, result, index, dim, fl, fm, q, t);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			recBB(source, result, index, dim, fm, fr, q, t);
		}

		index.up(dim);
	}
}

} // namespace detail

} // namespace sg
