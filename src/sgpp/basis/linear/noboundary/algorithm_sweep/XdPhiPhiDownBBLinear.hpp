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

#ifndef XDPHIPHIDOWNBBLINEAR_HPP
#define XDPHIPHIDOWNBBLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * Implementation of sweep operator (): 1D Down for
 * Bilinearform \f$\int_{x} x \frac{\partial \phi(x)}{x} \phi(x) dx\f$
 */
class XdPhiPhiDownBBLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the GridStorage Object
	GridStorage* storage;
	/// Pointer to the bounding box Obejct
	BoundingBox* boundingBox;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	XdPhiPhiDownBBLinear(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~XdPhiPhiDownBBLinear();

	/**
	 * This operations performs the calculation of down in the direction of dimension <i>dim</i>
	 * on a grid with fix Dirichlet 0 boundary condition
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

protected:

	/**
	 * recursive function for the calculation of Down without Bounding Box support
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr);

	/**
	 * recursive function for the calculation of Down wit Bounding Box support
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 * @param q interval width in the current dimension <i>dim</i>
	 * @param t interval offset in current dimension <i>dim</i>
	 */
	void recBB(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr, double q, double t);
};

} // namespace detail

} // namespace sg

#endif /* XDPHIPHIDOWNBBLINEAR_HPP */
