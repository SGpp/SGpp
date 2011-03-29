/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef PHIPHIDOWNBBLINEARSTRETCHED_HPP
#define PHIPHIDOWNBBLINEARSTRETCHED_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * Implementation of sweep operator (): 1D Down for
 * Bilinearform \f$\int_{x} \phi(x) \phi(x) dx\f$
 */
class PhiPhiDownBBLinearStretched
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the GridStorage Object
	GridStorage* storage;
	/// Pointer to the bounding box Obejct
	Stretching* stretching;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	PhiPhiDownBBLinearStretched(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~PhiPhiDownBBLinearStretched();

	/**
	 * This operations performs the calculation of down in the direction of dimension <i>dim</i>
	 * on a grid with Dirichlet 0 boundary conditions.
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

protected:

	/**
	 * recursive function for the calculation of Down with Stretching
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr);

};

} // namespace detail

} // namespace sg

#endif /* PHIPHIDOWNBBLINEARSTRETCHED_HPP */
