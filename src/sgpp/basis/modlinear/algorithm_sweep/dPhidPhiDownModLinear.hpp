/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DPHIDPHIDOWNMODLINEAR_HPP
#define DPHIDPHIDOWNMODLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"
using namespace sg::base;

namespace sg
{

namespace detail
{

/**
 * Implementation of sweep operator (): 1D Down for
 * Bilinearform \f$\int_{x} \frac{\partial \phi(x)}{\partial x} \frac{\partial \phi(x)}{\partial x} dx\f$
 * on mod-linear grids
 */
class dPhidPhiDownModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	/// Pointer to GridStorage object
	GridStorage* storage;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	dPhidPhiDownModLinear(GridStorage* storage);

	/**
	 * Destructor
	 */
	~dPhidPhiDownModLinear();

	/**
	 * This operations performs the calculation of downGradient in the direction of dimension <i>dim</i>
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

protected:

	/**
	 * recursive function for the calculation of downGradient
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param f function value in the middle
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double f);
};

}

}

#endif /* DPHIDPHIDOWNMODLINEAR_HPP */
