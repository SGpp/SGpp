/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef STDUPDOWN_HPP
#define STDUPDOWN_HPP

#include "grid/GridStorage.hpp"

#include "operation/common/OperationMatrix.hpp"

#include "data/DataVector.hpp"

#ifdef USEOMPTHREE
#include <omp.h>
#endif

namespace sg
{

/**
 * Implements a standard Up/Down Schema without any operation dim.
 *
 * @version $HEAD$
 */
class StdUpDown: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	StdUpDown(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~StdUpDown();


	virtual void mult(DataVector& alpha, DataVector& result);

protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;

#ifndef USEOMPTHREE
	/**
	 * Recursive procedure for updown
	 *
	 * @param dim the current dimension
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown(DataVector& alpha, DataVector& result, size_t dim);
#endif

#ifdef USEOMPTHREE
	/**
	 * Recursive procedure for updown, parallel version using OpenMP 3
	 *
	 * @param dim the current dimension
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown_parallel(DataVector& alpha, DataVector& result, size_t dim);
#endif

	/**
	 * 1D up Operation
	 *
	 * @param dim dimension in which to apply the up-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void up(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D down Operation
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void down(DataVector& alpha, DataVector& result, size_t dim) = 0;
};

}

#endif /* STDUPDOWN_HPP */
