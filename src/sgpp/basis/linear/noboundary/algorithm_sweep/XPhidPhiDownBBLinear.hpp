/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef XPHIDPHIDOWNBBLINEAR_HPP
#define XPHIDPHIDOWNBBLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * down-operation in dimension dim. for use with sweep
 */
class XPhidPhiDownBBLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the GridStorage Object
	GridStorage* storage;
	/// Pointer to the bounding box Obejct
	BoundingBox* boundingBox;
	/// width of the interval in dimension
	double q;
	/// intervals offset in dimension
	double t;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	XPhidPhiDownBBLinear(GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()), q(1.0), t(0.0)
	{
	}

	/**
	 * Destructor
	 */
	~XPhidPhiDownBBLinear()
	{
	}

	/**
	 * This operations performs the calculation of down in the direction of dimension <i>dim</i>
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * On level zero the getfixDirechletBoundaries of the storage object evaluated
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		q = boundingBox->getIntervalWidth(dim);
		t = boundingBox->getIntervalOffset(dim);

		bool useBB = false;

		if (q != 1.0 || t != 0.0)
		{
			useBB = true;
		}

		if (useBB)
		{
			recBB(source, result, index, dim, 0.0, 0.0);
		}
		else
		{
			rec(source, result, index, dim, 0.0, 0.0);
		}
	}

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
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		size_t seq = index.seq();

		double alpha_value = source[seq];

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double helper = (1.0/pow(2.0, static_cast<int>(l+1))) * (static_cast<double>(i));

		// integration
		result[seq] = (  ( (fr-fl) * (helper) ) - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l)))) * alpha_value)) );    // diagonal entry

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

	/**
	 * recursive function for the calculation of Down wit Bounding Box support
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 */
	void recBB(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		size_t seq = index.seq();

		double alpha_value = source[seq];

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double helper = (1.0/pow(2.0, static_cast<int>(l+1))) * (q * static_cast<double>(i));

		// integration
		result[seq] = (  ( (fr-fl) * (helper + (0.5*t)) )
							  - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l))) * q) * alpha_value)) );    // diagonal entry

		// dehierarchisation
		double fm = ((fl+fr)/2.0) + alpha_value;

		if(!index.hint())
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				recBB(source, result, index, dim, fl, fm);
			}

			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				recBB(source, result, index, dim, fm, fr);
			}

			index.up(dim);
		}
	}
};

} // namespace detail

} // namespace sg

#endif /* XPHIDPHIDOWNBBLINEAR_HPP */
