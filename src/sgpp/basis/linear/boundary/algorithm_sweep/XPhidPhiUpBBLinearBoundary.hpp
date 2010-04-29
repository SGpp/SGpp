/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef XPHIDPHIUPBBLINEARBOUNDARY_HPP
#define XPHIDPHIUPBBLINEARBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp"

namespace sg
{

namespace detail
{

/**
 * up-operation in dimension dim. for use with sweep
 */
class XPhidPhiUpBBLinearBoundary : public XPhidPhiUpBBLinear
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	XPhidPhiUpBBLinearBoundary(GridStorage* storage) : XPhidPhiUpBBLinear(storage)
	{
	}

	/**
	 * Destructor
	 */
	~XPhidPhiUpBBLinearBoundary()
	{
	}

	/**
	 * This operations performs the calculation of up in the direction of dimension <i>dim</i>
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		this->q = this->boundingBox->getIntervalWidth(dim);
		this->t = this->boundingBox->getIntervalOffset(dim);

		bool useBB = false;

		if (this->q != 1.0 || this->t != 0.0)
		{
			useBB = true;
		}

		// get boundary values
		double fl = 0.0;
		double fr = 0.0;

		if(useBB)
		{
			if(!index.hint())
			{
				index.top(dim);

				if(!this->storage->end(index.seq()))
				{
					recBB(source, result, index, dim, fl, fr);
				}

				index.left_levelzero(dim);
			}

			size_t seq_left;
			size_t seq_right;

			// left boundary
			seq_left = index.seq();

			// right boundary
			index.right_levelzero(dim);
			seq_right = index.seq();

			// check boundary conditions
			if (this->boundingBox->hasDirichletBoundaryLeft(dim))
			{
				result[seq_left] = 0.0; // source[seq_left];
			}
			else
			{
				// up
				//////////////////////////////////////
				result[seq_left] = fl;

				result[seq_left] += (source[seq_right] * (((1.0/6.0)*this->q) + (0.5*this->t)));
			}

			if (this->boundingBox->hasDirichletBoundaryRight(dim))
			{
				result[seq_right] = 0.0; //source[seq_right];
			}
			else
			{
				result[seq_right] = fr;
			}

			index.left_levelzero(dim);
		}
		else
		{
			if(!index.hint())
			{
				index.top(dim);

				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, fl, fr);
				}

				index.left_levelzero(dim);
			}

			size_t seq_left;
			size_t seq_right;

			// left boundary
			seq_left = index.seq();

			// right boundary
			index.right_levelzero(dim);
			seq_right = index.seq();

			// check boundary conditions
			if (this->boundingBox->hasDirichletBoundaryLeft(dim))
			{
				result[seq_left] = 0.0; // source[seq_left];
			}
			else
			{
				// up
				//////////////////////////////////////
				result[seq_left] = fl;

				result[seq_left] += (source[seq_right] * (1.0/6.0));
			}

			if (this->boundingBox->hasDirichletBoundaryRight(dim))
			{
				result[seq_right] = 0.0; //source[seq_right];
			}
			else
			{
				result[seq_right] = fr;
			}

			index.left_levelzero(dim);
		}
	}
};

} // namespace detail

} // namespace sg

#endif /* XPHIDPHIUPBBLINEARBOUNDARY_HPP */
