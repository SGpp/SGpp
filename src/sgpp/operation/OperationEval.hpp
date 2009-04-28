/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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

#ifndef OPERATIONEVAL_HPP
#define OPERATIONEVAL_HPP

#include "data/DataVector.h"

namespace sg
{

/**
 * Operation the evaluate the function that is applied the current Sparse Grid at a given point
 */
class OperationEval
{
public:
	/**
	 * Constructor
	 */
	OperationEval() {}

	/**
	 * Destructor
	 */
	virtual ~OperationEval() {}

	/**
	 * Evaluates the grid's function at a given point
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 * @param point the coordinates of the evaluation point
	 */
	virtual double eval(DataVector& alpha, std::vector<double>& point) = 0;

	/**
	 * Evaluates the grid's function at a given point
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 * @param point the coordinates of the evaluation point
	 */
	virtual double eval(DataVector& alpha, DataVector& point)
	{
		std::vector<double> p;
		for(size_t i = 0; i < point.getDim(); i++)
		{
			p.push_back(point[i]);
		}
		return eval(alpha, p);
	}

	/**
	 * Evaluates the grid's function at a given point and tests the class
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 * @param data the coordinates of the evaluation points
	 * @param classes DataVector the holds the class information
	 */
	virtual double test(DataVector& alpha, DataVector& data, DataVector& classes) = 0;
};

}

#endif /* OPERATIONEVAL_HPP */
