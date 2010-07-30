/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONEVAL_HPP
#define OPERATIONEVAL_HPP

#include "data/DataVector.hpp"

#ifdef WINDOWS
#pragma warning(disable: 4267)
#endif

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
};

}

#endif /* OPERATIONEVAL_HPP */
