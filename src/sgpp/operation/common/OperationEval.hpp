/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONEVAL_HPP
#define OPERATIONEVAL_HPP

#include "data/DataVector.hpp"

#ifdef WINDOWS
#pragma warning(disable: 4267)
#endif

namespace sg
{

/**
 * Operation that evaluates the current sparse grid function defined
 * by the coefficient vector @em alpha at a given point.
 *
 * @todo (pflueged) Use eval(DataVector& alpha, DataVector& point) as default
 */
class OperationEval
{
public:
	/**
	 * Default constructor.
	 */
	OperationEval() {}

	/**
	 * Destructor
	 */
	virtual ~OperationEval() {}

	/**
	 * Evaluates the sparse grid function at a given point.
	 *
	 * @param alpha The coefficients of the sparse grid's basis functions
	 * @param point The coordinates of the evaluation point
	 */
	virtual double eval(DataVector& alpha, std::vector<double>& point) = 0;

	/**
	 * Evaluates the sparse grid function at a given point.
	 *
	 * @param alpha The coefficients of the sparse grid's basis functions
	 * @param point The coordinates of the evaluation point
	 */
	virtual double eval(DataVector& alpha, DataVector& point)
	{
		std::vector<double> p;
		for(size_t i = 0; i < point.getSize(); i++)
		{
			p.push_back(point[i]);
		}
		return eval(alpha, p);
	}
};

}

#endif /* OPERATIONEVAL_HPP */
