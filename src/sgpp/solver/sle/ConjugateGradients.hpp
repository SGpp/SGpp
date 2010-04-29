/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CONJUGATEGRADIENTS_HPP
#define CONJUGATEGRADIENTS_HPP

#include "solver/SLESolver.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

#include <iostream>

namespace sg
{

class ConjugateGradients : public SLESolver
{
private:


public:
	/**
	 * Std-Constructor
	 */
	ConjugateGradients(size_t imax, double epsilon);

	/**
	 * Std-Destructor
	 */
	virtual ~ConjugateGradients();

	virtual void solve(OperationMatrix& SystemMatrix, DataVector& alpha, DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0);

	// Define functions for observer pattern in python

	/**
	 * function that signals the start of the CG method (used in python)
	 */
	virtual void starting();

	/**
	 * function that signals the start of the calculation of the CG method (used in python)
	 */
	virtual void calcStarting();

	/**
	 * function that signals that one iteration step of the CG method has been completed (used in python)
	 */
	virtual void iterationComplete();

	/**
	 * function that signals the finish of the cg method (used in python)
	 */
	virtual void complete();
};

}

#endif /* CONJUGATEGRADIENTS_HPP */
