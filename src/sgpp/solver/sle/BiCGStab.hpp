/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BICGSTAB_HPP
#define BICGSTAB_HPP

#include "solver/SLESolver.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

#include <iostream>
using namespace sg::base;

namespace sg
{

class BiCGStab : public SLESolver
{
private:


public:
	/**
	 * Std-Constructor
	 */
	BiCGStab(size_t imax, double epsilon);

	/**
	 * Std-Destructor
	 */
	virtual ~BiCGStab();

	/**
	 * max_threashold is ignored in this solver
	 *
	 * Reference:
	 * http://www.iue.tuwien.ac.at/phd/heinreichsberger/node70.html
	 * http://www.numerik.math.tu-graz.ac.at/kurse/lgs/SIMNET6.pdf
	 * http://netlib.org
	 */
	virtual void solve(OperationMatrix& SystemMatrix, DataVector& alpha, DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0);
};

}

#endif /*BICGSTAB_HPP */
