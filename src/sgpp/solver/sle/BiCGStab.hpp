/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#ifndef BICGSTAB_HPP
#define BICGSTAB_HPP

#include "solver/SLESolver.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

#include <iostream>

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
