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

#include "solver/euler/ExplicitEuler.hpp"
#include "solver/cg/ConjugateGradients.hpp"

namespace sg
{

ExplicitEuler::ExplicitEuler(size_t imax, double timestepSize) : ODESolver(imax, timestepSize)
{
	this->residuum = 0.0;
}

ExplicitEuler::~ExplicitEuler()
{
}

void ExplicitEuler::solve(OperationMatrix& SystemMatrix, DataVector& alpha, bool verbose)
{
	// expliziter Euler
	////////////////////
	DataVector temp(alpha.getSize());

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		temp.setAll(0.0);

		SystemMatrix.mult(alpha, temp);

		alpha.axpy(this->myEpsilon, temp);

		//std::cout << temp.toString() << std::endl;
		//std::cout << alpha.toString() << std::endl;
		//std::cout << this->myEpsilon << std::endl;
	}

	// Heun ??
	////////////////////
	/*DataVector temp(alpha.getSize());
	DataVector alphaEuler(alpha.getSize());
	DataVector alphaHeun(alpha.getSize());

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		double r = alpha.get(1);
		double l = alpha.get(0);
		temp.setAll(0.0);
		alphaEuler = alpha;
		alphaHeun = alpha;

		// Euler step:
		SystemMatrix.mult(alphaEuler, temp);
		alphaEuler.axpy(this->myEpsilon, temp);
		alphaEuler.mult(static_cast<double>(i+1)*(this->myEpsilon));

		alphaHeun.mult(static_cast<double>(i)*(this->myEpsilon));
		alphaHeun.add(alphaEuler);
		alphaHeun.mult(this->myEpsilon*0.5);

		alpha.add(alphaHeun);

		alpha.set(0,l);
		alpha.set(1,r);

		//std::cout << temp.toString() << std::endl;
		//std::cout << alpha.toString() << std::endl;
		//std::cout << this->myEpsilon << std::endl;
	}*/

	// Crank Nicholson???
	//////////////////////
	/*DataVector temp(alpha.getSize());
	DataVector cgresult(alpha.getSize());
	DataVector rhs(alpha.getSize());
    ConjugateGradients myCG(this->nMaxIterations, this->myEpsilon);

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		temp.setAll(0.0);
		rhs = alpha;

		SystemMatrix.mult(alpha, temp);
		//std::cout << temp.toString() << std::endl;
		//std::cout << this->myEpsilon << std::endl;
		rhs.axpy(this->myEpsilon, temp);

		alpha.setAll(0.0);
		//std::cout << alpha.toString() << std::endl;
		//std::cout << rhs.toString() << std::endl;

	    // slove the system of linear equations, SystemMatrix is wrong Matrix ??!!!!!
	    myCG.solve(SystemMatrix, alpha, rhs, false, true, -1.0);
	}*/
}

}
