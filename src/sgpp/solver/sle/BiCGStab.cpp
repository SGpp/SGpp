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

#include "solver/sle/BiCGStab.hpp"
#include <cmath>

namespace sg
{

BiCGStab::BiCGStab(size_t imax, double epsilon) : SLESolver(imax, epsilon)
{
}

BiCGStab::~BiCGStab()
{
}

void BiCGStab::solve(OperationMatrix& SystemMatrix, DataVector& alpha, DataVector& b, bool reuse, bool verbose, double max_threshold)
{
	this->nIterations = 1;
	double epsilonSqd = this->myEpsilon* this->myEpsilon;

	if (reuse == false)
	{
		// Choose x0
		alpha.setAll(0.0);
	}

	//Calculate r0
	DataVector r(alpha.getSize());
	SystemMatrix.mult(alpha, r);
	r.sub(b);

	double delta_0 = r.dotProduct(r)*epsilonSqd;
	double delta = 0.0;

	if (verbose == true)
	{
		std::cout <<  "delta_0 " << delta_0 << std::endl;
	}

	//Choose r0 as r
	DataVector rZero(r);
	// Set p as r0
	DataVector p(rZero);

	double rho = rZero.dotProduct(r);
	double rho_new = 0.0;
	double sigma = 0.0;
	double a = 0.0;
	double omega = 0.0;
	double beta = 0.0;

	DataVector s(alpha.getSize());
	DataVector v(alpha.getSize());
	DataVector w(alpha.getSize());

	s.setAll(0.0);
	v.setAll(0.0);
	w.setAll(0.0);

	while (this->nIterations < this->nMaxIterations)
	{
		// s  = Ap
		s.setAll(0.0);
		SystemMatrix.mult(p, s);

		//std::cout << "s " << s.get(0) << " " << s.get(1)  << std::endl;

		sigma = s.dotProduct(rZero);
		if (fabs(sigma) == 0.0)
		{
			break;
		}

		a = rho/sigma;

		// w = r - a*s
		w = r;
		w.axpy((-1.0)*a, s);

		// v = Aw
		v.setAll(0.0);
		SystemMatrix.mult(w, v);

		//std::cout << "v " << v.get(0) << " " << v.get(1)  << std::endl;

		omega = v.dotProduct(w) / v.dotProduct(v);

		// x = x - a*p - omega*w
		alpha.axpy((-1.0)*a, p);
		alpha.axpy((-1.0)*omega, w);

		// r = r - a*s - omega*v
		r.axpy((-1.0)*a, s);
		r.axpy((-1.0)*omega, v);

		rho_new = r.dotProduct(rZero);

		delta = r.dotProduct(r);

		if (verbose == true)
		{
			std::cout << "delta: " << delta << std::endl;
		}
		this->residuum = delta;

		// Stop in case of better accurency
		if (delta < delta_0)
		{
			break;
		}

		beta = rho_new/rho * a/omega;
		rho = rho_new;

		// p = r + beta*(p - omega*s)
		p.axpy((-1.0)*omega, s);
		p.mult(beta);
		p.add(r);

		this->nIterations++;
	}
}

}
