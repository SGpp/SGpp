/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)       		 */
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

#ifndef LINEARSYSTEM_HPP
#define LINEARSYSTEM_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/OperationB.hpp"
#include "operation/OperationMatrix.hpp"

namespace sg
{

class LinearSystem
{
	public:
		/**
		 * Constructor
		 * 
		 * @param data DataVector points
		 * @param result DataVector values
		 * @param bOberation OperationB B operator
		 * @param cOperation OperationMatrix C operator
		 * @param grid Grid grid
		 * @param l double regression parameter
		 */
		LinearSystem(DataVector& data, DataVector& result, OperationB& bOperation, OperationMatrix& cOperation, Grid& grid, double l);
	
	
		/**
		 * Apply linear system on the vector
		 * 
		 * @param vector DataVector [in] vector the linear system should be applied to
		 * @param result DataVector [out] vector for results
		 */
		void apply(DataVector& vector,  DataVector& result);
	
		/**
		 * Calculate the right side of the linear system
		 * 
		 * @param result DataVector [out] vector for result
		 */
		void getRightHandSide(DataVector& result);
		
		/**
		 * Returns the number of points on the grid
		 * 
		 * @return int number of points on the grid
		 */
		int getNumGridPoints();
		
		/**
		 * Return number of point in the dataset
		 * 
		 * @return int number of point in the dataset
		 */
		int getNumInputPoints();
		
	private:
		DataVector* data;
		DataVector* y;
		OperationB* bOperation;
		OperationMatrix* cOperation;
		double l;
		Grid* grid;
//		int numInputPoints;
		
};

}

#endif /* LINEARSYSTEM_HPP */