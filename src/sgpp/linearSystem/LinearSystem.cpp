/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)       */
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

#include "linearSystem/LinearSystem.hpp"
#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/OperationB.hpp"
#include "operation/OperationMatrix.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

	LinearSystem::LinearSystem(DataVector& d, DataVector& y, OperationB& bOp, OperationMatrix& cOp, Grid& grid, double l)
	{
		this->data = &d;
		this->y = &y;
		this->bOperation = &bOp;
		this->cOperation = &cOp;
		this->grid = &grid;
		this->l = l;
	}
	
	void LinearSystem::apply(DataVector& vector, DataVector& result)
	{
		try{
			size_t size = (*this->data).getSize();
			DataVector* temp = new  DataVector(size);
	        int M = this->data->getSize();
	    
	        this->bOperation->multTranspose(vector, *this->data, *temp);
	        this->bOperation->mult(*temp, *this->data, result);
	
	        
	        temp = new DataVector(vector.getSize());
	        this->cOperation->mult(vector, *temp);
	        result.axpy(M * this->l, *temp);
		}
		catch(int ex){
			std::cout<<"Exception #"<<ex<<"catched"<<std::endl;
		}
        
        /*
         * @todo: impelment C Operators for (khakhutv)
         * - ratio
         * - levelsum
         * - energy
         * - copy
         * - pseudounit 
         */
	}
	
	void LinearSystem::getRightHandSide(DataVector& result)
	{
        this->bOperation->mult(*(this->y), *(this->data), result);
	}
	
	int LinearSystem::getNumGridPoints()
	{
		return this->grid->getStorage()->size();
	}
	
	int LinearSystem::getNumInputPoints()
	{
		return this->data->getSize();
	}
	

}