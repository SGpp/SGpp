/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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

#include "sgpp.hpp"

#include "basis/basis.hpp"

#include "basis/linearboundary/operation/datadriven/OperationTestLinearBoundary.hpp"

#include "data/DataVector.hpp"

namespace sg
{

double OperationTestLinearBoundary::test(DataVector& alpha, DataVector& data, DataVector& classes)
{
	linearboundaryBase<unsigned int, unsigned int> base;
	return test_dataset_boundary(this->storage, base, alpha, data, classes);
}

double OperationTestLinearBoundary::testWithCharacteristicNumber(DataVector& alpha, DataVector& data, DataVector& classes, DataVector& charaNumbers)
{
	linearboundaryBase<unsigned int, unsigned int> base;
	return test_dataset_boundaryWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}

