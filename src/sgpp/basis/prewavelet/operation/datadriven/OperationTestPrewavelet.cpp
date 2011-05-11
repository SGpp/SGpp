/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#include "basis/prewavelet/operation/datadriven/OperationTestPrewavelet.hpp"

#include "data/DataVector.hpp"

#include "data/DataMatrix.hpp"
using namespace sg::base;

namespace sg
{
namespace datadriven
{

double OperationTestPrewavelet::test(DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	prewavelet_base<unsigned int, unsigned int> base;
	return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestPrewavelet::testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	prewavelet_base<unsigned int, unsigned int> base;
	return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestPrewavelet::testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	prewavelet_base<unsigned int, unsigned int> base;
	return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes, charaNumbers);
}

}
}

