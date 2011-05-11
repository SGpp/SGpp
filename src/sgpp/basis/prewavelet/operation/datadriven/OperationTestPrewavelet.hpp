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

#ifndef OPERATIONTESTPREWAVELET_HPP
#define OPERATIONTESTPREWAVELET_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"
using namespace sg::base;

namespace sg
{
namespace datadriven
{

/**
 * This class implements OperationTest for a grids with prewavelet basis ansatzfunctions without boundaries
 */
class OperationTestPrewavelet : public OperationTest
{
public:
	/**
	 * Construtor of OperationTestPrewavelet
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationTestPrewavelet(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestPrewavelet() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
};

}
}

#endif /* OPERATIONTESTPREWAVELET_HPP */
