/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef OPERATIONCONVERTPREWAVELET_HPP
#define OPERATIONCONVERTPREWAVELET_HPP

#include "base/operation/OperationConvert.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg
{
namespace base
{

/**
 *
 */
class OperationConvertPrewavelet : public OperationConvert
{
public:
	/**
	 * Constructor of OperationHierarchisationPrewavelet
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationConvertPrewavelet(GridStorage* storage, GridStorage* shadowstorage) :
		storage(storage),shadowstorage(shadowstorage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationConvertPrewavelet()
	{
	}

	virtual void doConvertToLinear(DataVector& alpha);
	virtual void doConvertFromLinear(DataVector& alpha);

protected:
	/// Pointer to the grid's GridStorage object
	GridStorage* storage;
	GridStorage* shadowstorage;
};

}
}

#endif /* OPERATIONCONVERTPREWAVELET_HPP */
