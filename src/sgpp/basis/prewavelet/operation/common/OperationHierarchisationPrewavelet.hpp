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

#ifndef OPERATIONHIERARCHISATIONPREWAVELET_HPP
#define OPERATIONHIERARCHISATIONPREWAVELET_HPP

#include "operation/common/OperationHierarchisation.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{
namespace base
{

/**
 * Hierarchisation on sparse grid with prewavelets and no boundary. Please note, that there is
 * no efficient way to directly calculate the hierarchical surpluses for the prewavelet base.
 * But there is a fast and efficient way to transform hierarchical surpluses of the normal
 * linear ansatzfunctions into the prewavelet base and vice versa. Thus, we use the normal
 * hierarchisation of the linear basis and afterwards transform the resulting Vector into
 * prewavelets (see ConvertLinearToPrewavelet.hpp). For the Dehierarchisation, this process is
 * reversed (see ConvertPrewaveletToLinear.hpp)
 */
class OperationHierarchisationPrewavelet: public OperationHierarchisation
{
public:
	/**
	 * Constructor of OperationHierarchisationPrewavelet
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationHierarchisationPrewavelet(GridStorage* storage,GridStorage* shadowStorage) :
		storage(storage), shadowStorage(shadowStorage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationHierarchisationPrewavelet()
	{
	}

	virtual void doHierarchisation(DataVector& node_values);
	virtual void doDehierarchisation(DataVector& alpha);

protected:
	/// Pointer to the grid's GridStorage object
	GridStorage* storage;
	GridStorage* shadowStorage;

	void expandGrid();

	void shrinkGrid();
};

}
}

#endif /* OPERATIONHIERARCHISATIONPREWAVELET_HPP */
