/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONHIERARCHISATIONMODLINEAR_HPP
#define OPERATIONHIERARCHISATIONMODLINEAR_HPP

#include "operation/common/OperationHierarchisation.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * Hierarchisation on sparse grid, mod linear case
 *
 * @version $HEAD$
 */
class OperationHierarchisationModLinear : public OperationHierarchisation
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationHierarchisationModLinear(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationHierarchisationModLinear() {}

	virtual void doHierarchisation(DataVector& node_values);
	virtual void doDehierarchisation(DataVector& alpha);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONHIERARCHISATIONMODLINEAR_HPP */
