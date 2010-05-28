/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONHIERARCHISATION_HPP
#define OPERATIONHIERARCHISATION_HPP

#include "data/DataVector.hpp"

namespace sg
{

/**
 * This class implements the hierarchisation and dehierarchisation on the sparse grid
 */
class OperationHierarchisation
{
public:
	/**
	 * Constructor
	 */
	OperationHierarchisation() {}

	/**
	 * Destructor
	 */
	virtual ~OperationHierarchisation() {}

	/**
	 * Implements the hierarchisation on a sprase grid
	 *
	 * @param node_values the functions values in the node base
	 */
	virtual void doHierarchisation(DataVector& node_values) = 0;

	/**
	 * Implements the dehierarchisation on a sprase grid
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 */
	virtual void doDehierarchisation(DataVector& alpha) = 0;
};

}

#endif /* OPERATIONHIERARCHISATION_HPP */
