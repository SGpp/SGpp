/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)


#include "operation/common/OperationMatrix.hpp"

#include "data/DataVector.hpp"

namespace sg
{

/**
 * Implementation of identity Operation for all kinds of grids
 */
class OperationIdentity: public OperationMatrix
{
public:
	/**
	 * Construtor of OperationIdentity
	 */
	OperationIdentity()
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationIdentity() {}

	void mult(DataVector& alpha, DataVector& result)
	{
		result = alpha;
	}
};

}
