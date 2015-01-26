/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
#ifndef OPERATIONROSENBLATTTRANSFORMATION_HPP
#define OPERATIONROSENBLATTTRANSFORMATION_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

/**
 * Sampling on all dimensions
 */

class OperationRosenblattTransformation {
public:
	OperationRosenblattTransformation() {
	}
	virtual ~OperationRosenblattTransformation() {
	}

	/**
	 * Rosenblatt Transformation with mixed starting dimension
	 *
	 * @param alpha Coefficient vector for current grid
	 * @param points Input DataMatrix (rows: # of samples, columns: # of dims)
	 * @param pointscdf Output DataMatrix (rows: # of samples, columns: # of dims)
	 */
	virtual void doTransformation(base::DataVector* alpha,
			base::DataMatrix* points, base::DataMatrix* pointscdf) = 0;

	/**
	 * Rosenblatt Transformation with fixed starting dimension
	 *
	 * @param alpha Coefficient vector for current grid
	 * @param points Input DataMatrix (rows: # of samples, columns: # of dims)
	 * @param pointscdf Output DataMatrix (rows: # of samples, columns: # of dims)
	 * @param dim_start starting dimension
	 */
	virtual void doTransformation(base::DataVector* alpha,
			base::DataMatrix* points, base::DataMatrix* pointscdf,
			size_t dim_start) = 0;
};

}
}
#endif /* OPERATIONROSENBLATTTRANSFORMATION_HPP */
