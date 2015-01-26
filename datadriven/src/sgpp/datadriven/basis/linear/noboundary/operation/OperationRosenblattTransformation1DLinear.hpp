// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONROSENBLATTTRANSFORMATION1DLINEAR_HPP
#define OPERATIONROSENBLATTTRANSFORMATION1DLINEAR_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/OperationTransformation1D.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {
class OperationRosenblattTransformation1DLinear: public OperationTransformation1D {
protected:
	base::Grid* grid;
public:
	OperationRosenblattTransformation1DLinear(base::Grid* grid);
	virtual ~OperationRosenblattTransformation1DLinear();

	/**
	 * Rosenblatt Transformation 1D
	 * @param alpha1d
	 * @param coord1d
	 * @return
	 */
	double doTransformation1D(base::DataVector* alpha1d, double coord1d);
};

}
}

#endif /* OPERATIONROSENBLATTTRANSFORMATION1DLINEAR_HPP */