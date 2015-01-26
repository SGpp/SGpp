/* ****************************************************************************
 * Copyright (C) 2012 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author: Emily Mo-Hellenbrand
#ifndef OPERATIONROSENBLATTTRANSFORMATION1DLINEAR_HPP
#define OPERATIONROSENBLATTTRANSFORMATION1DLINEAR_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/OperationTransformation1D.hpp>

namespace sg {
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
