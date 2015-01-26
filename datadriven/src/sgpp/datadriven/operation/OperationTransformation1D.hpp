/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
#ifndef OPERATIONTRANSFORMATION1D_HPP
#define OPERATIONTRANSFORMATION1D_HPP

#include "base/grid/Grid.hpp"
#include <cstring>

namespace sg {
namespace datadriven {

/**
 * Sample 1D Probability Density Function
 */

class OperationTransformation1D {
public:
	OperationTransformation1D() {
	}
	virtual ~OperationTransformation1D() {
	}

	/**
	 * Transform 1d
	 * @param alpha1d
	 * @param coord1d
	 * @return
	 */
	virtual double doTransformation1D(base::DataVector* alpha1d, double coord1d) = 0;
};

}
}
#endif /* OPERATIONTRANSFORMATION1D_HPP */
