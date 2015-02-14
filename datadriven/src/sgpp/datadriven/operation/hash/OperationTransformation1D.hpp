// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONTRANSFORMATION1D_HPP
#define OPERATIONTRANSFORMATION1D_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
	virtual float_t doTransformation1D(base::DataVector* alpha1d, float_t coord1d) = 0;
};

}
}
#endif /* OPERATIONTRANSFORMATION1D_HPP */