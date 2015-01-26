/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include <sgpp/misc/pde/basis/linear/noboundary/operation/OperationLaplaceEnhancedLinear.hpp>

#include <sgpp/misc/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedDownBBLinear.hpp>
#include <sgpp/misc/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace pde {

OperationLaplaceEnhancedLinear::OperationLaplaceEnhancedLinear(SGPP::base::GridStorage* storage) :
		UpDownOneOpDimEnhanced(storage) {
}

OperationLaplaceEnhancedLinear::OperationLaplaceEnhancedLinear(SGPP::base::GridStorage* storage,
		SGPP::base::DataVector& coef) :
		UpDownOneOpDimEnhanced(storage, coef) {
}

OperationLaplaceEnhancedLinear::~OperationLaplaceEnhancedLinear() {
}

void OperationLaplaceEnhancedLinear::up(SGPP::base::DataMatrix& alpha, SGPP::base::DataMatrix& result, size_t dim) {
	LaplaceEnhancedUpBBLinear func(this->storage);
	SGPP::base::sweep<LaplaceEnhancedUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceEnhancedLinear::down(SGPP::base::DataMatrix& alpha, SGPP::base::DataMatrix& result, size_t dim) {
	LaplaceEnhancedDownBBLinear func(this->storage);
	SGPP::base::sweep<LaplaceEnhancedDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
