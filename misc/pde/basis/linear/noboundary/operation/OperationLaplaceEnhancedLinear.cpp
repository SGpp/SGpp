/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include "misc/pde/basis/linear/noboundary/operation/OperationLaplaceEnhancedLinear.hpp"

#include "misc/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedDownBBLinear.hpp"
#include "misc/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
namespace pde {

OperationLaplaceEnhancedLinear::OperationLaplaceEnhancedLinear(sg::base::GridStorage* storage) :
		UpDownOneOpDimEnhanced(storage) {
}

OperationLaplaceEnhancedLinear::OperationLaplaceEnhancedLinear(sg::base::GridStorage* storage,
		sg::base::DataVector& coef) :
		UpDownOneOpDimEnhanced(storage, coef) {
}

OperationLaplaceEnhancedLinear::~OperationLaplaceEnhancedLinear() {
}

void OperationLaplaceEnhancedLinear::up(sg::base::DataMatrix& alpha, sg::base::DataMatrix& result, size_t dim) {
	LaplaceEnhancedUpBBLinear func(this->storage);
	sg::base::sweep<LaplaceEnhancedUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceEnhancedLinear::down(sg::base::DataMatrix& alpha, sg::base::DataMatrix& result, size_t dim) {
	LaplaceEnhancedDownBBLinear func(this->storage);
	sg::base::sweep<LaplaceEnhancedDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}
