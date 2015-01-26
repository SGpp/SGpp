/* *****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/misc/pde/basis/linear/boundary/operation/OperationLaplaceEnhancedLinearBoundary.hpp>

#include <sgpp/misc/pde/basis/linear/boundary/algorithm_sweep/LaplaceEnhancedDownBBLinearBoundary.hpp>
#include <sgpp/misc/pde/basis/linear/boundary/algorithm_sweep/LaplaceEnhancedUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationLaplaceEnhancedLinearBoundary::OperationLaplaceEnhancedLinearBoundary(SGPP::base::GridStorage* storage) : UpDownOneOpDimEnhanced(storage) {
    }

    OperationLaplaceEnhancedLinearBoundary::OperationLaplaceEnhancedLinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef) : UpDownOneOpDimEnhanced(storage, coef) {
    }

    OperationLaplaceEnhancedLinearBoundary::~OperationLaplaceEnhancedLinearBoundary() {
    }

    void OperationLaplaceEnhancedLinearBoundary::up(SGPP::base::DataMatrix& alpha, SGPP::base::DataMatrix& result, size_t dim) {
      LaplaceEnhancedUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<LaplaceEnhancedUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceEnhancedLinearBoundary::down(SGPP::base::DataMatrix& alpha, SGPP::base::DataMatrix& result, size_t dim) {
      LaplaceEnhancedDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<LaplaceEnhancedDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
