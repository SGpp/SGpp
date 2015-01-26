/******************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/operation/OperationEllipticPDESolverSystemFreeBoundaries.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

namespace sg {
  namespace pde {

    OperationEllipticPDESolverSystemFreeBoundaries::OperationEllipticPDESolverSystemFreeBoundaries(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs) {
    }

    OperationEllipticPDESolverSystemFreeBoundaries::~OperationEllipticPDESolverSystemFreeBoundaries() {
    }

    void OperationEllipticPDESolverSystemFreeBoundaries::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      applyLOperator(alpha, result);
    }

    sg::base::DataVector* OperationEllipticPDESolverSystemFreeBoundaries::generateRHS() {
      return this->rhs;
    }

  }
}

