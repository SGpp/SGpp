/******************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/operation/OperationEllipticPDESolverSystemFreeBoundaries.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationEllipticPDESolverSystemFreeBoundaries::OperationEllipticPDESolverSystemFreeBoundaries(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs) {
    }

    OperationEllipticPDESolverSystemFreeBoundaries::~OperationEllipticPDESolverSystemFreeBoundaries() {
    }

    void OperationEllipticPDESolverSystemFreeBoundaries::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      applyLOperator(alpha, result);
    }

    SGPP::base::DataVector* OperationEllipticPDESolverSystemFreeBoundaries::generateRHS() {
      return this->rhs;
    }

  }
}

