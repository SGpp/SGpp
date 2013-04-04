/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "base/exception/algorithm_exception.hpp"
#include "pde/operation/PdeOpFactory.hpp"

using namespace sg::op_factory;

namespace sg {
  namespace pde {

    PoissonEquationEllipticPDESolverSystemDirichlet::PoissonEquationEllipticPDESolverSystemDirichlet(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
      this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
      this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
    }

    PoissonEquationEllipticPDESolverSystemDirichlet::~PoissonEquationEllipticPDESolverSystemDirichlet() {
      delete this->Laplace_Complete;
      delete this->Laplace_Inner;
    }

    void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      Laplace_Inner->mult(alpha, result);
    }

    void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      Laplace_Complete->mult(alpha, result);
    }

  }
}
