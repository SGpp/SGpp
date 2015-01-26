/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

using namespace SGPP::op_factory;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    PoissonEquationEllipticPDESolverSystemDirichlet::PoissonEquationEllipticPDESolverSystemDirichlet(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
      this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
      this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
    }

    PoissonEquationEllipticPDESolverSystemDirichlet::~PoissonEquationEllipticPDESolverSystemDirichlet() {
      delete this->Laplace_Complete;
      delete this->Laplace_Inner;
    }

    void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      Laplace_Inner->mult(alpha, result);
    }

    void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      Laplace_Complete->mult(alpha, result);
    }

  }
}
