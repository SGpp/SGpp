// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETVECTORIZEDMPI_HPP
#define POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETVECTORIZEDMPI_HPP

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystemDirichlet.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

/**
 * This class uses OperationEllipticPDESolverSystemDirichlet
 * to define a solver system for the Poission Equation.
 *
 * For the mult-routine only the Laplace-Operator is required
 */
class PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : public
  SGPP::pde::OperationEllipticPDESolverSystemDirichlet {
 protected:
  std::unique_ptr<SGPP::base::OperationMatrix> Laplace_Inner;
  std::unique_ptr<SGPP::base::OperationMatrix> Laplace_Complete;

  void applyLOperatorComplete(SGPP::base::DataVector& alpha,
                              SGPP::base::DataVector& result);

  void applyLOperatorInner(SGPP::base::DataVector& alpha,
                           SGPP::base::DataVector& result);

 public:
  /**
   * Constructor
   *
   * @param SparseGrid reference to a sparse grid on which the Poisson Equation should be solved
   * @param rhs the right hand side for solving the elliptic PDE
   */
  PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(
    SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs);

  /**
   * Destructor
   */
  virtual ~PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI();
};

}
}

#endif /* POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETVECTORIZEDMPI_HPP */
