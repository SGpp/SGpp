// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETPARALLELMPI_HPP
#define POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETPARALLELMPI_HPP

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystemDirichlet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

/**
 * This class uses sgpp::pde::OperationEllipticPDESolverSystemDirichlet
 * to define a solver system for the Poission Equation.
 *
 * For the mult-routine only the Laplace-Operator is required
 *
 * There is a parallelization over all operators, required
 * to solve the poisson equation.
 */
class PoissonEquationEllipticPDESolverSystemDirichletParallelMPI
    : public sgpp::pde::OperationEllipticPDESolverSystemDirichlet {
 protected:
  std::unique_ptr<sgpp::base::OperationMatrix> Laplace_Inner;
  std::unique_ptr<sgpp::base::OperationMatrix> Laplace_Complete;

  void applyLOperatorComplete(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  void applyLOperatorInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 public:
  /**
   * Constructor
   *
   * @param SparseGrid reference to a sparse grid on which the Poisson Equation should be solved
   * @param rhs the right hand side for solving the elliptic PDE
   */
  PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(sgpp::base::Grid& SparseGrid,
                                                             sgpp::base::DataVector& rhs);

  /**
   * Destructor
   */
  virtual ~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI();

  void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  sgpp::base::DataVector* generateRHS();
};
}  // namespace parallel
}  // namespace sgpp

#endif /* POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETPARALLELMPI_HPP */
