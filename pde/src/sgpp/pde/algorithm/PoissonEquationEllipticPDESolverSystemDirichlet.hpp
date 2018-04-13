// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP
#define POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystemDirichlet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * This class uses OperationEllipticPDESolverSystemDirichlet
 * to define a solver system for the Poission Equation.
 *
 * For the mult-routine only the Laplace-Operator is required
 */
class PoissonEquationEllipticPDESolverSystemDirichlet
    : public OperationEllipticPDESolverSystemDirichlet {
 protected:
  sgpp::base::OperationMatrix* Laplace_Inner;
  sgpp::base::OperationMatrix* Laplace_Complete;

  void applyLOperatorComplete(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  void applyLOperatorInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 public:
  /**
   * Constructor
   *
   * @param SparseGrid reference to a sparse grid on which the Poisson Equation should be solved
   * @param rhs the right hand side for solving the elliptic PDE
   */
  PoissonEquationEllipticPDESolverSystemDirichlet(sgpp::base::Grid& SparseGrid,
                                                  sgpp::base::DataVector& rhs);

  /**
   * Destructor
   */
  virtual ~PoissonEquationEllipticPDESolverSystemDirichlet();
};
}  // namespace pde
}  // namespace sgpp

#endif /* POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP */
