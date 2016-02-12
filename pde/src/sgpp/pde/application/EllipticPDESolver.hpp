// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ELLIPTICPDESOLVER_HPP
#define ELLIPTICPDESOLVER_HPP

#include <sgpp/pde/application/PDESolver.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

/**
 * This class extends the PDESolver with functions that are needed to
 * solve elliptic PDEs
 *
 */
class EllipticPDESolver : public PDESolver {
 public:
  /**
   * Std-Constructor of the solver
   */
  EllipticPDESolver();

  /**
   * Std-Destructor of the solver
   */
  virtual ~EllipticPDESolver();

  /**
   * abstract method to solve an elliptic PDE. All solver of elliptic PDEs
   * have to implement this method.
   *
   * @param alpha the coefficients of the Sparse Gird's basis functions will be in this vector after
   * solving
   * @param rhs the right hand side of the SLE
   * @param maxCGIterations the maximum of interation in the CG solver
   * @param epsilonCG the epsilon used in the CG
   * @param verbose enables verbose output during solving
   */
  virtual void solvePDE(SGPP::base::DataVector& alpha, SGPP::base::DataVector& rhs,
                        size_t maxCGIterations, float_t epsilonCG, bool verbose = false) = 0;
};
}  // namespace pde
}  // namespace SGPP

#endif /* ELLIPTICPDESOLVER_HPP */
