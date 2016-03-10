// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONELLIPTICPDESOLVERSYSTEMFREEBOUNDARIES_HPP
#define OPERATIONELLITPICPDESOLVERSYSTEMFREEBOUNDARIES_HPP

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystem.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Defines a System that is used to solve elliptic partial
 * differential equations. So an instance of this class has to pass to
 * any SLE Solver used in SGpp, here degrees of freedom exists on
 * the boundaries!
 *
 * \f$L \vec{u} = rhs\f$
 *
 * L: space discretization (L-Operator)
 * rhs: right hand sider)
 *
 */
class OperationEllipticPDESolverSystemFreeBoundaries : public OperationEllipticPDESolverSystem {
 protected:
  /**
   * applies the PDE's system matrix, on complete grid - with boundaries
   *
   * @param alpha the coefficients of the sparse grid's ansatzfunctions
   * @param result reference to the sgpp::base::DataVector into which the result is written
   */
  virtual void applyLOperator(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) = 0;

 public:
  /**
   * Constructor
   *
   * @param SparseGrid the grid, for which the system should be solved
   * @param rhs the right hand side of the corresponding system
   */
  OperationEllipticPDESolverSystemFreeBoundaries(sgpp::base::Grid& SparseGrid,
                                                 sgpp::base::DataVector& rhs);

  /**
   * Destructor
   */
  virtual ~OperationEllipticPDESolverSystemFreeBoundaries();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual sgpp::base::DataVector* generateRHS();
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONELLITPTICPDESOLVERMATRIXFREEBOUNDARIES_HPP */
