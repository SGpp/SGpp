// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP
#define OPERATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystem.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/base/grid/common/DirichletGridConverter.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace pde {

/**
 * Defines a System that is used to solve elliptic partial
 * differential equations. So an instance of this class has to pass to
 * any SLE Solver used in SGpp.
 *
 * \f$L \vec{u} = rhs\f$
 *
 * L: space discretization (L-Operator)
 * rhs: right hand side
 *
 * This class is a specialized version of OperationEllipticPDESolverSystem which
 * exploits Dirichlet boundary conditions. Since there are no degrees of freedom
 * on on the boundaries the iterative solver (CG or BiCGSTAB) has only to take
 * inner grid points into account.
 *
 * The inner grid is constructed during the constructor call!
 */
class OperationEllipticPDESolverSystemDirichlet : public OperationEllipticPDESolverSystem {
 protected:
  /// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
  sgpp::base::DataVector* alpha_inner;
  /// Routine to modify the boundaries/inner points of the grid
  sgpp::base::DirichletUpdateVector* BoundaryUpdate;
  /// Class that allows a simple conversion between a grid with and a without boundary points
  sgpp::base::DirichletGridConverter* GridConverter;
  /// Pointer to the inner grid object
  sgpp::base::Grid* InnerGrid;
  /// rhs for the inner grid
  sgpp::base::DataVector* rhs_inner;

  /**
   * applies the PDE's system matrix, on complete grid - with boundaries
   *
   * @param alpha the coefficients of the sparse grid's ansatzfunctions
   * @param result reference to the sgpp::base::DataVector into which the result is written
   */
  virtual void applyLOperatorComplete(sgpp::base::DataVector& alpha,
                                      sgpp::base::DataVector& result) = 0;

  /**
   * applies the PDE's system matrix, on inner grid only
   *
   * @param alpha the coefficients of the sparse grid's ansatzfunctions
   * @param result reference to the sgpp::base::DataVector into which the result is written
   */
  virtual void applyLOperatorInner(sgpp::base::DataVector& alpha,
                                   sgpp::base::DataVector& result) = 0;

 public:
  /**
   * Constructor
   *
   * @param SparseGrid the grid, for which the system should be solved
   * @param rhs the right hand side of the corresponding system
   */
  OperationEllipticPDESolverSystemDirichlet(sgpp::base::Grid& SparseGrid,
                                            sgpp::base::DataVector& rhs);

  /**
   * Destructor
   */
  virtual ~OperationEllipticPDESolverSystemDirichlet();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual sgpp::base::DataVector* generateRHS();

  /**
   * gets a pointer to the sparse grids coefficients used in the CG method to solve
   * one timestep. This is useful because (direchlet) boundaries can be skipped when
   * solving the system.
   *
   * @return alpha vector for CG method
   */
  virtual sgpp::base::DataVector* getGridCoefficientsForCG();

  /**
   * Gets the solution for the complete grid
   *
   * @param Solution sgpp::base::DataVector that must have a dimension equal to the bound's grid
   * dimension, the result is written to Solution
   * @param SolutionInner Solution on the inner grid
   */
  virtual void getSolutionBoundGrid(sgpp::base::DataVector& Solution,
                                    sgpp::base::DataVector& SolutionInner);
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONELLIPTICPDESOLVERMATRIXDIRICHLET_HPP */
