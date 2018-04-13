// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP
#define OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP

#include <sgpp/solver/operation/hash/OperationParabolicPDESolverSystem.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/base/grid/common/DirichletGridConverter.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Defines a System that is used to solve parabolic partial
 * differential equations. So an instance of this class has to pass to
 * any ODE Solver used in SGpp.
 *
 * \f$A \dot{u} = L \vec{u}\f$
 *
 * A: mass matrix
 * L: space discretization (L-Operator)
 *
 * This class defines an elliptic problem in every timestep which is solved
 * using an iterative SLE solver, that solving step is integrated in the
 * ODE Solver.
 *
 * This class is a specialized version of OperationParabolicPDESolverSystem which
 * exploit Dirichlet boundary conditions. Hence there are no degrees of freedom
 * on on the boundaries the iterative solver (CG or BiCGSTAB) has only to take
 * inner grid points into account.
 */
class OperationParabolicPDESolverSystemDirichlet
    : public sgpp::solver::OperationParabolicPDESolverSystem {
 protected:
  /// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
  sgpp::base::DataVector* alpha_inner;
  /// Routine to modify the boundaries/inner points of the grid
  sgpp::base::DirichletUpdateVector* BoundaryUpdate;
  /// Class that allows a simple conversion between a grid with and a without boundary points
  sgpp::base::DirichletGridConverter* GridConverter;
  /// Pointer to the inner grid object
  sgpp::base::Grid* InnerGrid;

  /**
   * applies the PDE's mass matrix, on complete grid - with boundaries
   *
   * @param alpha the coefficients of the sparse grid's ansatzfunctions
   * @param result reference to the sgpp::base::DataVector into which the result is written
   */
  virtual void applyMassMatrixComplete(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result) = 0;

  /**
   * applies the PDE's system matrix, on complete grid - with boundaries
   *
   * @param alpha the coefficients of the sparse grid's ansatzfunctions
   * @param result reference to the sgpp::base::DataVector into which the result is written
   */
  virtual void applyLOperatorComplete(sgpp::base::DataVector& alpha,
                                      sgpp::base::DataVector& result) = 0;

  /**
   * applies the PDE's mass matrix, on inner grid only
   *
   * @param alpha the coefficients of the sparse grid's ansatzfunctions
   * @param result reference to the sgpp::base::DataVector into which the result is written
   */
  virtual void applyMassMatrixInner(sgpp::base::DataVector& alpha,
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
   */
  OperationParabolicPDESolverSystemDirichlet();

  /**
   * Destructor
   */
  virtual ~OperationParabolicPDESolverSystemDirichlet();

  /**
   * Multiplicates a vector with the matrix
   *
   * @param alpha sgpp::base::DataVector that contains the ansatzfunctions' coefficients
   * @param result sgpp::base::DataVector into which the result of the space discretization
   * operation is stored
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  /**
   * generates the right hand side of the system
   *
   * @return returns the rhs
   */
  virtual sgpp::base::DataVector* generateRHS();

  virtual sgpp::base::DataVector* getGridCoefficientsForCG();
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP */
