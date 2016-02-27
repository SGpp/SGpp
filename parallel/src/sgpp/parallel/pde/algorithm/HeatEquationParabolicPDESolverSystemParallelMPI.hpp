// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELMPI_HPP
#define HEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELMPI_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace parallel {

/**
 * This class implements the ParabolicPDESolverSystem for the
 * Heat Equation parallelized with MPI.
 */
class HeatEquationParabolicPDESolverSystemParallelMPI
    : public sgpp::pde::OperationParabolicPDESolverSystemDirichlet {
 private:
  /// the heat coefficient
  double a;
  /// the Laplace Operation (Stiffness Matrix), on boundary grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpLaplaceBound;
  /// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpMassBound;
  /// the Laplace Operation (Stiffness Matrix), on inner grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpLaplaceInner;
  /// the LTwoDotProduct Operation (Mass Matrix), on inner grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpMassInner;

  void applyMassMatrixComplete(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  void applyLOperatorComplete(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  void applyMassMatrixInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  void applyLOperatorInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param alpha the sparse grid's coefficients
   * @param a the heat coefficient
   * @param TimestepSize the size of one timestep used in the ODE Solver
   * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for
   * explicit Euler,
   *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
   */
  HeatEquationParabolicPDESolverSystemParallelMPI(sgpp::base::Grid& SparseGrid,
                                                  sgpp::base::DataVector& alpha, double a,
                                                  double TimestepSize,
                                                  std::string OperationMode = "ExEul");

  /**
   * Std-Destructor
   */
  virtual ~HeatEquationParabolicPDESolverSystemParallelMPI();

  void finishTimestep();

  void coarsenAndRefine(bool isLastTimestep = false);

  void startTimestep();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual sgpp::base::DataVector* generateRHS();
};
}  // namespace parallel
}  // namespace sgpp

#endif /* HEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELMPI_HPP */
