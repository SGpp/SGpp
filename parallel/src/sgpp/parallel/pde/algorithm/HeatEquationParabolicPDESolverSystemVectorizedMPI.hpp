// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HEATEQUATIONPARABOLICPDESOLVERSYSTEMVECTORIZEDMPI_HPP
#define HEATEQUATIONPARABOLICPDESOLVERSYSTEMVECTORIZEDMPI_HPP

#include <sgpp/parallel/pde/operation/OperationParabolicPDESolverSystemDirichletCombined.hpp>
#include <sgpp/parallel/pde/operation/OperationParabolicPDEMatrixCombined.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace parallel {

/**
 * This class implements the ParabolicPDESolverSystem for the
 * Heat Equation using vectorized and iterative operators
 */
class HeatEquationParabolicPDESolverSystemVectorizedMPI
    : public sgpp::parallel::OperationParabolicPDESolverSystemDirichletCombined {
 private:
  /// the heat coefficient
  double a;
  /// the Laplace Operation, on boundary grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpLaplaceBound;
  /// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpLTwoBound;
  /// the Laplace Operation, on Inner grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpLaplaceInner;
  /// the LTwoDotProduct Operation (Mass Matrix), on Inner grid
  std::unique_ptr<sgpp::base::OperationMatrix> OpLTwoInner;
  /// the combination of the LTwoDotProduct Operation (Mass Matrix) and the Laplace Operation, on
  /// Inner grid
  std::unique_ptr<sgpp::parallel::OperationParabolicPDEMatrixCombined> OpLTwoDotLaplaceInner;
  /// the combination of the LTwoDotProduct Operation (Mass Matrix) and the Laplace Operation, on
  /// Bound grid
  std::unique_ptr<sgpp::parallel::OperationParabolicPDEMatrixCombined> OpLTwoDotLaplaceBound;

  virtual void applyMassMatrixComplete(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result);

  virtual void applyLOperatorComplete(sgpp::base::DataVector& alpha,
                                      sgpp::base::DataVector& result);

  virtual void applyMassMatrixInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual void applyLOperatorInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual void applyMassMatrixLOperatorInner(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result);

  virtual void applyMassMatrixLOperatorBound(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result);

  virtual void setTimestepCoefficientInner(double timestep_coefficient);

  virtual void setTimestepCoefficientBound(double timestep_coefficient);

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param alpha the sparse grid's coefficients
   * @param a the heat coefficient
   * @param TimestepSize the size of one timestep used in the ODE Solver
   * @param OperationMode specifies in which solver this matrix is used, valid values are:
   *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
   */
  HeatEquationParabolicPDESolverSystemVectorizedMPI(sgpp::base::Grid& SparseGrid,
                                                    sgpp::base::DataVector& alpha, double a,
                                                    double TimestepSize,
                                                    std::string OperationMode = "ImEul");

  /**
   * Std-Destructor
   */
  virtual ~HeatEquationParabolicPDESolverSystemVectorizedMPI();

  virtual void finishTimestep();

  virtual void coarsenAndRefine(bool isLastTimestep = false);

  void startTimestep();
};
}  // namespace parallel
}  // namespace sgpp

#endif /* HEATEQUATIONPARABOLICPDESOLVERSYSTEMVECTORIZEDMPI_HPP */
