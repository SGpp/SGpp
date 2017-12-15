// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp_optimization.hpp>

#include <functional>
#include <random>
#include <vector>
#ifdef USE_CGAL
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_solution.h>
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
typedef CGAL::Quadratic_program_solution<ET> Solution;
typedef CGAL::Quadratic_program_from_iterators<
 double**,                                           //for A
 double*,                                                 // for b
 CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
 bool*,                                                // for fl
 double*,                                                 // for l
 bool*,                                                // for fu
 double*,                                                 // for u
 double**,                                                // for D
 double*>                                                 // for c
Program;
#endif

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::GridType;
using sgpp::base::OperationEval;
using sgpp::base::OperationMatrix;
using sgpp::base::OperationMultipleEval;
using sgpp::datadriven::DensitySystemMatrix;
using sgpp::optimization::WrapperScalarFunction;
using sgpp::optimization::WrapperScalarFunctionGradient;
using sgpp::optimization::WrapperVectorFunction;
using sgpp::optimization::WrapperVectorFunctionGradient;
using sgpp::optimization::optimizer::AugmentedLagrangian;
using sgpp::optimization::optimizer::LogBarrier;
using sgpp::optimization::optimizer::SquaredPenalty;

void randu(DataVector& rvar, std::mt19937& generator) {
  std::normal_distribution<double> distribution(0.5, 0.1);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
    // make sure the sample is in the unit cube
    while (rvar[j] < 0 || rvar[j] > 1) rvar[j] = distribution(generator);
  }
}

void createSamples(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();
  std::mt19937 generator(seedValue);
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

void checkPositive(Grid& grid, const DataVector& alpha) {
  GridStorage* gridStorage = &grid.getStorage();
  size_t storage_size = gridStorage->getSize();
  size_t dims = grid.getDimension();
  GridPoint gp;
  DataVector coords(dims);
  OperationEval* op_eval = sgpp::op_factory::createOperationEval(grid);
  for (size_t i = 0; i < storage_size; i++) {
    gp = gridStorage->getPoint(i);
    gp.getStandardCoordinates(coords);
    // if (op_eval->eval(alpha, coords) < 0) {
    std::cout << "f(" << coords.toString() << ")=" << op_eval->eval(alpha, coords) << std::endl;
    // }
  }
}

void solve(DataMatrix& samples, sgpp::base::RegularGridConfiguration& gridConfig, double lambda) {
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createGrid(gridConfig));
  GridStorage* gridStorage = &(grid->getStorage());
  grid->getGenerator().regular(gridConfig.level_);
  size_t storage_size = gridStorage->getSize();
  std::cout << "storage size:" << storage_size << std::endl;
  size_t numSamples = samples.getNrows();
  size_t dims = samples.getNcols();
  OperationMatrix* A_op = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  OperationMultipleEval* B_op = sgpp::op_factory::createOperationMultipleEval(*grid, samples);
  OperationMatrix* C_op = sgpp::op_factory::createOperationLaplace(*grid);
  DataVector b(storage_size);
  std::cout << "b size:" << b.getSize() << std::endl;
  DataVector* b_ptr = &b;
  DensitySystemMatrix AlambC(A_op, B_op, C_op, lambda, numSamples);
  DensitySystemMatrix* AlambC_ptr = &AlambC;
  OperationEval* op_eval = sgpp::op_factory::createOperationEval(*grid);
  DataVector* shift_vec = new DataVector(storage_size, 20.0);
  AlambC_ptr->generateb(b);

  // define function to optimize
  // ||(A+lambda*C)*alpha - b||^2
  std::function<double(const DataVector&)> residual_norm =
      [storage_size, AlambC_ptr, b_ptr, shift_vec](const DataVector& alpha) -> double {
    DataVector resultVec(storage_size);
    DataVector alpha_cpy(alpha);
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    alpha_cpy.mult(40.0);
    alpha_cpy.sub(*shift_vec);
    // std::cout << "res alpha:" << alpha_cpy.toString() << std::endl;
    // std::cout << "numSamples: " << numSamples << std::endl;
    AlambC_ptr->mult(alpha_cpy, resultVec);
    resultVec.sub(*b_ptr);
    // std::cout << "return res:" << resultVec.dotProduct(resultVec) << std::endl;
    return resultVec.dotProduct(resultVec);
  };

  // define gradient for residual_norm
  // = 2*((A+lambda*C)*alpha - b).T * (A+lambda*C)
  // A and C symmetric: = (2*(A+lambda*C) * ((A+lambda*C)*alpha - b)).T
  std::function<double(const DataVector&, DataVector&)> residual_norm_grad =
      [storage_size, AlambC_ptr, b_ptr, shift_vec](const DataVector& alpha,
                                                   DataVector& gradient) -> double {
    // std::cout << "res norm grad" << std::endl;
    // std::cout << "res grad alpha:" << alpha.toString() << std::endl;
    DataVector alpha_cpy(alpha);
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    DataVector rightResult(storage_size);
    alpha_cpy.mult(40.0);
    alpha_cpy.sub(*shift_vec);
    // std::cout << "res grad alpha:" << alpha.toString() << std::endl;
    AlambC_ptr->mult(alpha_cpy, rightResult);
    rightResult.sub(*b_ptr);  // = ((A+lambda*C)*alpha - b)
    AlambC_ptr->mult(rightResult, gradient);
    gradient.mult(2.0);  // gradient f
    // std::cout << "return grad:" << rightResult.dotProduct(rightResult) << std::endl;
    // std::cout << "grad:" << gradient.toString() << std::endl;
    return rightResult.dotProduct(rightResult);  // return function value f
  };
  // define inequality constraint
  // for alpha: eval at all grid points >= 0
  // i.e. PSI*alpha >= 0 but the optimizers take constraints of the form g(x) <= 0
  // so we return -(PSI*alpha)
  std::function<void(const DataVector&, DataVector&)> in_equ =
      [op_eval, storage_size, numSamples, gridStorage, dims, shift_vec](
          const DataVector& alpha, DataVector& result) -> void {
    // std::cout << "in equ" << std::endl;
    // std::cout << "alpha:" << alpha.toString() << std::endl;
    DataVector alpha_cpy(alpha);
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    alpha_cpy.mult(40.0);
    alpha_cpy.sub(*shift_vec);
    GridPoint gp;
    for (size_t i = 0; i < storage_size; i++) {
      gp = gridStorage->getPoint(i);
      DataVector coords(dims);
      gp.getStandardCoordinates(coords);
      result[i] = -op_eval->eval(alpha_cpy, coords);
    }
    // std::cout << "return in equ:" << result.toString() << std::endl;
  };


  // precalc jacobi matrix
  DataMatrix jacobi(storage_size, storage_size);
  DataMatrix* jacobi_ptr = &jacobi;
  DataVector tmp_alpha(storage_size);
  GridPoint gp;
  DataVector coords(dims);
  for (size_t i = 0; i < storage_size; i++) {
    gp = gridStorage->getPoint(i);
    gp.getStandardCoordinates(coords);
    for (size_t j = 0; j < storage_size; j++) {
      tmp_alpha.setAll(0.0);
      tmp_alpha.set(j, 1.0);
      double val = op_eval->eval(tmp_alpha, coords);
      jacobi.set(i, j, -40. * val);
    }
  }

  // define inequality constraint gradient
  // = -PSI
  std::function<void(const DataVector&, DataVector&, DataMatrix&)> in_equ_grad =
      [op_eval, storage_size, gridStorage, dims, jacobi_ptr, shift_vec](
          const DataVector& alpha, DataVector& result, DataMatrix& gradient) -> void {
    // std::cout << "in equ grad" << std::endl;
    // std::cout << "alpha:" << alpha.toString() << std::endl;
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    DataVector alpha_cpy(alpha);
    alpha_cpy.mult(40.0);
    alpha_cpy.sub(*shift_vec);
    GridPoint gp;
    DataVector coords(dims);
    for (size_t i = 0; i < storage_size; i++) {
      // std::cout << "i:" << i << std::endl;
      gp = gridStorage->getPoint(i);
      gp.getStandardCoordinates(coords);
      // std::cout << "coords:" << coords.toString() << std::endl;
      result[i] = -op_eval->eval(alpha_cpy, coords);
      // std::cout << "result[i]:" << result[i] << std::endl;
    }
    gradient = *jacobi_ptr;
    // std::cout << "return in equ gradient:" << gradient.toString() << std::endl;
    // std::cout << "return in equ :" << result.toString() << std::endl;
  };

  std::function<void(const DataVector&, DataVector&)> equ =
      [](const DataVector& alpha, DataVector& result) -> void { result.setAll(0.0); };

  std::function<void(const DataVector&, DataVector&, DataMatrix&)> equ_grad =
      [](const DataVector& alpha, DataVector& result, DataMatrix& jacobi) -> void {
    result.setAll(0.0);
    jacobi.setAll(0.0);
  };
  // wrapping functions
  WrapperScalarFunction wrapped_res_norm(storage_size, residual_norm);
  WrapperScalarFunctionGradient wrapped_res_norm_grad(storage_size, residual_norm_grad);
  WrapperVectorFunction wrapped_in_equ(storage_size, storage_size, in_equ);
  WrapperVectorFunctionGradient wrapped_in_equ_grad(storage_size, storage_size, in_equ_grad);
  WrapperVectorFunction wrapped_equ(storage_size, storage_size, equ);
  WrapperVectorFunctionGradient wrapped_equ_grad(storage_size, storage_size, equ_grad);

  // optimization
  // LogBarrier optimizer(wrapped_res_norm, wrapped_res_norm_grad,
  // wrapped_in_equ, wrapped_in_equ_grad);

  // SquaredPenalty optimizer(wrapped_res_norm, wrapped_res_norm_grad, wrapped_in_equ,
                           // wrapped_in_equ_grad, wrapped_equ, wrapped_equ_grad,
                           // 20000, 1e-15, 1e-15);

  AugmentedLagrangian optimizer(wrapped_res_norm, wrapped_res_norm_grad,
                                wrapped_in_equ, wrapped_in_equ_grad,
                                wrapped_equ, wrapped_equ_grad,
                                20000, 1e-15, 1e-15);
  double al[17] = {18.90769281,  -9.09911171,  -9.0176581 ,   0.10785147,
                   -1.18044529,  -0.28708884,   0.15533845,  -9.38380473,
                   -9.00915268,  -0.02802101,  -2.24950203,  -1.58041974,
                   -0.09580067,   5.12791209,   4.54942521,   4.66378412,   4.32550775};
  DataVector cmp_alpha(al, 17);
  cmp_alpha.add(*shift_vec);
  cmp_alpha.mult(1./40.);
  optimizer.setStartingPoint(cmp_alpha);
  optimizer.optimize();
  DataVector best_alpha = optimizer.getOptimalPoint();
  // best_alpha.mult(40.);
  // best_alpha.sub(*shift_vec);

  std::cout << "Found optimal alpha:" << best_alpha.toString() << std::endl;
  std::cout << "Iterations:" << optimizer.getHistoryOfOptimalValues().getSize() << std::endl;
  // checkPositive(*grid, best_alpha);
  // checkPositive(*grid, cmp_alpha);
  std::cout << "res norm best:" << residual_norm(best_alpha) << std::endl;
  std::cout << "res norm cmp:" << residual_norm(cmp_alpha) << std::endl;
}

void solve_cgal(DataMatrix& samples, sgpp::base::RegularGridConfiguration& gridConfig, double lambda) {
  // CGAL documentation: https://doc.cgal.org/latest/QP_solver/index.html
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createGrid(gridConfig));
  GridStorage* gridStorage = &(grid->getStorage());
  grid->getGenerator().regular(gridConfig.level_);
  size_t storage_size = gridStorage->getSize();
  std::cout << "storage size:" << storage_size << std::endl;
  size_t numSamples = samples.getNrows();
  size_t dims = samples.getNcols();
  DataMatrix M(storage_size, storage_size);
  DataMatrix C(storage_size, storage_size);
  OperationMatrix* A_op = sgpp::op_factory::createOperationLTwoDotExplicit(&M, *grid);
  OperationMultipleEval* B_op = sgpp::op_factory::createOperationMultipleEval(*grid, samples);
  OperationMatrix* C_op = sgpp::op_factory::createOperationLaplaceExplicit(&C, *grid);
  DataVector b(storage_size);
  DataVector q(storage_size);
  DensitySystemMatrix AlambC(A_op, B_op, C_op, lambda, numSamples);
  OperationEval* op_eval = sgpp::op_factory::createOperationEval(*grid);
  AlambC.generateb(b);

  // setting up System matrix: M + lambda*C
  C.mult(lambda);
  M.add(C);

  // Skalaraproduct term of quadratic program: q.T x
  M.mult(b, q);

  // Quadratic program matrix P = M*M.T (M is symmetric)
  double** P_it = new double*[storage_size];
  for (size_t i = 0; i < M.getNcols(); i++) {
    DataVector col(storage_size);
    DataVector tmp(storage_size);
    M.getColumn(i, col);
    M.mult(col, tmp);
    P_it[i] = new double[storage_size];
    for (size_t j = i; j < storage_size; j++) {
      P_it[i][j] = tmp.get(j);
      P_it[j][i] = tmp.get(j);
      std::cout << P_it[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "---------------" << std::endl;

  // getting all grid points for interpolation matrix
  DataMatrix gridPoints(storage_size, dims);
  GridPoint gp;
  for (size_t i = 0; i < storage_size; i++) {
    gp = gridStorage->getPoint(i);
    DataVector coords(dims);
    gp.getStandardCoordinates(coords);
    gridPoints.setRow(i, coords);
  }

  // interpolation matrix (grid point values when multiplied with alpha-vec)
  OperationMultipleEval* G_op = sgpp::op_factory::createOperationMultipleEval(*grid, gridPoints);
  double** G_it = new double*[storage_size];
  for (size_t i = 0; i < storage_size; i++) {
    G_it[i] = new double[storage_size];
  }
  // G_it[i] should contain the i-th COLUMN of G
  for (size_t i = 0; i < storage_size; i++) {
    DataVector result(storage_size);
    DataVector alpha(storage_size, 0.0);
    alpha.set(i, 1.0);
    G_op->mult(alpha, result);
    for (size_t j = 0; j < storage_size; j++) {
      G_it[i][j] = result.get(j);
    }
  }

  for (size_t i = 0; i < storage_size; i++) {
    for (size_t j = 0; j < storage_size; j++) {
      std::cout << G_it[j][i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << q.toString() << std::endl;
  // constraint relation (i.e. greater than zero)
  CGAL::Const_oneset_iterator<CGAL::Comparison_result> r(CGAL::SMALLER);

  // lower upper bounds unused:
  bool* bounded = new bool[storage_size];
  for (size_t i = 0; i < storage_size; i++) {
    bounded[i] = false;
  }
  DataVector bounds(storage_size, 0.0);
  // define the quadratic Programm
  Program qp(storage_size, storage_size,  // size of problem
             G_it, b.getPointer(), r,     // constraints
             bounded, bounds.getPointer(), bounded, bounds.getPointer(),  // bounds
             P_it, q.getPointer()  // optimization goal
             );
  Solution s = CGAL::solve_quadratic_program(qp, ET());
  Solution::Variable_value_iterator it = s.variable_values_begin();
  DataVector best_alpha(storage_size);
  for (size_t i = 0; i < storage_size; i++) {
    best_alpha.set(i, to_double(*it));
    it++;
  }
  std::cout << "---------------------------" << std::endl;
  std::cout << "best alpha:" << best_alpha.toString() << std::endl;
  std::cout << "objective function:" << to_double(s.objective_value()) << std::endl;
}

int main(int argc, char** argv) {
  size_t d = 2;
  int level = 3;
  GridType gridType = sgpp::base::GridType::Linear;
  double lambda = 0.0;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = d;
  gridConfig.level_ = level;
  gridConfig.type_ = gridType;
  gridConfig.maxDegree_ = 5;
  DataMatrix trainSamples(1000, 2);
  createSamples(trainSamples, 1234567);
  std::cout << trainSamples.getNrows() << std::endl;
  // solve(trainSamples, gridConfig, lambda);
#ifdef USE_CGAL
  solve_cgal(trainSamples, gridConfig, lambda);
#endif
  return 0;
}
