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
using sgpp::optimization::optimizer::LogBarrier;
using sgpp::optimization::optimizer::SquaredPenalty;
using sgpp::optimization::optimizer::AugmentedLagrangian;

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
  AlambC_ptr->generateb(b);

  // define function to optimize
  // ||(A+lambda*C)*alpha - b||^2
  std::function<double(const DataVector&)> residual_norm =
      [storage_size, AlambC_ptr, b_ptr](const DataVector& alpha) -> double {
    DataVector resultVec(storage_size);
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    // std::cout << "numSamples: " << numSamples << std::endl;
    AlambC_ptr->mult(*alpha_ptr, resultVec);
    resultVec.sub(*b_ptr);
    // std::cout << "return res:" << resultVec.dotProduct(resultVec) << std::endl;
    return resultVec.dotProduct(resultVec);
  };

  // define gradient for residual_norm
  // = 2*((A+lambda*C)*alpha - b).T * (A+lambda*C)
  // A and C symmetric: = 2*(A+lambda*C) * ((A+lambda*C)*alpha - b)
  std::function<double(const DataVector&, DataVector&)> residual_norm_grad =
      [storage_size, AlambC_ptr, b_ptr](const DataVector& alpha, DataVector& gradient) -> double {
    std::cout << "res norm grad" << std::endl;
    std::cout << "alpha:" << alpha.toString() << std::endl;
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    DataVector rightResult(storage_size);
    AlambC_ptr->mult(*alpha_ptr, rightResult);
    rightResult.sub(*b_ptr);  // = ((A+lambda*C)*alpha - b)
    AlambC_ptr->mult(rightResult, gradient);
    gradient.mult(2.0);                          // gradient f
    std::cout << "return grad:" << rightResult.dotProduct(rightResult) << std::endl;
    std::cout << "grad:" << gradient.toString() << std::endl;
    return rightResult.dotProduct(rightResult);  // return function value f
  };

  // define inequality constraint
  // for alpha: eval at all grid points >= 0
  // i.e. PSI*alpha >= 0 but the optimizers take constraints of the form g(x) <= 0
  // so we return -(PSI*alpha)
  std::function<void(const DataVector&, DataVector&)> in_equ =
      [op_eval, storage_size, numSamples, gridStorage, dims](const DataVector& alpha,
                                                             DataVector& result) -> void {
    std::cout << "in equ" << std::endl;
    std::cout << "alpha:" << alpha.toString() << std::endl;
    GridPoint gp;
    for (size_t i = 0; i < storage_size; i++) {
      gp = gridStorage->getPoint(i);
      DataVector coords(dims);
      gp.getStandardCoordinates(coords);
      result[i] = op_eval->eval(alpha, coords);
    }
    std::cout << "in equ:" << result.toString() << std::endl;
  };

  // define inequality constraint gradient
  // = -PSI

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
      jacobi.set(i, j, val);
    }
  }

  std::function<void(const DataVector&, DataVector&, DataMatrix&)> in_equ_grad =
    [op_eval, storage_size, gridStorage, dims, jacobi_ptr](const DataVector& alpha, DataVector& result,
                                                           DataMatrix& gradient) -> void {
    std::cout << "in equ grad" << std::endl;
    std::cout << "alpha:" << alpha.toString() << std::endl;
    GridPoint gp;
    DataVector coords(dims);
    for (size_t i = 0; i < storage_size; i++) {
      std::cout << "i:" << i << std::endl;
      gp = gridStorage->getPoint(i);
      gp.getStandardCoordinates(coords);
      std::cout << "coords:" << coords.toString() << std::endl;
      result[i] = op_eval->eval(alpha, coords);
      std::cout << "result[i]:" << result[i] << std::endl;
    }
    gradient = *jacobi_ptr;
    std::cout << "return in equ gradient:" << gradient.toString() << std::endl;
    std::cout << "return in equ :" << result.toString() << std::endl;
  };

  std::function<void(const DataVector&, DataVector&)> equ =
    [](const DataVector& alpha, DataVector& result) -> void {
    result.setAll(0.0);
  };

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

  // SquaredPenalty optimizer(wrapped_res_norm, wrapped_res_norm_grad,
                           // wrapped_in_equ, wrapped_in_equ_grad,
                           // wrapped_equ, wrapped_equ_grad);

  AugmentedLagrangian optimizer(wrapped_res_norm, wrapped_res_norm_grad,
                           wrapped_in_equ, wrapped_in_equ_grad,
                           wrapped_equ, wrapped_equ_grad);
  optimizer.optimize();
  const DataVector best_alpha = optimizer.getOptimalPoint();
  std::cout << best_alpha.toString() << std::endl;
  checkPositive(*grid, best_alpha);
}

int main(int argc, char** argv) {
  size_t d = 2;
  int level = 3;
  GridType gridType = sgpp::base::GridType::Linear;
  double lambda = 0.0;

  sgpp::base::RegularGridConfiguration gridConfig;
  // sgpp::base::AdpativityConfiguration adaptivityConfig;
  // sgpp::solver::SLESolverConfiguration solverConfig;
  // sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  // sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
  // sgpp::datadriven::SGDEConfiguration sgdeConfig;
  gridConfig.dim_ = d;
  gridConfig.level_ = level;
  gridConfig.type_ = gridType;
  gridConfig.maxDegree_ = 5;
  // adaptivityConfig.numRefinements_ = 0;
  // adaptivityConfig.noPoints_ = 15;
  // solverConfig.threshold_ = 1e-15;
  // solverConfig.verbose_ = false;
  // regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;
  // crossvalidationConfig.enable_ = false;
  // sgdeConfig.makePositive_ = false;

  // sgpp::datadriven::LearnerSGDE learnerSGDE(gridConfig, adaptivityConfig, solverConfig,
  //                                           regularizationConfig, crossvalidationConfig,
  //                                           sgdeConfig);
  DataMatrix trainSamples(1000, 2);
  createSamples(trainSamples, 1234567);
  // learnerSGDE.initialize(trainSamples);
  // std::shared_ptr<Grid> grid(&learnerSGDE.getGrid());
  // DataVector alpha(learnerSGDE.getSurpluses());
  std::cout << trainSamples.getNrows() << std::endl;
  solve(trainSamples, gridConfig, lambda);
  // std::cout << trainSamples.toString() << std::endl;
  return 0;
}
