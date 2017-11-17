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

void solve(DataMatrix& samples, sgpp::base::RegularGridConfiguration gridConfig, double lambda) {
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createGrid(gridConfig));
  GridStorage* gridStorage = &grid->getStorage();
  size_t storage_size = gridStorage->getSize();
  size_t numSamples = samples.getNrows();
  size_t dims = samples.getNcols();
  OperationMatrix* A_op = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  OperationMultipleEval* B_op = sgpp::op_factory::createOperationMultipleEval(*grid, samples);
  OperationMatrix* C_op = sgpp::op_factory::createOperationLaplace(*grid);
  DataVector b(numSamples);
  DataVector* b_ptr = &b;
  DensitySystemMatrix AlambC(A_op, B_op, C_op, lambda, numSamples);
  DensitySystemMatrix* AlambC_ptr = &AlambC;
  OperationEval* op_eval = sgpp::op_factory::createOperationEval(*grid);
  AlambC_ptr->generateb(b);
  // define function to optimize
  // ||(A+lambda*C)*alpha - b||^2
  std::function<double(const DataVector&)> residual_norm =
      [numSamples, AlambC_ptr, b_ptr](const DataVector& alpha) -> double {
    DataVector alpha_copy(alpha);
    DataVector resultVec(numSamples);
    AlambC_ptr->mult(alpha_copy, resultVec);
    resultVec.sub(*b_ptr);
    return resultVec.dotProduct(resultVec);
  };

  // define gradient for residual_norm
  // = 2*((A+lambda*C)*alpha - b).T * (A+lambda*C)
  // A and C symmetric: = 2*(A+lambda*C) * ((A+lambda*C)*alpha - b)
  std::function<double(const DataVector&, DataVector&)> residual_norm_grad =
      [numSamples, AlambC_ptr, b_ptr](const DataVector& alpha, DataVector& gradient) -> double {
    DataVector alpha_copy(alpha);
    DataVector rightResult(numSamples);
    AlambC_ptr->mult(alpha_copy, rightResult);
    rightResult.sub(*b_ptr);  // = ((A+lambda*C)*alpha - b)
    AlambC_ptr->mult(rightResult, gradient);
    gradient.mult(2.0);                          // gradient f
    return rightResult.dotProduct(rightResult);  // return function value f
  };

  // define inequality constraint
  // for alpha: eval at all grid points >= 0
  // i.e. PSI*alpha >= 0
  std::function<void(const DataVector&, DataVector&)> in_equ =
      [op_eval, storage_size, numSamples, gridStorage, dims](const DataVector& alpha,
                                                             DataVector& result) -> void {
    GridPoint gp;
    for (size_t i = 0; i < storage_size; i++) {
      gp = gridStorage->getPoint(i);
      DataVector coords(dims);
      gp.getStandardCoordinates(coords);
      result[i] = op_eval->eval(alpha, coords);
    }
  };

  // define inequality constraint gradient
  // = PSI
  std::function<void(const DataVector&, DataVector&, DataMatrix&)> in_equ_grad =
      [op_eval, storage_size, gridStorage, dims](const DataVector& alpha, DataVector& result,
                                                 DataMatrix& jacobi) -> void {
    DataVector* alpha_ptr = &const_cast<DataVector&>(alpha);
    GridPoint gp;
    for (size_t i = 0; i < storage_size; i++) {
      gp = gridStorage->getPoint(i);
      DataVector coords(dims);
      DataVector tmp_alpha(storage_size);
      gp.getStandardCoordinates(coords);
      result[i] = op_eval->eval(*alpha_ptr, coords);
      // dumm gelöst
      for (size_t j = 0; j < storage_size; i++) {
        tmp_alpha.setAll(0.0);
        tmp_alpha.set(j, 1.0);
        double val = op_eval->eval(*alpha_ptr, coords);
        jacobi.set(i, j, val);
      }
    }
  };
  // wrapping functions
  WrapperScalarFunction wrapped_res_norm(storage_size, residual_norm);
  WrapperScalarFunctionGradient wrapped_res_norm_grad(storage_size, residual_norm_grad);
  WrapperVectorFunction wrapped_in_equ(storage_size, storage_size, in_equ);
  WrapperVectorFunctionGradient wrapped_in_equ_grad(storage_size, storage_size, in_equ_grad);
  // optimization
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
  // DataMatrix trainSamples;
  // createSamples(trainSamples);
  // learnerSGDE.initialize(trainSamples);
  // std::shared_ptr<Grid> grid(&learnerSGDE.getGrid());
  // DataVector alpha(learnerSGDE.getSurpluses());

  return 0;
}
