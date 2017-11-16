// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/globaldef.hpp>

#include <random>
#include <functional>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridType;
using sgpp::base::GridPoint;
using sgpp::base::OperationEval;
using sgpp::base::OperationMatrix;
using sgpp::base::OperationMultipleEval;
using sgpp::datadriven::DensitySystemMatrix;

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
  std::function<double(DataVector&)> residual_norm =
    [numSamples, AlambC_ptr, b_ptr](DataVector& alpha) -> double {
    DataVector resultVec(numSamples);
    AlambC_ptr->mult(alpha, resultVec);
    resultVec.sub(*b_ptr);
    return resultVec.dotProduct(resultVec);
  };

  // define gradient for residual_norm
  // = 2*((A+lambda*C)*alpha - b).T * (A+lambdaC)
  // std::function<double(DataVector&)> gradient_residual_norm =
    // [numSamples, AlambC_ptr, b_ptr](DataVector& alpha) -> double {
    // DataVector resultVec(numSamples);
    // AlambC_ptr->mult(alpha, resultVec);
    // resultVec -= *b_ptr;
    // return resultVec.dotProduct(resultVec);
  // }

  // define inequality constraint
  // for alpha: eval at all grid points >= 0
  // i.e. PSI*alpha >= 0
  std::function<DataVector(const DataVector&)> inEqu =
    [op_eval, numSamples, gridStorage, dims](const DataVector& alpha) -> DataVector {
    DataVector resultVec(gridStorage->getSize());
    GridPoint gp;
    for (size_t i = 0; i < gridStorage->getSize(); i++) {
      gp = gridStorage->getPoint(i);
      DataVector coords(dims);
      gp.getStandardCoordinates(coords);
      resultVec[i] = op_eval->eval(alpha, coords);
    }
    return resultVec;
  };

  // define inequality constraint gradient
  // = PSI
  std::function<void(const DataVector&, DataVector&, DataMatrix&)> inEquGrad =
    [op_eval, gridStorage, dims](const DataVector& alpha,
                                             DataVector& result, DataMatrix& jacobi) -> void {
    GridPoint gp;
    for (size_t i = 0; i < gridStorage->getSize(); i++) {
      gp = gridStorage->getPoint(i);
      DataVector coords(dims);
      DataVector tmp_alpha(gridStorage->getSize());
      gp.getStandardCoordinates(coords);
      result[i] = op_eval->eval(alpha, coords);
      // dumm gelöst
      for (size_t j = 0; j < gridStorage->getSize(); i++) {
        tmp_alpha.setAll(0.0);
        tmp_alpha.set(j, 1.0);
        double val = op_eval->eval(alpha, coords);
        jacobi.set(i, j, val);
      }
    }
  };

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
