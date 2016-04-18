// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGRESSIONLEARNER_H
#define REGRESSIONLEARNER_H

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include "sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp"
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class RegressionLearner {
 public:
  RegressionLearner(sgpp::base::RegularGridConfiguration gridConfig,
                    sgpp::base::AdpativityConfiguration adaptivityConfig,
                    sgpp::solver::SLESolverConfiguration solverConfig,
                    sgpp::datadriven::RegularizationConfiguration regularizationConfig);
  void train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes);
  sgpp::base::DataVector predict(sgpp::base::DataMatrix& data);
  double getMSE(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& y);

 private:
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;

  /// sparse grid object
  std::unique_ptr<sgpp::base::Grid> grid;
  /// the grid's coefficients
  sgpp::base::DataVector weights;
  std::unique_ptr<sgpp::base::OperationMatrix> opMatrix;

  void initializeGrid(sgpp::base::RegularGridConfiguration GridConfig);
  std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> createDMSystem(
      sgpp::base::DataMatrix& trainDataset);
  std::unique_ptr<sgpp::solver::SLESolver> createSolver();

  void trainStep(size_t curStep, sgpp::datadriven::DMSystemMatrixBase& DMSystem,
                 sgpp::solver::SLESolver& solver, sgpp::base::DataVector& classes);
  double getMSE(const sgpp::base::DataVector& y, sgpp::base::DataVector yPrediction);
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // REGRESSIONLEARNER_H
