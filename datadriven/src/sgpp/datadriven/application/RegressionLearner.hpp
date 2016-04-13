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
  void train();

 private:
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;

  /// sparse grid object
  std::unique_ptr<sgpp::base::Grid> grid;
  /// the grid's coefficients
  sgpp::base::DataVector weights;
  sgpp::base::DataMatrix samples;

  void initializeGrid(sgpp::base::RegularGridConfiguration GridConfig);
  std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> createDMSystem(
      sgpp::base::DataMatrix& trainDataset, double lambda);
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // REGRESSIONLEARNER_H
