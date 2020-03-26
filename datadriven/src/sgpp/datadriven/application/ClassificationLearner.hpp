// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * @brief The ClassificationLearner class
 * Solves a classification problem.
 */
class ClassificationLearner {
 public:
  /**
   * @brief ClassificationLearner
   * @param gridConfig
   * @param adaptConfig
   * @param solverConfig is the solver used during each adaptivity step
   * @param finalSolverConfig is the solver used to build the final model
   * @param regularizationConfig
   * @param terms is a vector that contains all desired interaction terms.
   * For example, if we want to include grid points that model an
   * interaction between the first and the second predictor, we would
   * include the vector [1,2] in terms.
   */
  ClassificationLearner(sgpp::base::RegularGridConfiguration gridConfig,
                        sgpp::base::AdaptivityConfiguration adaptConfig,
                        sgpp::solver::SLESolverConfiguration solverConfig,
                        sgpp::solver::SLESolverConfiguration finalSolverConfig,
                        sgpp::datadriven::RegularizationConfiguration regularizationConfig,
                        std::set<std::set<size_t>> terms);
  /**
   * @brief ClassificationLearner
   * @param gridConfig
   * @param adaptConfig
   * @param solverConfig is the solver used during each adaptivity step
   * @param finalSolverConfig is the solver used to build the final model
   * @param regularizationConfig
   */
  ClassificationLearner(sgpp::base::RegularGridConfiguration gridConfig,
                        sgpp::base::AdaptivityConfiguration adaptConfig,
                        sgpp::solver::SLESolverConfiguration solverConfig,
                        sgpp::solver::SLESolverConfiguration finalSolverConfig,
                        sgpp::datadriven::RegularizationConfiguration regularizationConfig);
  /**
   * @brief train fits a sparse grid regression model.
   * @param trainDataset is the design matrix
   * @param classes is the (continuous) target
   */
  void train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes);
  /**
   * @brief predict
   * @param data are observations
   * @return the predicted target for matrix data
   */
  sgpp::base::DataVector predict(sgpp::base::DataMatrix& data);
  /**
   * @brief predict
   * @param data are observations
   * @return the predicted target for matrix data and the certainty of the model,
   * that the observation belongs to the predicted class.
   */
  std::pair<sgpp::base::DataVector, sgpp::base::DataVector> predictWithCertainty(
      sgpp::base::DataMatrix& data);
  /**
   * @brief getGridSize
   * @return the grid size
   */
  size_t getGridSize() const;
  /**
   * @brief getAccuracy
   * @param data is the design matrix
   * @param y is the target
   * @return the accuracy of the prediction of the model for the matrix data
   */
  double getAccuracy(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& y);

 private:
  typedef double class_t;
  typedef std::pair<class_t, RegressionLearner> learner_t;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdaptivityConfiguration adaptConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::solver::SLESolverConfiguration finalSolverConfig;
  RegularizationConfiguration regularizationConfig;
  std::set<std::set<size_t>> terms;

  std::vector<learner_t> learners;
  std::set<class_t> uniqueClasses;

  std::vector<std::pair<class_t, sgpp::base::DataVector>> getPredictions(
      sgpp::base::DataMatrix& data);
  sgpp::base::DataVector generateYOneVsAll(const sgpp::base::DataVector& oldY, class_t classValue);
  double getAccuracy(const sgpp::base::DataVector& y, const sgpp::base::DataVector yPrediction);
};

}  // namespace datadriven
}  // namespace sgpp
