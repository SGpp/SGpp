// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CLASSIFICATIONLEARNER_HPP
#define CLASSIFICATIONLEARNER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>
#include <vector>
#include <set>
#include <utility>

namespace sgpp {
namespace datadriven {

class ClassificationLearner {
 public:
  ClassificationLearner(sgpp::base::RegularGridConfiguration gridConfig,
                        sgpp::base::AdpativityConfiguration adaptivityConfig,
                        sgpp::solver::SLESolverConfiguration solverConfig,
                        sgpp::datadriven::RegularizationConfiguration regularizationConfig);

  void train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes);
  sgpp::base::DataVector predict(sgpp::base::DataMatrix& data);
  size_t getGridSize() const;
  double getAccuracy(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& y);

 private:
  typedef double class_t;
  typedef std::pair<class_t, RegressionLearner> learner_t;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  RegularizationConfiguration regularizationConfig;

  std::vector<learner_t> learners;
  std::set<class_t> uniqueClasses;

  sgpp::base::DataVector generateYOneVsAll(const sgpp::base::DataVector& oldY, class_t classValue);
  double getAccuracy(const sgpp::base::DataVector& y, const sgpp::base::DataVector yPrediction);
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // CLASSIFICATIONLEARNER_HPP
