// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGD_HPP
#define LEARNERSGD_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class LearnerSGD : public sgpp::datadriven::Learner {
 public:
  LearnerSGD(sgpp::datadriven::RegularizationType& regularization, const bool isRegression,
             const bool isVerbose = true);

  /*
   * Implements stochastic gradient descent.
   *
   * @param trainDataset training dataset: x values
   * @param classes training dataset: y values
   * @param maxIterations stops after maxIterations
   * @param eps stop if alpha_i < eps for all i
   * @param lambda regularization factor
   * @param gamma step width
   * */
  virtual void train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes,
                     sgpp::base::RegularGridConfiguration& GridConfig, size_t maxIterations,
                     double eps, double lambda, double gamma);

  virtual ~LearnerSGD();

  sgpp::base::DataVector& getAlpha();
  sgpp::base::Grid& getGrid();

 private:
  int getRandom(int limit);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGD_HPP */
