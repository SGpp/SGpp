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


namespace SGPP {

namespace datadriven {

class LearnerSGD: public SGPP::datadriven::Learner {

 public:
  LearnerSGD(SGPP::datadriven::RegularizationType& regularization,
             const bool isRegression, const bool isVerbose = true);

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
  virtual void train(
    SGPP::base::DataMatrix& trainDataset,
    SGPP::base::DataVector& classes,
    SGPP::base::RegularGridConfiguration& GridConfig,
    size_t maxIterations,
    float_t eps,
    float_t lambda,
    float_t gamma
  );

  virtual ~LearnerSGD();

  SGPP::base::DataVector* getAlpha();
  SGPP::base::Grid* getGrid();

 private:
  int getRandom(int limit);

};
}
}

#endif /* LEARNERSGD_HPP */
