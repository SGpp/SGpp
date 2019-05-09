/* Copyright (C) 2008-today The SG++ project
 *
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Scorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun, Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Base class for supervised learning used to fit a model and quantify accuracy using a
 * #sgpp::datadriven::Metric with either testing or cross validation. Splits a dataset into testing
 * and training parts, trains the model and measures the accuracy.
 */
class Scorer {
 public:
  /**
   * Constructor
   *
   * @param metric  #sgpp::datadriven::Metric to to quantify approximation quality of a trained
   * model. Scorer will take ownership of this object.
   * puts a random seed.
   */
  explicit Scorer(Metric* metric);

  /**
   * Move constructor
   * @param rhs R-value reference to a scorer object to moved from.
   */
  Scorer(Scorer&& rhs) = default;

  /**
   * Copy assign operator
   * @param rhs const reference to the scorer object to copy from.
   * @return rerefernce to this with updated values.
   */
  Scorer& operator=(const Scorer& rhs) = default;

  /**
   * Move assign operator
   * @param rhs R-value reference to an a scorer object to move from.
   * @return rerefernce to this with updated values.
   */
  Scorer& operator=(Scorer&& rhs) = default;

  /**
   * Destructor.
   */
  ~Scorer() = default;

  /**
   * evaluate the accuracy on the test set using the #sgpp::datadriven::Metric.
   *
   * @param model model to be fitted based on the train dataset.
   * @param testDataset dataset used quantify accuracy using #sgpp::datadriven::Metric.
   * @return accuracy of the fit.
   */
  double test(ModelFittingBase& model, Dataset& testDataset);

 private:
  /**
   * evaluate the accuracy on the test set using the #sgpp::datadriven::Metric.
   * Uses ScaLAPACK to distribute evaluation (depends on model)
   *
   * @param model model to be fitted based on the train dataset.
   * @param testDataset dataset used quantify accuracy using #sgpp::datadriven::Metric.
   * @return accuracy of the fit.
   */
  double testDistributed(ModelFittingBase& model, Dataset& testDataset);

  /**
   * #sgpp::datadriven::Metric to be used to quantify accuracy of the fit.
   */
  std::unique_ptr<Metric> metric;
};

} /* namespace datadriven */
} /* namespace sgpp */
