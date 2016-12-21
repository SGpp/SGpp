/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidation.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Supervised learning with cross validation used to fit a model and quantify accuracy using a
 * #sgpp::datadriven::Metric.
 *
 * Splits a dataset into testing and training parts, trains the model and measures average accuracy
 * and standard deviation of the fits.
 */
class CrossValidation : public Scorer {
 public:
  /**
   * Constructor
   *
   * @param metric  #sgpp::datadriven::Metric to to quantify approximation quality of a trained
   * model. Scorer will take ownership of this object.
   * @param shuffling #sgpp::datadriven::ShufflingFunctor to rearrange samples of a dataset in the
   * desired manner, ready to be split into testing and training sets. Scorer will take ownership of
   * this object.
   * @param seed seed for randomization in #sgpp::datadriven::ShufflingFunctor. Default is -1 which
   * puts a random seed.
   * @param foldNumber amount of folds used for cross validation.
   */
  CrossValidation(Metric* metric, ShufflingFunctor* shuffling, int64_t seed = -1,
                  size_t foldNumber = 5);

  Scorer* clone() const override;

  /**
   * Train and test a model on a dataset and provide a score to quantify the approximation quality.
   * If multiple models are trained, calculate the standard deviation between the different fits.
   * @param model A model to be fitted on the training part of the dataset.
   * @param dataset Set of samples to use for fitting and testing the model.
   * @param stdDeviation return standard deviation between different runs.
   * @return average accuracy of all fits as calculated by the #metric provided.
   */
  virtual double calculateScore(ModelFittingBase& model, Dataset& dataset,
                                double* stdDeviation = nullptr);

 private:
  /**
   * amount of folds used for cross validation.
   */
  size_t foldNumber;
};

} /* namespace datadriven */
} /* namespace sgpp */
