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
#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>
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
   * @param shuffling #sgpp::datadriven::ShufflingFunctor to rearrange samples of a dataset in the
   * desired manner, ready to be split into testing and training sets. Scorer will take ownership of
   * this object.
   * @param seed seed for randomization in #sgpp::datadriven::ShufflingFunctor. Default is -1 which
   * puts a random seed.
   */
  Scorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed = -1);

  /**
   * Copy constructor
   * @param rhs const reference to the scorer object to copy from.
   */
  Scorer(const Scorer& rhs);

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
  Scorer& operator=(const Scorer& rhs);

  /**
   * Move assign operator
   * @param rhs R-value reference to an a scorer object to move from.
   * @return rerefernce to this with updated values.
   */
  Scorer& operator=(Scorer&& rhs) = default;

  /**
   * virtual destructor.
   */
  virtual ~Scorer() = default;

  /**
   * Polymorphic clone pattern
   * @return deep copy of this object. New object is owned by caller.
   */
  virtual Scorer* clone() const = 0;

  /**
   * Train and test a model on a dataset and provide a score to quantify the approximation quality.
   * If multiple models are trained, calculate the standard deviation between the different fits.
   * @param model A model to be fitted on the training part of the dataset.
   * @param dataset Set of samples to use for fitting and testing the model.
   * @param stdDeviation If multiple models are trained (e.g. for cross validation) calculate
   * standard deviation.
   * @return accuracy of the fit as calculated by the #metric provided.
   */
  virtual double calculateScore(ModelFittingBase& model, Dataset& dataset,
                                double* stdDeviation = nullptr) = 0;

 protected:
  /**
   * Helper method to generate an ordering for the samples of the dataset based on the shuffling
   * functor.
   * @param data: Dataset to be permuted.
   * @param randomizedIndices: vector with the same size as the dataset. Will be initialized with
   * contiguous values (0 -> vector.size()) and permuted in place
   */
  void randomizeIndices(const Dataset& data, std::vector<size_t>& randomizedIndices);

  /**
   * Split dataset into testing and training set.
   * @param fullDataset full dataset containing the samples to be split into testing and training
   * set.
   * @param trainDataset dataset where training samples will be stored. Needs to have the correct
   * size.
   * @param testDataset dataset where testing samples will be stored. Needs to have the correct
   * size.
   * @param randomizedIndices vector of permuted indices, describing in which order samples will be
   * read from the full dataset.
   * @param offset offset the testing set by the desired amount of samples. Used to generate testing
   * and training portions for cross validation. The samples skipped by the training set because of
   * the offset, will be assigned to the training set.
   */
  void splitSet(const Dataset& fullDataset, Dataset& trainDataset, Dataset& testDataset,
                const std::vector<size_t>& randomizedIndices, size_t offset = 0);

  /**
   * evaluate the accuracy on the test set using the #metric.
   *
   * @param model model to be fitted based on the train dataset.
   * @param testDataset dataset used quantify accuracy using #metric.
   * @return accuracy of the fit.
   */
  double test(ModelFittingBase& model, Dataset& testDataset);

  /**
   * Fit the model on the train dataset and evaluate the accuracy on the test set. Includes some
   * verbose output.
   *
   * @param model model to be fitted based on the train dataset
   * @param trainDataset dataset used for fitting the model.
   * @param testDataset dataset used quantify accuracy using #metric.
   * @return accuracy of the fit.
   */
  double train(ModelFittingBase& model, Dataset& trainDataset, Dataset& testDataset);

  /**
   * Fit the model on the train dataset and evaluate the accuracy on the test set. Includes some
   * verbose output.
   *
   * @param model model which is refined based on train dataset.
   * @param testDataset dataset used quantify accuracy using #metric.
   * @return accuracy of the fit after refinement.
   */
  double refine(ModelFittingBase& model, Dataset& testDataset);

  /**
   * #sgpp::datadriven::Metric to be used to quantify accuracy of the fit.
   */
  std::unique_ptr<Metric> metric;

  /**
   * #sgpp::datadriven::ShufflingFunctor used to rearrange samples of a dataset in the
   * desired manner, ready to be split into testing and training sets
   */
  std::unique_ptr<ShufflingFunctor> shuffling;
};

} /* namespace datadriven */
} /* namespace sgpp */
